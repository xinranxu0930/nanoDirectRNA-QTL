import pandas as pd
import pysam
import scipy
import argparse
import pickle
from collections import Counter
import os
import statsmodels.api as sm
import numpy as np
import ast

def get_reverse_complementary_sequence(seq):
    seqreverse = seq[::-1]
    transtable = str.maketrans('ATGCatgc','TACGtacg')
    finalseq = seqreverse.translate(transtable).upper()
    return finalseq

def count_haplotype(chrom, end, strand, bamfile, read_st_dict, snp_file_dict,threadsnum,base_minQ,read_minQ):
    # read_st_dict read对应的稳定性指标 readid:stability
    # snp_file_dict snp位点信息 # chrom_pos0: rsid;A1;A2;EAF
    snp_bases = {}
    # snp_bases 存放的是snp位点 [pos0]:{read1:"A",read2:"A",read3:"T"}
    with pysam.AlignmentFile(bamfile, "rb",threads=threadsnum) as samfile:
        for pileupcolumn in samfile.pileup(chrom, 0, end,min_base_quality=base_minQ,min_mapping_quality=read_minQ,stepper="samtools",max_depth=50000):
            base_pos = str(pileupcolumn.reference_pos)
            ## 处理snp位点的情况
            # 1判断位置是否是snp
            if f"{chrom}_{base_pos}" in snp_file_dict:
                # 2判断是否是杂合子位点
                seqs = [i.upper() for i in pileupcolumn.get_query_sequences()]  # 被覆盖位点处的所有碱基
                if len(set(seqs)) > 1:
                    snp_bases[base_pos] = {
                        i: j for i, j in zip(pileupcolumn.get_query_names(), seqs)
                    }
    if len(snp_bases) == 0:
        print(f'{chrom} {strand}中,没有snp被覆盖')
        return None
    haplotype_df = pd.DataFrame(columns=[
        "chrom","strand","snp_pos_1base",
        "rsID","A1","A2",
        "A1_STscore_l","A2_STscore_l",
        "EAF",
        ])
    result = start_get_haplotypes(chrom, strand, snp_bases, read_st_dict, snp_file_dict)
    if result is None:
        print(f'{chrom} {strand}中没有符合条件的snp')
        return None
    else:
        ith = 0
        for i in result:
            haplotype_df.loc[ith] = (
                chrom,strand,int(i[5])+1,
                i[4],i[0],i[1],
                i[2],i[3],
                i[6]
            )
            ith += 1
        return haplotype_df

def start_get_haplotypes(chrom,strand,snp_bases,read_st_dict,snp_file_dict):
    res_l = [] # 每一个元素都表示一个snp的结果 [a1,a2,A1_STscore_l, A2_STscore_l,snpID,snp_pos,eaf]
    # snp_bases [pos0]:{read1:"A",read2:"A",read3:"T"}
    # read_st_dict readid:stability
    for snp_pos, snp_base in snp_bases.items():
        readid_base_STscore = [] # 每一个新的snp都有一个新的readid_base_STscore，里面是覆盖这个snp的所有read的base和STscore
        for read,base in snp_base.items():
            if read not in read_st_dict:
                continue
            else:
                readid_base_STscore.append([base,read_st_dict[read]])
        res = get_haplotypes(readid_base_STscore, snp_file_dict[f"{chrom}_{snp_pos}"],strand)
        if res is not None:
            a1,a2,A1_STscore_l, A2_STscore_l,snpID,eaf = res
            res_l.append([a1,a2,A1_STscore_l, A2_STscore_l,snpID,snp_pos,eaf])
    if len(res_l) != 0:
        return res_l
    else:
        return None

def get_haplotypes(readid_base_STscore,snp_file_dict_res,strand):
    if len(readid_base_STscore) == 0:
        return None
    snpID = snp_file_dict_res.split(";")[0]
    eaf = float(snp_file_dict_res.split(";")[3])
    all_base_list = list(i[0] for i in readid_base_STscore)
    all_ST_list = list(i[1] for i in readid_base_STscore)
    all_base = list(set(all_base_list))
    A1,A2 = snp_file_dict_res.split(";")[1],snp_file_dict_res.split(";")[2]
    a1,a2 = "",""
    if strand == "-":
        A1 = get_reverse_complementary_sequence(A1)
        A2 = get_reverse_complementary_sequence(A2)
    for base in all_base:
        if base == A1:
            a1 = base
        elif base == A2:
            a2 = base
    if (a1 == "") and (a2 == ""):
        return None
    # 遍历all_ST_type列表
    A1_STscore_l, A2_STscore_l = [], []
    for i in range(len(all_base_list)):
        if all_base_list[i] == a1:
            A1_STscore_l.append(all_ST_list[i])
        elif all_base_list[i] == a2:
            A2_STscore_l.append(all_ST_list[i])
    return a1,a2,A1_STscore_l, A2_STscore_l,snpID,eaf

# 线性回归
def analyze_snp_stability(A1_stability, A2_stability,min_cov):
    if (len(A1_stability) >= min_cov) & (len(A2_stability) >= min_cov):
        combined_stability = A1_stability + A2_stability
        snp_labels = [0] * len(A1_stability) + [1] * len(A2_stability)  # 0 表示 A1，1 表示 A2
        df = pd.DataFrame({
            'stability': combined_stability,
            'snp': snp_labels
        })
        # 添加常数项（截距）
        X = sm.add_constant(df['snp'])
        y = df['stability']
        # 进行线性回归
        model = sm.OLS(y, X).fit()
        beta = model.params[1]  # 获取 SNP 的 beta 值
        p_value = model.pvalues[1]  # 获取 beta 值的 p 值
        return p_value,beta
    else:
        return None, None

def apply_linear_regression(row,min_cov):
    return analyze_snp_stability(row['A1_STscore_l'], row['A2_STscore_l'],min_cov)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Nanopore dapaect RNA data call RNA stability QTL(stqtl).')
    parser.add_argument("-b","--bam", type=str, help="bam file path")
    parser.add_argument("--snp_info", type=str, help="processed snp pkl file path")
    parser.add_argument("-p","--outdirpre", type=str, help="outdir and pre")
    parser.add_argument("--readST", type=str, help="read stability csv file")
    parser.add_argument("-c","--chrom", type=str, help="chromosome")
    parser.add_argument("-s","--strand", type=str, help="different strand processing")
    parser.add_argument("--geno_size", type=str, help="genome size file path")
    parser.add_argument("-t","--threads", type=int, default=4, help="threads number (default: 4)")
    parser.add_argument("--base_minQ", type=int, default=5, help="base min qscore(default=5)")
    parser.add_argument("--read_minQ", type=int, default=0, help="read min qscore(default=0)")
    parser.add_argument("--snp_min_cov", type=int, default=5, help="ref or alt min coverage(default=5)")
    args = parser.parse_args()

    output_path = f"{args.outdirpre}_haplotype_{args.chrom}_{args.strand}_tmp.csv"

    base_minQ = args.base_minQ-1 if args.base_minQ != 0 else 0
    read_minQ = args.read_minQ-1 if args.read_minQ != 0 else 0

    # csv2dict
    read_st = pd.read_csv(args.readST, usecols = [0,5])
    read_st_dict = read_st.set_index('readID')['stability_zscore'].to_dict() # readid:stability
    # txt2dict
    snp_info = pd.read_csv(args.snp_info, sep="\t")
    snp_info.columns = ["chrom","pos1","rsID","A1","A2","EAF"]
    snp_info["pos0"] = snp_info["pos1"].astype(int) - 1
    snp_info["k"] = snp_info["chrom"] + "_" + snp_info["pos0"].astype(str)
    snp_info["v"] = snp_info["rsID"] + ";" + snp_info["A1"] + ";" + snp_info["A2"] + ";" + snp_info["EAF"].astype(str)
    snp_dict = dict(zip(snp_info["k"], snp_info["v"])) # chrom_pos0: rsid;A1;A2;EAF

    geno_size_df = pd.read_csv(args.geno_size, sep="\t", header=None, names=["chrom","size"])
    geno_size_dict = dict(zip(geno_size_df["chrom"], geno_size_df["size"]))
    haplotype_df = pd.DataFrame(columns=[
        "chrom","strand","snp_pos_1base",
        "rsID","A1","A2",
        "A1_STscore_l","A2_STscore_l",
        "EAF",
        ])
    end = geno_size_dict[args.chrom]
    haplotype_df = pd.concat([haplotype_df, count_haplotype(args.chrom, end, args.strand, args.bam, read_st_dict, snp_dict, args.threads,base_minQ,read_minQ)],ignore_index=True)
    haplotype_df = haplotype_df[(haplotype_df['A1'] != '') & (haplotype_df['A2'] != '')] # df['A1']和df['A2']只要有一个是None 就删除这一行
    haplotype_df = haplotype_df.reset_index(drop=True)
    if len(haplotype_df) != 0:
        results = haplotype_df.apply(lambda row: apply_linear_regression(row, args.snp_min_cov), axis=1, result_type='expand')
        haplotype_df[['p_value', 'beta']] = results
        haplotype_df = haplotype_df[haplotype_df['p_value'].notna()]
        haplotype_df.to_csv(output_path, index=None)
        print(f"{output_path}已保存")