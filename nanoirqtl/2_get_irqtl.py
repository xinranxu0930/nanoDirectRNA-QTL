import pandas as pd
import pysam
import scipy
import argparse
import pickle
from collections import Counter
import os
from scipy.stats import fisher_exact
import numpy as np
import ast

def get_reverse_complementary_sequence(seq):
    seqreverse = seq[::-1]
    transtable = str.maketrans('ATGCatgc','TACGtacg')
    finalseq = seqreverse.translate(transtable).upper()
    return finalseq

def count_haplotype(chrom, end, strand, bamfile, read_isoform_dict, snp_file_dict,threadsnum,base_minQ,read_minQ):
    # read_isoform_dict read对应的isoform read:isoform
    # snp_file_dict snp位点信息 chrom_pos0: rsid;A1;A2;EAF
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
    haplotype_df = pd.DataFrame(
        columns=[
        "chrom","strand","snp_pos_1base",
        "rsID","A1","A2",
        "A1_isoform_l","A2_isoform_l",
        "isoform_id_l",
        "EAF",
        ])
    result = start_get_haplotypes(chrom, strand, snp_bases, read_isoform_dict, snp_file_dict)
    if result is None:
        print(f'{chrom} {strand}中没有符合条件的snp')
        return None
    else:
        ith = 0
        for i in result:
            haplotype_df.loc[ith] = (
                chrom,strand,int(i[7])+1,
                i[5],i[0],i[1],
                i[2],i[3],
                i[4],
                i[6]
            )
            ith += 1
        return haplotype_df

def start_get_haplotypes(chrom,strand,snp_bases,read_isoform_dict,snp_file_dict):
    res_l = [] # 每一个元素都表示一个snp的结果 [snp_genome_0base,A1,A2,A1_l,A2_l,isoform_l]
    # snp_bases [pos0]:{read1:"A",read2:"A",read3:"T"}
    # read_isoform_dict read:isoform
    for snp_pos, snp_base in snp_bases.items():
        readid_base_isoform = [] # 每一个新的snp都有一个新的readid_base_isoform，里面是覆盖这个snp的所有read的base和isoform
        for read,base in snp_base.items():
            if read not in read_isoform_dict:
                continue
            else:
                readid_base_isoform.append([base,read_isoform_dict[read]]) # 这里是某个snp处的read、base和isoform的对应情况，因为不需要知道具体的read是什么，只要知道base和isoform的对应关系就可以，所以放在list中
        res = get_haplotypes(readid_base_isoform, snp_file_dict[f"{chrom}_{snp_pos}"],strand)
        if res is not None:
            a1,a2,A1_count_l,A2_count_l,all_isoform_type,snpID,eaf = res
            res_l.append([a1,a2,A1_count_l,A2_count_l,all_isoform_type,snpID,eaf,snp_pos])
    if len(res_l) != 0:
        return res_l
    else:
        return None

def get_haplotypes(readid_base_isoform,snp_file_dict_res,strand):
    if len(readid_base_isoform) == 0:
        return None
    snpID = snp_file_dict_res.split(";")[0]
    eaf = float(snp_file_dict_res.split(";")[3])
    all_base_list = list(i[0] for i in readid_base_isoform)
    all_isoform_type = list(set(i[1] for i in readid_base_isoform))
    A1,A2 = snp_file_dict_res.split(";")[1],snp_file_dict_res.split(";")[2]
    a1,a2 = "",""
    if strand == "-":
        A1 = get_reverse_complementary_sequence(A1)
        A2 = get_reverse_complementary_sequence(A2)
    for base in all_base_list:
        if base == A1:
            a1 = base
        elif base == A2:
            a2 = base
    if (a1 == "") and (a2 == ""):
        return None
    # 遍历all_isoform_type列表
    A1_count_l, A2_count_l = [], []
    if a2 != "" and a1 != "":
        A1_count_l = [sum(1 for base_isoform in readid_base_isoform if base_isoform[0] == a1 and base_isoform[1] == isoform) for isoform in all_isoform_type]
        A2_count_l = [sum(1 for base_isoform in readid_base_isoform if base_isoform[0] == a2 and base_isoform[1] == isoform) for isoform in all_isoform_type]
    elif a2 == "":
        A1_count_l = [sum(1 for base_isoform in readid_base_isoform if base_isoform[0] == a1 and base_isoform[1] == isoform) for isoform in all_isoform_type]
        A2_count_l = [0 for _ in range(len(all_isoform_type))]
    else:
        A2_count_l = [sum(1 for base_isoform in readid_base_isoform if base_isoform[0] == a2 and base_isoform[1] == isoform) for isoform in all_isoform_type]
        A1_count_l = [0 for _ in range(len(all_isoform_type))]
    return a1,a2,A1_count_l,A2_count_l,all_isoform_type,snpID,eaf

# fisher
def simulate_fisher_exact_se(A1_isoform1, A2_isoform1, A1_isoform2, A2_isoform2, min_cov):
    if (A1_isoform1+A1_isoform2 >= min_cov) & (A1_isoform1+A1_isoform2 >= min_cov):
        observed = np.array([[A1_isoform1, A2_isoform1], [A1_isoform2, A2_isoform2]])
        if 0 in [A1_isoform1, A2_isoform1, A1_isoform2, A2_isoform2]:
            odds_ratio = ((A1_isoform1 + 0.75) * (A2_isoform2 + 0.75)) / ((A2_isoform1 + 0.75) * (A1_isoform2 + 0.75))
            p_value = fisher_exact(observed)[1]
        else:
            odds_ratio, p_value = fisher_exact(observed)
        if np.isinf(odds_ratio) or odds_ratio == 0 :
            return p_value, None
        beta = np.log(odds_ratio)
        return p_value, beta
    else:
        return None, None

def apply_simulate_fisher(row, min_cov):
    return simulate_fisher_exact_se(row['A1_isoform1'], row['A2_isoform1'], row['A1_isoform2'], row['A2_isoform2'],min_cov)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Nanopore direct RNA data call isoform rato qtl(irqtl).')
    parser.add_argument("-b","--bam", type=str, help="bam file path")
    parser.add_argument("--snp_info", type=str, help="processed snp pkl file path")
    parser.add_argument("-p","--outdirpre", type=str, help="outdir and pre")
    parser.add_argument("-r","--read_isoform_dict", type=str, help="read isoform pkl file path")
    parser.add_argument("-c","--chrom", type=str, help="chromosome")
    parser.add_argument("-s","--strand", type=str, help="different strand processing")
    parser.add_argument("--geno_size", type=str, help="genome size file path")
    parser.add_argument("-t","--threads", type=int, default=4, help="threads number (default: 4)")
    parser.add_argument("--base_minQ", type=int, default=5, help="base min qscore(default=5)")
    parser.add_argument("--read_minQ", type=int, default=0, help="read min qscore(default=0)")
    parser.add_argument("--snp_min_cov", type=int, default=5, help="ref or alt min coverage(default=5)")
    args = parser.parse_args()

    output_path = f"{args.outdirpre}_haplotype_{args.chrom}_{args.strand}_tmp.csv"

    base_minQ = args.base_minQ-1
    read_minQ = args.read_minQ-1

    with open(args.read_isoform_dict, 'rb') as file:
        read_isoform_dict = pickle.load(file) # readID: isoformID

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
        "A1_isoform_l","A2_isoform_l",
        "isoform_id_l",
        "EAF",
        ])
    end = geno_size_dict[args.chrom]
    haplotype_df = pd.concat([haplotype_df, count_haplotype(args.chrom, end, args.strand, args.bam, read_isoform_dict, snp_dict, args.threads, base_minQ, read_minQ)],ignore_index=True)
    # haplotype_df = haplotype_df[haplotype_df['A1'].notnull() & haplotype_df['A2'].notnull()]
    if len(haplotype_df) != 0:
        haplotype_df = haplotype_df[haplotype_df['isoform_id_l'].str.len() != 1] # 删除df['isoform_id_l']中长度为1的行
        haplotype_df = haplotype_df[(haplotype_df['A1'] != '') & (haplotype_df['A2'] != '')] # df['A1']和df['A2']只要有一个是None 就删除这一行
        df = haplotype_df.reset_index(drop=True)
        for i in df.index:
            if len(df.loc[i, "isoform_id_l"])>2:
                subset_df = pd.DataFrame([df.loc[i, "isoform_id_l"], df.loc[i, "A1_isoform_l"], df.loc[i, "A2_isoform_l"]])
                subset_df.columns = subset_df.iloc[0]
                subset_df = subset_df.iloc[1:].reset_index(drop=True)
                subset_df.index=['A1', 'A2']
                max_sum_column = subset_df.sum().astype(int).idxmax() # 找出列的总和最大的那个列
                new_data = {
                    max_sum_column: subset_df[max_sum_column],
                    'other': subset_df.drop(columns=[max_sum_column]).sum(axis=1)
                } # 保留最大列并将其他列合并
                grouped_df = pd.DataFrame(new_data)
                df.loc[i, "A1_isoform1"] = grouped_df.iloc[0][max_sum_column]
                df.loc[i, "A2_isoform1"] = grouped_df.iloc[1][max_sum_column]
                df.loc[i, "A1_isoform2"] = grouped_df.iloc[0]["other"]
                df.loc[i, "A2_isoform2"] = grouped_df.iloc[1]["other"]
                df.loc[i, "isoform_id_l"] = f'{max_sum_column},other'
            else:
                df.loc[i, "A1_isoform1"] = df.loc[i, "A1_isoform_l"][0]
                df.loc[i, "A2_isoform1"] = df.loc[i, "A2_isoform_l"][0]
                df.loc[i, "A1_isoform2"] = df.loc[i, "A1_isoform_l"][1]
                df.loc[i, "A2_isoform2"] = df.loc[i, "A2_isoform_l"][1]
                df.loc[i, "isoform_id_l"] = ','.join(df.loc[i, "isoform_id_l"])
        del df['A1_isoform_l'], df['A2_isoform_l']
        results = df.apply(lambda row: apply_simulate_fisher(row, args.snp_min_cov), axis=1, result_type='expand')
        df[['p_value','beta']] = results
        df = df[df['p_value'].notna()]
        df.to_csv(output_path, index=None)
        print(f"{output_path}已保存")