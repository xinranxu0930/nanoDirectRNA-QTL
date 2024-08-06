import pandas as pd
import pysam
import scipy
import argparse
import pickle
import os
import numpy as np


def get_reverse_complementary_sequence(seq):
    seqreverse = seq[::-1]
    transtable = str.maketrans('ATGCatgc','TACGtacg')
    finalseq = seqreverse.translate(transtable).upper()
    return finalseq

def count_haplotype(chrom, end, strand, bamfile, snp_file_dict, threadsnum, base_minQ, read_minQ):
    # snp_file_dict snp位点信息 chrom_pos0: rsid;A1;A2;MAF
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
        "A1_num","A2_num",
        "MAF",
        ])
    result = start_get_haplotypes(chrom, strand, snp_bases, snp_file_dict)
    if result is None:
        print(f'{chrom} {strand}中没有符合条件的snp')
        return None
    else:
        ith = 0
        for i in result:
            haplotype_df.loc[ith] = (
                chrom,strand,int(i[6])+1,
                i[4],i[0],i[1],
                i[2],i[3],
                i[5]
            )
            ith += 1
        return haplotype_df

def start_get_haplotypes(chrom,strand,snp_bases,snp_file_dict):
    res_l = [] # 每一个元素都表示一个snp的结果 [a1,a2,A1_num,A2_num,snpID,maf,snp_pos]
    # snp_bases [pos0]:{read1:"A",read2:"A",read3:"T"}
    for snp_pos, snp_base in snp_bases.items():
        snp_base_count = {} # 统计每个snp处每个碱基的数量 {A:10,T:20}
        for k, v in snp_base.items():
            if v in snp_base_count:
                snp_base_count[v] += 1
            else:
                snp_base_count[v] = 1
        if len(snp_base_count) < 2:
            continue
        res = get_haplotypes(snp_base_count, snp_file_dict[f"{chrom}_{snp_pos}"],strand)
        if res is not None:
            a1,a2,A1_num,A2_num,snpID,maf = res
            res_l.append([a1,a2,A1_num,A2_num,snpID,maf,snp_pos])
    if len(res_l) != 0:
        return res_l
    else:
        return None

def get_haplotypes(snp_base_count,snp_file_dict_res,strand):
    snpID = snp_file_dict_res.split(";")[0]
    all_base = list(set(snp_base_count.keys()))
    A1,A2 = snp_file_dict_res.split(";")[1],snp_file_dict_res.split(";")[2]
    maf = float(snp_file_dict_res.split(";")[3])
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
    if a1 != "" and a2 != "":
        A1_num,A2_num = int(snp_base_count[a1]),int(snp_base_count[a2])
    elif a2 == "":
        A1_num,A2_num  = int(snp_base_count[a1]),0
    else:
        A1_num,A2_num  = 0,int(snp_base_count[a2])
    return a1,a2,A1_num,A2_num,snpID,maf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Nanopore direct RNA data call ASE expression qtl.')
    parser.add_argument("-b","--bam", type=str, help="bam file path")
    parser.add_argument("--snp_info", type=str, help="processed snp pkl file path")
    parser.add_argument("-p","--outdirpre", type=str, help="outdir and pre")
    parser.add_argument("-c","--chrom", type=str, help="chromosome")
    parser.add_argument("-s","--strand", type=str, help="different strand processing")
    parser.add_argument("--geno_size", type=str, help="genome size file path")
    parser.add_argument("-t","--threads", type=int, default=4, help="threads number (default: 4)")
    parser.add_argument("--base_minQ", type=int, default=5, help="base min qscore(default=5)")
    parser.add_argument("--read_minQ", type=int, default=0, help="read min qscore(default=0)")
    args = parser.parse_args()

    base_minQ = args.base_minQ-1 if args.base_minQ != 0 else 0
    read_minQ = args.read_minQ-1 if args.read_minQ != 0 else 0

    output_path = f"{args.outdirpre}_haplotype_{args.chrom}_{args.strand}_tmp.csv"

    snp_info = pd.read_csv(args.snp_info, sep="\t")
    snp_info.columns = ["chrom","pos1","rsID","A1","A2","MAF"]
    snp_info["pos0"] = snp_info["pos1"].astype(int) - 1
    snp_info["k"] = snp_info["chrom"] + "_" + snp_info["pos0"].astype(str)
    snp_info["v"] = snp_info["rsID"] + ";" + snp_info["A1"] + ";" + snp_info["A2"] + ";" + snp_info["MAF"].astype(str)
    snp_dict = dict(zip(snp_info["k"], snp_info["v"])) # chrom_pos0: rsid;A1;A2;MAF

    geno_size_df = pd.read_csv(args.geno_size, sep="\t", header=None, names=["chrom", "size"])
    end = geno_size_df.loc[geno_size_df['chrom'] == args.chrom, 'size'].iloc[0]
    haplotype_df = pd.DataFrame(columns=[
        "chrom","strand","snp_pos_1base",
        "rsID","A1","A2",
        "A1_num","A2_num",
        "MAF",
        ])

    haplotype_df = pd.concat([haplotype_df, count_haplotype(args.chrom, end, args.strand, args.bam, snp_dict, args.threads,base_minQ,read_minQ)],ignore_index=True)
    if len(haplotype_df) != 0:
        haplotype_df = haplotype_df[(haplotype_df['A1'] != '') & (haplotype_df['A2'] != '')] # df['A1']和df['A2']只要有一个是None 就删除这一行
        df = haplotype_df.reset_index(drop=True)

        df.to_csv(output_path, index=None)
        print(f"{output_path}已保存")