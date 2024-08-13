import pandas as pd
import pysam
import scipy
import argparse
import pickle
from collections import Counter
import os
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import fisher_exact


def create_m6A_dict(row):
    if row['methy_rate'] != 100:
        return f'{row["pos_0base"]}', row["methy_rate"]

def get_reverse_complementary_sequence(seq):
    seqreverse = seq[::-1]
    transtable = str.maketrans('ATGCatgc','TACGtacg')
    finalseq = seqreverse.translate(transtable).upper()
    return finalseq

## call m6Aqtl
def count_haplotype(chrom, start, end, strand, bamfile, methy_dict, m6A_id_dict, snp_file_dict,base_minQ,read_minQ):
    # methy_dict 甲基化位点信息 pos0:methy_rate
    # m6A_id_dict 含有甲基化的read chrom_pos0_strand: read1;read2;...
    # snp_file_dict snp位点信息 chrom_pos0: rsid;A1;A2;EAF
    m6A_bases,snp_bases,snp_m6A_bases = {},{},{}
    # m6A_bases 存放目标区域中甲基化位点处所有的read id，这里同时包括甲基化的read和非甲基化的read，和m6A_id_dict联合使用
    # 后两个 XXX_bases[pos0]:{read1:"A",read2:"A",read3:"T"}
    # snp_bases 存放的是snp位点
    # snp_m6A_bases 存放的位点同时是snp和甲基化位点
    with pysam.AlignmentFile(bamfile, "rb",threads=10) as samfile:
        for pileupcolumn in samfile.pileup(chrom, start, end,min_base_quality=base_minQ,min_mapping_quality=read_minQ,stepper="samtools",max_depth=50000):
            base_pos = str(pileupcolumn.reference_pos)
            ## 先处理甲基化位点的情况
            if base_pos in methy_dict:
                m6A_bases[base_pos] = pileupcolumn.get_query_names()
            ## 处理snp位点的情况
            # 1判断位置是否是snp
            if f"{chrom}_{base_pos}" in snp_file_dict:
                # 2判断是否是杂合子位点
                seqs = [i.upper() for i in pileupcolumn.get_query_sequences()]  # 被覆盖位点处的所有碱基
                if len(set(seqs)) > 1:
                    # 3判断snp是否是m6A位点
                    if base_pos in m6A_bases:
                        snp_m6A_bases[base_pos] = {
                            i: j for i, j in zip(pileupcolumn.get_query_names(), seqs)
                        }
                    else:
                        snp_bases[base_pos] = {
                            i: j for i, j in zip(pileupcolumn.get_query_names(), seqs)
                        }
    if len(m6A_bases) == 0 and (len(snp_m6A_bases)==0 or len(snp_bases) == 0):
        return None
    haplotype_df = pd.DataFrame(columns=[
        "chrom","strand","m6A_pos_1base","snp_pos_1base",
        "rsID","A1","A2",
        "A1_A","A2_A","A1_m6A","A2_m6A",
        "ism6A?","methy_rate",
        "EAF"
        ])
    ith = 0
    bases_list = [snp_bases, snp_m6A_bases]
    ism6A_list = ["", "Yes"]
    for bases, ism6A in zip(bases_list, ism6A_list):
        result = start_get_haplotypes(chrom, strand, m6A_bases, bases, m6A_id_dict, snp_file_dict)
        if result is None:
            continue
        else:
            for i in result:
                haplotype_df.loc[ith] = (
                    chrom,strand,int(i[0]) + 1,int(i[1]) + 1,
                    i[8],i[2],i[3],
                    i[4],i[5],i[6],i[7],
                    ism6A,methy_dict[str(i[0])],
                    i[9],
                )
                ith += 1
    if len(haplotype_df) != 0:
        return haplotype_df
    else:
        return None

def start_get_haplotypes(chrom,strand,m6A_bases,snp_bases,m6A_id_dict,snp_file_dict):
    res_l = [] # 每一个元素都表示一个位点的结果 [m6A_pos0,snp_0base,ref,alt,A1_A,A2_A,A1_m6A,A2_m6A]
    # snp_bases [pos0]:{read1:"A",read2:"A",read3:"T"}
    # m6A_bases [A1_pos]:[readid1,readid12,...] 位点处所有的read 甲基化 非甲基化都有
    # m6A_id_dict chrom_pos0_strand: read1;read2;... 所有甲基化的read
    for m6A_pos, m6A_allreadid_list in m6A_bases.items():
        for snp_pos, snp_base in snp_bases.items():
            m6A_readid = set(m6A_id_dict[f"{chrom}_{m6A_pos}_{strand}"].split(";")) # 某个点只是甲基化的所有read
            res = get_haplotypes(snp_base, m6A_allreadid_list, m6A_readid, snp_file_dict[f"{chrom}_{snp_pos}"],strand)
            if res is not None:
                a1,a2,a1_A,a2_A,a1_m6A,a2_m6A,snpID,eaf = res
                res_l.append([m6A_pos,snp_pos,a1,a2,a1_A,a2_A,a1_m6A,a2_m6A,snpID,eaf])
    if len(res_l) != 0:
        return res_l
    else:
        return None

# 现在是确定了snp-m6A，输入snp处的readid-base；m6A出所有readid；m6A处甲基化的readid；该snp的rsid ref alt；strand
def get_haplotypes(snp_base, m6A_allreadid_list, m6A_readid,snp_file_dict_res,strand):
    intersect_read_set = set(snp_base.keys()) & set(m6A_allreadid_list)  # 两个字典共有的read id的集合，待分析的所有read
    if len(intersect_read_set) == 0:
        return None
    in_m6A_readid_set = intersect_read_set & m6A_readid # 含有甲基化的read
    not_in_m6A_readid_set = intersect_read_set - m6A_readid # 含有正常A的read
    # most_2base = [item[0] for item in Counter(list(snp_base.values())).most_common(2)]
    all_base = list(set(snp_base.values()))
    snpID = snp_file_dict_res.split(";")[0]
    eaf = float(snp_file_dict_res.split(";")[3])
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
    a1_read = set(k for k, v in snp_base.items() if v == a1)
    a2_read = set(k for k, v in snp_base.items() if v == a2)
    a1_A = len(a1_read & not_in_m6A_readid_set)
    a2_A = len(a2_read & not_in_m6A_readid_set)
    a1_m6A = len(a1_read & in_m6A_readid_set)
    a2_m6A = len(a2_read & in_m6A_readid_set)
    return a1,a2,a1_A,a2_A,a1_m6A,a2_m6A,snpID,eaf

# fisher
def simulate_fisher_exact_se(A1_A, A2_A, A1_m6A, A2_m6A, num_simulations=100):
    observed = np.array([[A1_A, A2_A], [A1_m6A, A2_m6A]])
    if 0 in [A1_A, A2_A, A1_m6A, A2_m6A]:
        odds_ratio = ((A1_A + 0.75) * (A2_m6A + 0.75)) / ((A2_A + 0.75) * (A1_m6A + 0.75))
        p_value = fisher_exact(observed)[1]
    else:
        odds_ratio, p_value = fisher_exact(observed)
    n = int(A1_A+A2_A+A1_m6A+A2_m6A)
    if np.isinf(odds_ratio) or odds_ratio == 0 or n == 0:
        return p_value, n, None, None
    beta = np.log(odds_ratio)
    # Simulate to estimate the standard error
    odds_ratios = np.zeros(num_simulations)
    valid_simulations = 0
    for i in range(num_simulations):
        # Simulate contingency table based on margins
        n1 = A1_A + A2_A
        n2 = A1_m6A + A2_m6A
        # Ensure n2/n and (1 - n2/n) are valid probabilities
        p1 = np.clip(n2/n, 1e-10, 1-1e-10)  # Avoid probabilities exactly 0 or 1
        p2 = 1 - p1
        # Simulate cell counts under null hypothesis of independence
        simulated_table = np.array([[np.random.binomial(n1, p1), np.random.binomial(n1, p2)],
                                    [np.random.binomial(n - n1, p1), np.random.binomial(n - n1, p2)]])
        # Calculate odds ratio for simulated table
        sim_odds_ratio, _ = fisher_exact(simulated_table)
        if not np.isinf(sim_odds_ratio) and sim_odds_ratio != 0:
            odds_ratios[valid_simulations] = sim_odds_ratio
            valid_simulations += 1
    if valid_simulations > 2:
        odds_ratio_se = np.nanstd(odds_ratios[:valid_simulations], ddof=1)
        sd = odds_ratio_se * np.sqrt(n)
    else:
        sd = None
    return p_value, n, beta, sd

def simulate_fisher_exact(A1_A, A2_A, A1_m6A, A2_m6A,min_cov):
    p_value = None
    beta = None
    se_log_or = None
    if (A1_A+A1_m6A >= min_cov) & (A2_A+A2_m6A >= min_cov):
        A1_A += 0.75 if A1_A == 0 else 0
        A2_A += 0.75 if A2_A == 0 else 0
        A1_m6A += 0.75 if A1_m6A == 0 else 0
        A2_m6A += 0.75 if A2_m6A == 0 else 0
        observed = np.array([[A1_A, A2_A], [A1_m6A, A2_m6A]])
        odds_ratio = (A1_A * A2_m6A) / (A2_A * A1_m6A)
        se_log_or = np.sqrt(1/A1_A + 1/A2_m6A + 1/A2_A + 1/A1_m6A)
        p_value = fisher_exact(observed)[1]
        if np.isinf(odds_ratio) or odds_ratio == 0:
            beta = None
        else:
            beta = np.log(odds_ratio)
    return p_value, beta, se_log_or

def apply_simulate_fisher(row,min_cov):
    return simulate_fisher_exact(row['A1_A'], row['A2_A'], row['A1_m6A'], row['A2_m6A'],min_cov)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Nanopore direct RNA data call m6Aqtl.')
    parser.add_argument("-b","--bam", type=str, help="bam file path")
    parser.add_argument("--snp_info", type=str, help="processed snp pkl file path")
    parser.add_argument("-p","--outdirpre", type=str, help="outdir and pre")
    parser.add_argument("-m","--methy", type=str, help="reault methy pos file")
    parser.add_argument("-r","--read_m6A_dict", type=str, help="read m6A pkl file path")
    parser.add_argument("-c","--chrom", type=str, help="chromosome")
    parser.add_argument("-s","--strand", type=str, help="different strand processing")
    parser.add_argument("--geno_size", type=str, help="genome size file path")
    parser.add_argument("--base_minQ", type=int, default=5, help="base min qscore(default=5)")
    parser.add_argument("--read_minQ", type=int, default=0, help="read min qscore(default=0)")
    parser.add_argument("--snp_min_cov", type=int, default=5, help="ref or alt min coverage(default=5)")
    args = parser.parse_args()

    output_path = f"{args.outdirpre}_haplotype_{args.chrom}_{args.strand}_tmp.csv"

    base_minQ = args.base_minQ-1 if args.base_minQ != 0 else 0
    read_minQ = args.read_minQ-1 if args.read_minQ != 0 else 0

    snp_info = pd.read_csv(args.snp_info, sep="\t")
    snp_info.columns = ["chrom","pos1","rsID","A1","A2","EAF"]
    snp_info["pos0"] = snp_info["pos1"].astype(int) - 1
    snp_info["k"] = snp_info["chrom"] + "_" + snp_info["pos0"].astype(str)
    snp_info["v"] = snp_info["rsID"] + ";" + snp_info["A1"] + ";" + snp_info["A2"] + ";" + snp_info["EAF"].astype(str)
    snp_dict = dict(zip(snp_info["k"], snp_info["v"])) # chrom_pos0: rsid;A1;A2;EAF

    ## 读甲基化位点文件并转化成字典
    methy_df = pd.read_csv(args.methy)
    methy_df = methy_df[(methy_df['strand']==args.strand)&(methy_df['chrom']==args.chrom)]
    if len(methy_df) == 0:
        print(f'{args.chrom} 没有m6A位点')
        exit()
    methy_df['pos_0base'] = methy_df['pos_1base']-1
    methy_dict = dict(x for x in methy_df.apply(create_m6A_dict, axis=1).tolist() if x is not None)

    with open(args.read_m6A_dict, 'rb') as file:
        m6A_id_dict = pickle.load(file) # chrom_pos0_strand: read1;read2;...
    geno_size_df = pd.read_csv(args.geno_size, sep="\t", header=None, names=["chrom","size"])
    geno_size_dict = dict(zip(geno_size_df["chrom"], geno_size_df["size"]))
    haplotype_df = pd.DataFrame(columns=[
        "chrom","strand","m6A_pos_1base","snp_pos_1base",
        "rsID","A1","A2",
        "A1_A","A2_A","A1_m6A","A2_m6A",
        "ism6A?","methy_rate",
        "EAF"
        ])

    start,end = 0,int(geno_size_dict[args.chrom])
    haplotype_df = pd.concat([haplotype_df, count_haplotype(args.chrom, start, end, args.strand, args.bam, methy_dict, m6A_id_dict, snp_dict,base_minQ,read_minQ)],ignore_index=True)
    haplotype_df = haplotype_df[(haplotype_df['A1'] != '') & (haplotype_df['A2'] != '')]
    if len(haplotype_df) != 0:
        df = haplotype_df.sort_values(by=['chrom', 'm6A_pos_1base', 'snp_pos_1base'])
        df = df.reset_index(drop=True)
        results = df.apply(lambda row: apply_simulate_fisher(row, args.snp_min_cov), axis=1, result_type='expand')
        df[['p_value','beta', 'SE']] = results
        df = df[df['p_value'].notna()]
        df.to_csv(output_path, index=None)
        print(f"{output_path}已保存")
