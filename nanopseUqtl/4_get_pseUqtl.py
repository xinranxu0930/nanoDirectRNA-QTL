import pandas as pd
import pysam
import argparse
import pickle
from collections import Counter
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import fisher_exact


def create_pseU_dict(row):
    if row['mod_rate'] != 100:
        return f'{row["pos_0base"]}', row["mod_rate"]

def get_reverse_complementary_sequence(seq):
    seqreverse = seq[::-1]
    transtable = str.maketrans('ATGCatgc','TACGtacg')
    finalseq = seqreverse.translate(transtable).upper()
    return finalseq

## call pseUqtl
def count_haplotype(chrom, start, end, strand, bamfile, mod_dict, pseU_id_dict, snp_file_dict,base_minQ,read_minQ):
    # mod_dict 甲基化位点信息 pos0:mod_rate
    # pseU_id_dict 含有甲基化的read chrom_pos0_strand: read1;read2;...
    # snp_file_dict snp位点信息 chrom_pos0: rsid;ref;alt
    pseU_bases,snp_bases,snp_pseU_bases = {},{},{}
    # pseU_bases 存放目标区域中甲基化位点处所有的read id，这里同时包括甲基化的read和非甲基化的read，和pseU_id_dict联合使用
    # 后两个 XXX_bases[pos0]:{read1:"A",read2:"A",read3:"T"}
    # snp_bases 存放的是snp位点
    # snp_pseU_bases 存放的位点同时是snp和甲基化位点
    with pysam.AlignmentFile(bamfile, "rb",threads=10) as samfile:
        for pileupcolumn in samfile.pileup(chrom, start, end,min_base_quality=base_minQ,min_mapping_quality=read_minQ,stepper="samtools",max_depth=50000):
            base_pos = str(pileupcolumn.reference_pos)
            ## 先处理甲基化位点的情况
            if base_pos in mod_dict:
                pseU_bases[base_pos] = pileupcolumn.get_query_names()
            ## 处理snp位点的情况
            # 1判断位置是否是snp
            if f"{chrom}_{base_pos}" in snp_file_dict:
                # 2判断是否是杂合子位点
                seqs = [i.upper() for i in pileupcolumn.get_query_sequences()]  # 被覆盖位点处的所有碱基
                if len(set(seqs)) > 1:
                    # 3判断snp是否是pseU位点
                    if base_pos in pseU_bases:
                        snp_pseU_bases[base_pos] = {
                            i: j for i, j in zip(pileupcolumn.get_query_names(), seqs)
                        }
                    else:
                        snp_bases[base_pos] = {
                            i: j for i, j in zip(pileupcolumn.get_query_names(), seqs)
                        }
    if len(pseU_bases) == 0 and (len(snp_pseU_bases)==0 or len(snp_bases) == 0):
        return None
    haplotype_df = pd.DataFrame(columns=[
        "chrom","strand","pseU_pos_1base","snp_pos_1base",
        "rsID","A1","A2",
        "A1_U","A2_U","A1_pseU","A2_pseU",
        "ispseU?","mod_rate",
        "EAF",
        ])
    ith = 0
    bases_list = [snp_bases, snp_pseU_bases]
    ispseU_list = ["", "Yes"]
    for bases, ispseU in zip(bases_list, ispseU_list):
        result = start_get_haplotypes(chrom, strand, pseU_bases, bases, pseU_id_dict, snp_file_dict)
        if result is None:
            continue
        else:
            for i in result:
                haplotype_df.loc[ith] = (
                    chrom,strand,int(i[0]) + 1,int(i[1]) + 1,
                    i[8],i[2],i[3],
                    i[4],i[5],i[6],i[7],
                    ispseU,mod_dict[str(i[0])],
                    i[9]
                )
                ith += 1
    if len(haplotype_df) != 0:
        return haplotype_df
    else:
        return None

def start_get_haplotypes(chrom,strand,pseU_bases,snp_bases,pseU_id_dict,snp_file_dict):
    res_l = [] # 每一个元素都表示一个位点的结果 [pseU_pos0,snp_0base,ref,alt,A1_U,A2_U,A1_pseU,A2_pseU]
    # snp_bases [pos0]:{read1:"A",read2:"A",read3:"T"}
    # pseU_bases [base_pos]:[readid1,readid12,...] 位点处所有的read 甲基化 非甲基化都有
    # pseU_id_dict chrom_pos0_strand: read1;read2;... 所有甲基化的read
    for pseU_pos, pseU_allreadid_list in pseU_bases.items():
        for snp_pos, snp_base in snp_bases.items():
            pseU_readid = set(pseU_id_dict[f"{chrom}_{pseU_pos}_{strand}"].split(";")) # 某个点只是甲基化的所有read
            res = get_haplotypes(snp_base, pseU_allreadid_list, pseU_readid, snp_file_dict[f"{chrom}_{snp_pos}"],strand)
            if res is not None:
                a1,a2,a1_U,a2_U,a1_pseU,a2_pseU,snpID,eaf = res
                res_l.append([pseU_pos,snp_pos,a1,a2,a1_U,a2_U,a1_pseU,a2_pseU,snpID,eaf])
    if len(res_l) != 0:
        return res_l
    else:
        return None

# 现在是确定了snp-pseU，输入snp处的readid-base；pseU出所有readid；pseU处甲基化的readid；该snp的rsid ref alt；strand
def get_haplotypes(snp_base, pseU_allreadid_list, pseU_readid,snp_file_dict_res,strand):
    intersect_read_set = set(snp_base.keys()) & set(pseU_allreadid_list)  # 两个字典共有的read id的集合，待分析的所有read
    if len(intersect_read_set) == 0:
        return None
    in_pseU_readid_set = intersect_read_set & pseU_readid # 含有甲基化的read
    not_in_pseU_readid_set = intersect_read_set - pseU_readid # 含有正常A的read
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
    a1_U = len(a1_read & not_in_pseU_readid_set)
    a2_U = len(a2_read & not_in_pseU_readid_set)
    a1_pseU = len(a1_read & in_pseU_readid_set)
    a2_pseU = len(a2_read & in_pseU_readid_set)
    return a1,a2,a1_U,a2_U,a1_pseU,a2_pseU,snpID,eaf

# fisher
def simulate_fisher_exact_se(A1_U, A2_U, A1_pseU, A2_pseU, num_simulations=100):
    observed = np.array([[A1_U, A2_U], [A1_pseU, A2_pseU]])
    if 0 in [A1_U, A2_U, A1_pseU, A2_pseU]:
        odds_ratio = ((A1_U + 0.75) * (A2_pseU + 0.75)) / ((A2_U + 0.75) * (A1_pseU + 0.75))
        p_value = fisher_exact(observed)[1]
    else:
        odds_ratio, p_value = fisher_exact(observed)
    n = int(A1_U+A2_U+A1_pseU+A2_pseU)
    if np.isinf(odds_ratio) or odds_ratio == 0 or n == 0:
        return p_value, n, None, None
    beta = np.log(odds_ratio)
    # Simulate to estimate the standard error
    odds_ratios = np.zeros(num_simulations)
    valid_simulations = 0
    for i in range(num_simulations):
        # Simulate contingency table based on margins
        n1 = A1_U + A2_U
        n2 = A1_pseU + A2_pseU
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

def simulate_fisher_exact(A1_U, A2_U, A1_pseU, A2_pseU, min_cov):
    p_value = None
    beta = None
    se_log_or = None
    if (A1_U+A1_pseU >= min_cov) & (A2_U+A2_pseU >= min_cov):
        A1_U += 0.75 if A1_U == 0 else 0
        A2_U += 0.75 if A2_U == 0 else 0
        A1_pseU += 0.75 if A1_pseU == 0 else 0
        A2_pseU += 0.75 if A2_pseU == 0 else 0
        observed = np.array([[A1_U, A2_U], [A1_pseU, A2_pseU]])
        odds_ratio = (A1_U * A2_pseU) / (A2_U * A1_pseU)
        se_log_or = np.sqrt(1/A1_U + 1/A2_pseU + 1/A2_U + 1/A1_pseU)
        p_value = fisher_exact(observed)[1]
        if np.isinf(odds_ratio) or odds_ratio == 0:
            beta = None
        else:
            beta = np.log(odds_ratio)
    return p_value, beta, se_log_or

def apply_simulate_fisher(row,min_cov):
    return simulate_fisher_exact(row['A1_U'], row['A2_U'], row['A1_pseU'], row['A2_pseU'],min_cov)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Nanopore direct RNA data call pseUqtl.')
    parser.add_argument("-b","--bam", type=str, help="bam file path")
    parser.add_argument("--snp_info", type=str, help="processed snp pkl file path")
    parser.add_argument("-p","--outdirpre", type=str, help="outdir and pre")
    parser.add_argument("-m","--mod", type=str, help="reault mod pos file")
    parser.add_argument("-r","--read_pseU_dict", type=str, help="read pseU pkl file path")
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
    mod_df = pd.read_csv(args.mod)
    mod_df = mod_df[(mod_df['strand']==args.strand)&(mod_df['chrom']==args.chrom)]
    if len(mod_df) == 0:
        print(f'{args.chrom} 没有pseU位点')
        exit()
    mod_df['pos_0base'] = mod_df['pos_1base']-1
    mod_dict = dict(x for x in mod_df.apply(create_pseU_dict, axis=1).tolist() if x is not None)

    with open(args.read_pseU_dict, 'rb') as file:
        pseU_id_dict = pickle.load(file) # chrom_pos0_strand: read1;read2;...
    geno_size_df = pd.read_csv(args.geno_size, sep="\t", header=None, names=["chrom","size"])
    geno_size_dict = dict(zip(geno_size_df["chrom"], geno_size_df["size"]))
    haplotype_df = pd.DataFrame(columns=[
        "chrom","strand","pseU_pos_1base","snp_pos_1base",
        "rsID","A1","A2",
        "A1_U","A2_U","A1_pseU","A2_pseU",
        "ispseU?","mod_rate",
        "EAF",
        ])

    start,end = 0,int(geno_size_dict[args.chrom])
    haplotype_df = pd.concat([haplotype_df, count_haplotype(args.chrom, start, end, args.strand, args.bam, mod_dict, pseU_id_dict, snp_dict, base_minQ, read_minQ)],ignore_index=True)
    haplotype_df = haplotype_df[(haplotype_df['A1'] != '') & (haplotype_df['A2'] != '')]
    if len(haplotype_df) != 0:
        df = haplotype_df.sort_values(by=['chrom', 'pseU_pos_1base', 'snp_pos_1base'])
        df = df.reset_index(drop=True)
        results = df.apply(lambda row: apply_simulate_fisher(row, args.snp_min_cov), axis=1, result_type='expand')
        df[['p_value','beta', 'SE']] = results
        df = df[df['p_value'].notna()]
        df.to_csv(output_path, index=None)
        print(f"{output_path}已保存")
