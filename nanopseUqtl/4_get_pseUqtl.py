import pandas as pd
import pysam
import argparse
import pickle
from collections import Counter
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests


def create_pseU_dict(row):
    if row['mod_rate'] != 100:
        return row["pos_0base"], row["mod_rate"]


def count_haplotype(chrom, start, end, strand, bamfile, mod_dict, pseU_id_dict, snp_file_dict, base_minQ, read_minQ, threadsnum):
    # mod_dict 甲基化位点信息 int(pos0):mod_rate
    # pseU_id_dict 含有甲基化的read chrom_pos0_strand: read1;read2;...
    # snp_file_dict snp位点信息 int(pos0): rsid;A1;A2;EAF
    pseU_all_reads,snp_bases,snp_pseU_bases = {},{},{}
    # pseU_all_reads {pos0}:[read1,read2,read3]
    # snp_bases,snp_pseU_bases {pos0}:{read1:"A",read2:"A",read3:"T"}
    with pysam.AlignmentFile(bamfile, "rb",threads=threadsnum) as samfile:
        for pileupcolumn in samfile.pileup(chrom, start, end, min_base_quality=base_minQ, min_mapping_quality=read_minQ, stepper="samtools", max_depth=5000000):
            base_pos = pileupcolumn.reference_pos
            # 首先获取pseU处的所有read，因为存在pseU处是N的情况，所以这里要筛选一下
            if base_pos in mod_dict:
                seqs = [i.upper() for i in pileupcolumn.get_query_sequences()]
                pseU_all_reads[base_pos] = [
                    i for i, j in zip(pileupcolumn.get_query_names(), seqs) if j
                ]
            # 再获取所有的SNP信息
            ## 1判断位置是否是snp
            if base_pos in snp_file_dict:
                ## 2判断是否是杂合子位点
                seqs = [i.upper() for i in pileupcolumn.get_query_sequences()]
                if len(set([seq for seq in seqs if seq])) > 1:
                    ## 3判断snp是否是pseU位点
                    if base_pos in mod_dict:
                        snp_pseU_bases[base_pos] = {
                            i: j for i, j in zip(pileupcolumn.get_query_names(), seqs) if j
                        }
                    else:
                        snp_bases[base_pos] = {
                            i: j for i, j in zip(pileupcolumn.get_query_names(), seqs) if j
                        }
    if len(snp_pseU_bases) == 0 and len(snp_bases) == 0:
        print(f"{chrom} {strand} 没有杂合的SNP位点")
        return None
    haplotype_df = pd.DataFrame(columns=[
        "chrom","strand","snp_pos_1base","pseU_pos_1base",
        "rsID","A1","A2","EAF",
        "A1_U","A2_U","A1_pseU","A2_pseU",
        "ispseU?","mod_rate"
        ])
    ith = 0
    bases_list = [snp_bases, snp_pseU_bases]
    ispseU_list = ["", "Yes"]
    for bases, ispseU in zip(bases_list, ispseU_list):
        result = start_get_haplotypes(chrom, strand, pseU_all_reads, bases, pseU_id_dict, snp_file_dict)
        if result is None:
            continue
        else:
            for i in result:
                haplotype_df.loc[ith] = (
                    chrom,strand,i[0] + 1,i[1] + 1,
                    i[8],i[2],i[3],i[9],
                    i[4],i[5],i[6],i[7],
                    ispseU,mod_dict[i[1]]
                )
                ith += 1
    if len(haplotype_df) != 0:
        return haplotype_df
    else:
        print(f"{chrom} {strand} 没有A1、A2和U、paeU对应的read，无法创建列联表")
        return None


def start_get_haplotypes(chrom, strand, pseU_all_reads, snp_bases, pseU_id_dict, snp_file_dict):
    res_l = []
    # snp_bases [pos0]:{read1:"A",read2:"A",read3:"T"}
    # pseU_id_dict chrom_pos0_strand: read1;read2;... 所有甲基化的read
    for snp_pos, snp_base in snp_bases.items():
        for pseU_pos, pseU_all_reads_l in pseU_all_reads.items():
            pseU_readid = set(pseU_id_dict[f"{chrom}_{str(pseU_pos)}_{strand}"].split(";"))
            unpseU_readid = set(pseU_all_reads_l)- pseU_readid
            res = get_haplotypes(snp_base, pseU_readid, unpseU_readid, snp_file_dict[snp_pos])
            if res is not None:
                A1,A2,A1_U,A2_U,A1_pseU,A2_pseU,snpID,eaf = res
                res_l.append([snp_pos,pseU_pos,A1,A2,A1_U,A2_U,A1_pseU,A2_pseU,snpID,eaf])
    if len(res_l) != 0:
        return res_l
    else:
        return None


def get_haplotypes(snp_base, pseU_readid, unpseU_readid, snp_file_dict_res):
    intersect_pseU_read_set = set(snp_base.keys()) & pseU_readid  # SNP覆盖的所有被修饰的read
    intersect_unpseU_read_set = set(snp_base.keys()) & unpseU_readid  # SNP覆盖的所有没有修饰的read
    if (len(intersect_pseU_read_set) == 0) or (len(intersect_unpseU_read_set) == 0):
        return None
    snpID, A1, A2, eaf = snp_file_dict_res.split(";")[:4]
    eaf = float(eaf)
    A1_read = set(k for k, v in snp_base.items() if v == A1)
    A2_read = set(k for k, v in snp_base.items() if v == A2)
    if A1_read and A2_read:
        A1_U = len(A1_read & intersect_unpseU_read_set)
        A2_U = len(A2_read & intersect_unpseU_read_set)
        A1_pseU = len(A1_read & intersect_pseU_read_set)
        A2_pseU = len(A2_read & intersect_pseU_read_set)
        return A1,A2,A1_U,A2_U,A1_pseU,A2_pseU,snpID,eaf
    else:
        return None


def analyze_snp_methylation_bayes(A1_U, A2_U, A1_pseU, A2_pseU, n_samples=10000000, random_state=42):
    # if A1_U>0 and A2_U>0 and A1_pseU>0 and A2_pseU>0 and (A1_U+A2_U+A1_pseU+A2_pseU)>100:
    values = [A1_U, A2_U, A1_pseU, A2_pseU]
    zero_count = values.count(0)
    if zero_count <= 1:
        # 使用Beta分布作为先验和后验
        post_alpha_ref = 1 + A1_pseU
        post_beta_ref = 1 + A1_U
        post_alpha_alt = 1 + A2_pseU
        post_beta_alt = 1 + A2_U
        # 设置随机种子
        rng = np.random.default_rng(random_state)
        # 从后验分布中抽样
        theta_ref_samples = stats.beta.rvs(post_alpha_ref, post_beta_ref, size=n_samples, random_state=rng)
        theta_alt_samples = stats.beta.rvs(post_alpha_alt, post_beta_alt, size=n_samples, random_state=rng)
        # 计算效应大小
        effect_size = np.mean(theta_ref_samples - theta_alt_samples)
        # 计算 beta（标准化效应大小）
        pooled_sd = np.sqrt((np.var(theta_ref_samples) + np.var(theta_alt_samples)) / 2)
        beta = effect_size / pooled_sd
        # 计算 beta 的标准误
        se = np.std(theta_ref_samples - theta_alt_samples) / pooled_sd
        # 计算 z 分数
        z_score = beta / se
        # 计算传统p值（双侧）
        p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
        return p_value, z_score, beta, se
    else:
        return None, None, None, None


def process_all_data(df):
    results = df.apply(lambda row: analyze_snp_methylation_bayes(row['A1_U'], row['A2_U'], row['A1_pseU'], row['A2_pseU']), axis=1)
    df['p_value'],df['z_score'],df['beta'],df['SE'] = zip(*results)
    df = df[(df['p_value'].notna())]
    df['FDR'] = multipletests(df['p_value'], method='fdr_bh')[1]
    return df



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Nanopore direct RNA data call pseUQTL.')
    parser.add_argument("-b","--bam", type=str, help="bam file path")
    parser.add_argument("--snp_info", type=str, help="processed snp pkl file path")
    parser.add_argument("-o","--outdirpre", type=str, help="outdir and pre")
    parser.add_argument("-m","--mod", type=str, help="result mod pos file")
    parser.add_argument("-r","--read_pseU_dict", type=str, help="read pseU pkl file path")
    parser.add_argument("-c","--chrom", type=str, help="chromosome")
    parser.add_argument("-s","--strand", type=str, help="different strand processing")
    parser.add_argument("--geno_size", type=str, help="genome size file path")
    parser.add_argument("-t","--threads", type=int, default=4, help="threads number (default: 4)")
    parser.add_argument("--base_minQ", type=int, default=5, help="base min qscore(default=5)")
    parser.add_argument("--read_minQ", type=int, default=0, help="read min qscore(default=0)")
    args = parser.parse_args()

    output_path = f"{args.outdirpre}_haplotype_{args.chrom}_{args.strand}_tmp.csv"

    base_minQ = args.base_minQ-1 if args.base_minQ != 0 else 0
    read_minQ = args.read_minQ-1 if args.read_minQ != 0 else 0

    snp_info = pd.read_csv(args.snp_info, sep="\t")
    snp_info.columns = ["chrom","pos1","rsID","A1","A2","EAF"]
    snp_info = snp_info[snp_info['chrom'] == args.chrom]
    snp_info["pos0"] = snp_info["pos1"].astype(int) - 1
    snp_info["v"] = snp_info["rsID"] + ";" + snp_info["A1"] + ";" + snp_info["A2"] + ";" + snp_info["EAF"].astype(str)
    snp_dict = dict(zip(snp_info["pos0"], snp_info["v"])) # pos0: rsid;A1;A2;EAF

    mod_df = pd.read_csv(args.mod)
    mod_df = mod_df[(mod_df['strand']==args.strand)&(mod_df['chrom']==args.chrom)]
    if len(mod_df) == 0:
        print(f'{args.chrom} {args.strand} 没有pseU位点')
        exit()
    mod_df['pos_0base'] = mod_df['pos_1base']-1
    mod_dict = dict(x for x in mod_df.apply(create_pseU_dict, axis=1).tolist() if x is not None) # pos_0base:mod_rate

    with open(args.read_pseU_dict, 'rb') as file:
        pseU_id_dict = pickle.load(file) # chrom_pos0_strand: read1;read2;...
    geno_size_df = pd.read_csv(args.geno_size, sep="\t", header=None, names=["chrom","size"])
    geno_size_dict = dict(zip(geno_size_df["chrom"], geno_size_df["size"]))
    haplotype_df = pd.DataFrame(columns=[
        "chrom","strand","snp_pos_1base","pseU_pos_1base",
        "rsID","A1","A2","EAF",
        "A1_U","A2_U","A1_pseU","A2_pseU",
        "ispseU?","mod_rate",
        ])

    start,end = 0,int(geno_size_dict[args.chrom])
    haplotype_df = pd.concat([haplotype_df, count_haplotype(args.chrom, start, end, args.strand, args.bam, mod_dict, pseU_id_dict, snp_dict, base_minQ, read_minQ, args.threads)],ignore_index=True)
    haplotype_df = haplotype_df[(haplotype_df['A1'] != '') & (haplotype_df['A2'] != '')]
    if len(haplotype_df) != 0:
        df = haplotype_df.sort_values(by=['chrom', 'snp_pos_1base', 'pseU_pos_1base'])
        df = process_all_data(df)
        df = df.reset_index(drop=True)
        df.to_csv(output_path, index=None)
        print(f"{output_path}已保存")
