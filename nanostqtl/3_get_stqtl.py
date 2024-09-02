import pandas as pd
import pysam
import argparse
from collections import Counter
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy import stats
import statsmodels.api as sm


def count_haplotype(chrom, end, strand, bamfile, read_st_dict, snp_file_dict,threadsnum,base_minQ,read_minQ):
    # read_st_dict read对应的稳定性指标 readid:stability
    # snp_file_dict snp位点信息 # pos0: rsid;A1;A2;EAF
    snp_bases = {}
    # snp_bases 存放的是snp位点 [pos0]:{read1:"A",read2:"A",read3:"T"}
    with pysam.AlignmentFile(bamfile, "rb",threads=threadsnum) as samfile:
        for pileupcolumn in samfile.pileup(chrom, 0, end,min_base_quality=base_minQ,min_mapping_quality=read_minQ,stepper="samtools",max_depth=5000000):
            base_pos = pileupcolumn.reference_pos
            ## 处理snp位点的情况
            # 1判断位置是否是snp
            if base_pos in snp_file_dict:
                # 2判断是否是杂合子位点
                seqs = [i.upper() for i in pileupcolumn.get_query_sequences()]  # 被覆盖位点处的所有碱基
                if len(set([seq for seq in seqs if seq])) > 1:
                    snp_bases[base_pos] = {
                        i: j for i, j in zip(pileupcolumn.get_query_names(), seqs)
                    }
    if len(snp_bases) == 0:
        print(f'{chrom} {strand}中,没有杂合的snp被read覆盖')
        return None
    haplotype_df = pd.DataFrame(columns=[
        "chrom","strand","snp_pos_1base",
        "rsID","A1","A2","EAF",
        "A1_STscore_l","A2_STscore_l",
        ])
    result = start_get_haplotypes(snp_bases, read_st_dict, snp_file_dict)
    if result is None:
        print(f'{chrom} {strand}中没有符合条件的snp')
        return None
    else:
        ith = 0
        for i in result:
            haplotype_df.loc[ith] = (
                chrom,strand,int(i[6])+1,
                i[4],i[0],i[1],i[5],
                i[2],i[3],
            )
            ith += 1
        return haplotype_df


def start_get_haplotypes(snp_bases,read_st_dict,snp_file_dict):
    res_l = []
    # snp_bases [pos0]:{read1:"A",read2:"A",read3:"T"}
    # read_st_dict readid:stability
    for snp_pos, snp_base in snp_bases.items():
        # 筛选出snp覆盖的read有稳定性分数的read
        snp_base_filter = {key: value for key, value in snp_base.items() if key in read_st_dict}
        if len(snp_base_filter) == 0:
            continue
        res = get_haplotypes(snp_base_filter, read_st_dict, snp_file_dict[snp_pos])
        if res is not None:
            a1,a2,A1_STscore_l, A2_STscore_l,snpID,eaf = res
            res_l.append([a1,a2,A1_STscore_l, A2_STscore_l,snpID,eaf,snp_pos])
    if len(res_l) != 0:
        return res_l
    else:
        return None


def get_haplotypes(snp_base,read_TSS_dict,snp_file_dict_res):
    snpID, A1, A2, eaf = snp_file_dict_res.split(";")[:4]
    eaf = float(eaf)
    A1_read = set(k for k, v in snp_base.items() if v == A1)
    A2_read = set(k for k, v in snp_base.items() if v == A2)
    if A1_read and A2_read:
        A1_STscore_l = [read_TSS_dict[read_id] for read_id in A1_read if read_id in read_TSS_dict]
        A2_STscore_l = [read_TSS_dict[read_id] for read_id in A2_read if read_id in read_TSS_dict]
        if len(A1_STscore_l) == 0 or len(A2_STscore_l) == 0:
            return None
        return A1,A2,A1_STscore_l, A2_STscore_l,snpID,eaf
    else:
        return None


# 线性回归
def analyze_snp_stability(A1_stability, A2_stability):
    total_count = len(A1_stability) + len(A2_stability)
    A1_count = len(A1_stability)
    A2_count = len(A2_stability)
    if (total_count <= 20) or (A1_count == 0) or (A2_count == 0) or (min(A1_count, A2_count) / max(A1_count, A2_count) < 0.1):
        return None, None, None
    stability_score = A1_stability + A2_stability
    SNP_type = [1] * A1_count + [0] * A2_count
    df = pd.DataFrame({
        'stability_score': stability_score,
        'SNP_type': SNP_type
    })
    # 添加常数项
    X = sm.add_constant(df['SNP_type'])
    # 拟合线性回归模型
    model = sm.OLS(df['stability_score'], X).fit()
    # 提取关键结果
    beta = model.params['SNP_type']  # 回归系数
    se = model.bse['SNP_type']  # 标准误
    p_value = model.pvalues['SNP_type']
    return p_value, beta, se


def process_all_data(df):
    results = df.apply(
        lambda row: analyze_snp_stability(
            row["A1_STscore_l"], row["A2_STscore_l"]
        ),
        axis=1,
    )
    (
        df["pvalue"],
        df["beta"],
        df["SE"],
    ) = zip(*results)
    df = df[(df["pvalue"].notna())]
    if len(df) == 0:
        return None
    df["FDR"] = multipletests(df["pvalue"], method="fdr_bh")[1]
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Nanopore dapaect RNA data call RNA stability QTL(stqtl).')
    parser.add_argument("-b","--bam", type=str, help="bam file path")
    parser.add_argument("--snp_info", type=str, help="processed snp pkl file path")
    parser.add_argument("-o","--outdirpre", type=str, help="outdir and pre")
    parser.add_argument("--readST", type=str, help="read stability csv file")
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

    # csv2dict
    read_st = pd.read_csv(args.readST, usecols = [0,5])
    read_st_dict = read_st.set_index('readID')['stability_zscore'].to_dict() # readid:stability
    # txt2dict
    snp_info = pd.read_csv(args.snp_info, sep="\t")
    snp_info.columns = ["chrom","pos1","rsID","A1","A2","EAF"]
    snp_info = snp_info[snp_info['chrom'] == args.chrom]
    snp_info["pos0"] = snp_info["pos1"].astype(int) - 1
    snp_info["v"] = snp_info["rsID"] + ";" + snp_info["A1"] + ";" + snp_info["A2"] + ";" + snp_info["EAF"].astype(str)
    snp_dict = dict(zip(snp_info["pos0"], snp_info["v"])) # pos0: rsid;A1;A2;EAF

    geno_size_df = pd.read_csv(args.geno_size, sep="\t", header=None, names=["chrom","size"])
    geno_size_dict = dict(zip(geno_size_df["chrom"], geno_size_df["size"]))
    haplotype_df = pd.DataFrame(columns=[
        "chrom","strand","snp_pos_1base",
        "rsID","A1","A2","EAF",
        "A1_STscore_l","A2_STscore_l",
        ])
    end = geno_size_dict[args.chrom]
    haplotype_df = pd.concat([haplotype_df, count_haplotype(args.chrom, end, args.strand, args.bam, read_st_dict, snp_dict, args.threads, base_minQ, read_minQ)],ignore_index=True)
    if len(haplotype_df) != 0:
        df = haplotype_df.sort_values(by=["chrom", "snp_pos_1base"])
        df = process_all_data(df)
        if df is None:
            print(
                f"{args.chrom} {args.strand}中的SNP没有足够的read覆盖度，无法进行统计"
            )
            exit()
        df = df.reset_index(drop=True)
        df.to_csv(output_path, index=None)
        print(f"{output_path}已保存")