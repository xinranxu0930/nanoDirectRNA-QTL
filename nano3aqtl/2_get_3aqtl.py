import pandas as pd
import pysam
import scipy
import argparse
from collections import Counter
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy import stats


def count_haplotype(
    chrom,
    end,
    strand,
    bamfile,
    read_APA_dict,
    snp_file_dict,
    threadsnum,
    base_minQ,
    read_minQ,
):
    # read_APA_dict read对应的APA read:APA
    # snp_file_dict snp位点信息 int(pos0): rsid;A1;A2;EAF
    snp_bases = {}
    # snp_bases 存放的是snp位点 [pos0]:{read1:"A",read2:"A",read3:"T"}
    with pysam.AlignmentFile(bamfile, "rb", threads=threadsnum) as samfile:
        for pileupcolumn in samfile.pileup(
            chrom,
            0,
            end,
            min_base_quality=base_minQ,
            min_mapping_quality=read_minQ,
            stepper="samtools",
            max_depth=5000000,
        ):
            base_pos = pileupcolumn.reference_pos
            ## 处理snp位点的情况
            # 1判断位置是否是snp
            if base_pos in snp_file_dict:
                # 2判断是否是杂合子位点
                seqs = [i.upper() for i in pileupcolumn.get_query_sequences()]
                if len(set([seq for seq in seqs if seq])) > 1:
                    snp_bases[base_pos] = {
                        i: j for i, j in zip(pileupcolumn.get_query_names(), seqs)
                    }
    if len(snp_bases) == 0:
        print(f"{chrom} {strand}中,没有杂合的snp被read覆盖")
        return None
    haplotype_df = pd.DataFrame(
        columns=[
            "chrom",
            "strand",
            "snp_pos_1base",
            "rsID",
            "A1",
            "A2",
            "EAF",
            "A1_APA_count",
            "A2_APA_count",
            "APA_id",
        ]
    )
    result = start_get_haplotypes(snp_bases, read_APA_dict, snp_file_dict)
    if result is None:
        print(f"{chrom} {strand}中没有符合条件的snp")
        return None
    else:
        ith = 0
        for i in result:
            haplotype_df.loc[ith] = (
                chrom,
                strand,
                i[7] + 1,
                i[5],
                i[0],
                i[1],
                i[6],
                i[2],
                i[3],
                i[4],
            )
            ith += 1
        return haplotype_df


def start_get_haplotypes(snp_bases, read_APA_dict, snp_file_dict):
    res_l = []
    # snp_bases [pos0]:{read1:"A",read2:"A",read3:"T"}
    # read_APA_dict read:APA
    for snp_pos, snp_base in snp_bases.items():
        # 筛选出snp覆盖的read中有APA分类的read
        snp_base_filter = {
            key: value for key, value in snp_base.items() if key in read_APA_dict
        }
        if len(snp_base_filter) == 0:
            continue
        # res = get_haplotypes_fourfold(snp_base_filter, read_APA_dict, snp_file_dict[snp_pos])
        res = get_haplotypes_Multidimensional(
            snp_base_filter, read_APA_dict, snp_file_dict[snp_pos]
        )
        if res is not None:
            A1, A2, A1_count_l, A2_count_l, all_APA_type, snpID, eaf = res
            res_l.append(
                [A1, A2, A1_count_l, A2_count_l, all_APA_type, snpID, eaf, snp_pos]
            )
    if len(res_l) != 0:
        return res_l
    else:
        return None


def get_haplotypes_fourfold(snp_base, read_APA_dict, snp_file_dict_res):
    snpID, A1, A2, eaf = snp_file_dict_res.split(";")[:4]
    eaf = float(eaf)
    A1_read = set(k for k, v in snp_base.items() if v == A1)
    A2_read = set(k for k, v in snp_base.items() if v == A2)
    if A1_read and A2_read:
        A1_apa = Counter(
            read_APA_dict[read_id] for read_id in A1_read if read_id in read_APA_dict
        )
        A2_apa = Counter(
            read_APA_dict[read_id] for read_id in A2_read if read_id in read_APA_dict
        )
        all_apa_types = set(A1_apa.keys()) | set(A2_apa.keys())
        data = {
            "A1": [A1_apa.get(apa, 0) for apa in all_apa_types],
            "A2": [A2_apa.get(apa, 0) for apa in all_apa_types],
        }
        df = pd.DataFrame(data, index=list(all_apa_types))
        if len(df) == 1:
            return None
        elif len(df) > 2:
            total_counts = df.sum(axis=1)
            top_apa = total_counts.idxmax()
            new_df = pd.DataFrame(
                {
                    "A1": [df.loc[top_apa, "A1"], df.drop(top_apa).sum()["A1"]],
                    "A2": [df.loc[top_apa, "A2"], df.drop(top_apa).sum()["A2"]],
                },
                index=[top_apa, "others"],
            )
            return (
                A1,
                A2,
                new_df["A1"].tolist(),
                new_df["A2"].tolist(),
                df.index.tolist(),
                snpID,
                eaf,
            )
        else:
            return (
                A1,
                A2,
                df["A1"].tolist(),
                df["A2"].tolist(),
                df.index.tolist(),
                snpID,
                eaf,
            )


def get_haplotypes_Multidimensional(snp_base, read_APA_dict, snp_file_dict_res):
    snpID, A1, A2, eaf = snp_file_dict_res.split(";")[:4]
    eaf = float(eaf)
    A1_read = set(k for k, v in snp_base.items() if v == A1)
    A2_read = set(k for k, v in snp_base.items() if v == A2)
    if A1_read and A2_read:
        A1_apa = Counter(
            read_APA_dict[read_id] for read_id in A1_read if read_id in read_APA_dict
        )
        A2_apa = Counter(
            read_APA_dict[read_id] for read_id in A2_read if read_id in read_APA_dict
        )
        all_apa_types = set(A1_apa.keys()) | set(A2_apa.keys())
        data = {
            "A1": [A1_apa.get(apa, 0) for apa in all_apa_types],
            "A2": [A2_apa.get(apa, 0) for apa in all_apa_types],
        }
        df = pd.DataFrame(data, index=list(all_apa_types))
        if len(df) == 1:
            return None
        else:
            return (
                A1,
                A2,
                df["A1"].tolist(),
                df["A2"].tolist(),
                df.index.tolist(),
                snpID,
                eaf,
            )


def analyze_snp_methylation_bayes_fourfold(
    A1_APA_count, A2_APA_count, n_samples=10000000, random_state=42
):
    A1_APA1, A2_APA1, A1_APA2, A2_APA2 = (
        A1_APA_count[0],
        A2_APA_count[0],
        A1_APA_count[1],
        A2_APA_count[1],
    )
    if (
        A1_APA1 > 0
        and A2_APA1 > 0
        and A1_APA2 > 0
        and A2_APA2 > 0
        and (A1_APA1 + A2_APA1 + A1_APA2 + A2_APA2) > 20
    ):
        # 使用Beta分布作为先验和后验
        post_alpha_ref = 1 + A1_APA2
        post_beta_ref = 1 + A1_APA1
        post_alpha_alt = 1 + A2_APA2
        post_beta_alt = 1 + A2_APA1
        # 设置随机种子
        rng = np.random.default_rng(random_state)
        # 从后验分布中抽样
        theta_ref_samples = stats.beta.rvs(
            post_alpha_ref, post_beta_ref, size=n_samples, random_state=rng
        )
        theta_alt_samples = stats.beta.rvs(
            post_alpha_alt, post_beta_alt, size=n_samples, random_state=rng
        )
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


def process_all_data_fourfold(df):
    results = df.apply(
        lambda row: analyze_snp_methylation_bayes_fourfold(
            row["A1_APA_count"], row["A2_APA_count"]
        ),
        axis=1,
    )
    df["p_value"], df["z_score"], df["beta"], df["SE"] = zip(*results)
    df = df[(df["p_value"].notna())]
    if len(df) == 0:
        return None
    df["FDR"] = multipletests(df["p_value"], method="fdr_bh")[1]
    return df


def analyze_snp_methylation_bayes_Multidimensional(
    A1_counts, A2_counts, APA_types, n_samples=10000000, random_state=42
):
    A1_counts = np.array(A1_counts)
    A2_counts = np.array(A2_counts)
    # 检查总样本量和非零类别数
    total_count = np.sum(A1_counts) + np.sum(A2_counts)
    non_zero_categories = np.sum((A1_counts > 0) & (A2_counts > 0))
    if (total_count > 20) & (non_zero_categories >= 2):
        # 设置无信息先验 Dirichlet分布
        Dirichlet_prior = np.ones_like(A1_counts)
        # 计算后验参数
        A1_posterior = Dirichlet_prior + A1_counts
        A2_posterior = Dirichlet_prior + A2_counts
        # 设置随机种子
        rng = np.random.default_rng(random_state)
        # 从后验分布中抽样
        A1_samples = rng.dirichlet(A1_posterior, size=n_samples)
        A2_samples = rng.dirichlet(A2_posterior, size=n_samples)
        # 计算效应大小（各类别比例的差异）
        effect_sizes = A1_samples - A2_samples
        # 计算平均效应大小
        mean_effect = np.mean(effect_sizes, axis=0)
        # 计算合并标准差
        pooled_sd = np.sqrt(
            (np.var(A1_samples, axis=0) + np.var(A2_samples, axis=0)) / 2
        )
        # 标准化效应大小 (beta)
        beta = mean_effect / pooled_sd
        # 计算标准误
        se = np.std(effect_sizes, axis=0) / pooled_sd
        # 计算 Z 分数
        z_score = beta / se
        # 计算 p 值（双侧）
        p_values = 2 * (1 - stats.norm.cdf(np.abs(z_score)))
        # 得到overall_p_value对应的结果,取p值的中位数
        signal_p_value = np.min(p_values)
        signal_z_score = z_score[np.argmin(p_values)]
        signal_beta = beta[np.argmin(p_values)]
        signal_se = se[np.argmin(p_values)]
        signal_APA = APA_types[np.argmin(p_values)]
        return (
            signal_p_value,
            signal_z_score,
            signal_beta,
            signal_se,
            signal_APA,
            list(p_values),
        )
    else:
        return None, None, None, None, None, None


def process_all_data_Multidimensional(df):
    results = df.apply(
        lambda row: analyze_snp_methylation_bayes_Multidimensional(
            row["A1_APA_count"], row["A2_APA_count"], row["APA_id"]
        ),
        axis=1,
    )
    (
        df["signal_pvalue"],
        df["signal_z_score"],
        df["signal_beta"],
        df["signal_SE"],
        df["signal_APA_type"],
        df["All_pvalues"],
    ) = zip(*results)
    df = df[(df["signal_pvalue"].notna())]
    if len(df) == 0:
        return None
    df["FDR"] = multipletests(df["signal_pvalue"], method="fdr_bh")[1]
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Nanopore direct RNA data call 3'aQTL."
    )
    parser.add_argument("-b", "--bam", type=str, help="bam file path")
    parser.add_argument("--snp_info", type=str, help="processed snp pkl file path")
    parser.add_argument("-o", "--outdirpre", type=str, help="outdir and pre")
    parser.add_argument(
        "-f", "--read_overlap_file", type=str, help="read2apadb overlap file path"
    )
    parser.add_argument("-c", "--chrom", type=str, help="chromosome")
    parser.add_argument("-s", "--strand", type=str, help="different strand processing")
    parser.add_argument("--geno_size", type=str, help="genome size file path")
    parser.add_argument(
        "-t", "--threads", type=int, default=4, help="threads number (default: 4)"
    )
    parser.add_argument(
        "--base_minQ", type=int, default=5, help="base min qscore(default=5)"
    )
    parser.add_argument(
        "--read_minQ", type=int, default=0, help="read min qscore(default=0)"
    )
    args = parser.parse_args()

    output_path = f"{args.outdirpre}_haplotype_{args.chrom}_{args.strand}_tmp.csv"

    base_minQ = args.base_minQ - 1 if args.base_minQ != 0 else 0
    read_minQ = args.read_minQ - 1 if args.read_minQ != 0 else 0

    overlap_df = pd.read_csv(
        args.read_overlap_file,
        sep="\t",
        header=None,
        names=["readid", "polyA_gene"],
        usecols=[3, 9],
    )
    df_dict = dict(
        zip(overlap_df["readid"], overlap_df["polyA_gene"])
    )  # readid:polyA_gene

    snp_info = pd.read_csv(args.snp_info, sep="\t")
    snp_info.columns = ["chrom", "pos1", "rsID", "A1", "A2", "EAF"]
    snp_info = snp_info[snp_info["chrom"] == args.chrom]
    snp_info["pos0"] = snp_info["pos1"].astype(int) - 1
    snp_info["v"] = (
        snp_info["rsID"]
        + ";"
        + snp_info["A1"]
        + ";"
        + snp_info["A2"]
        + ";"
        + snp_info["EAF"].astype(str)
    )
    snp_dict = dict(zip(snp_info["pos0"], snp_info["v"]))  # pos0: rsid;A1;A2;EAF

    geno_size_df = pd.read_csv(
        args.geno_size, sep="\t", header=None, names=["chrom", "size"]
    )
    geno_size_dict = dict(zip(geno_size_df["chrom"], geno_size_df["size"]))
    haplotype_df = pd.DataFrame(
        columns=[
            "chrom",
            "strand",
            "snp_pos_1base",
            "rsID",
            "A1",
            "A2",
            "EAF",
            "A1_APA_count",
            "A2_APA_count",
            "APA_id",
        ]
    )
    end = geno_size_dict[args.chrom]
    haplotype_df = pd.concat(
        [
            haplotype_df,
            count_haplotype(
                args.chrom,
                end,
                args.strand,
                args.bam,
                df_dict,
                snp_dict,
                args.threads,
                base_minQ,
                read_minQ,
            ),
        ],
        ignore_index=True,
    )
    if len(haplotype_df) != 0:
        df = haplotype_df.sort_values(by=["chrom", "snp_pos_1base"])
        # df = process_all_data_fourfold(df)
        df = process_all_data_Multidimensional(df)
        if df is None:
            print(
                f"{args.chrom} {args.strand}中的SNP没有足够的read覆盖度，无法进行统计"
            )
            exit()
        df = df.reset_index(drop=True)
        df.to_csv(output_path, index=None)
        print(f"{output_path}已保存")
