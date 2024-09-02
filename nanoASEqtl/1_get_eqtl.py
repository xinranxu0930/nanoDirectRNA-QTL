import pandas as pd
import pysam
import argparse
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests


def count_haplotype(
    chrom, end, strand, bamfile, snp_file_dict, threadsnum, base_minQ, read_minQ
):
    # snp_file_dict snp位点信息 pos0: rsid;A1;A2;EAF
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
                seqs = [
                    i.upper() for i in pileupcolumn.get_query_sequences()
                ]  # 被覆盖位点处的所有碱基
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
            "A1_num",
            "A2_num",
            "p_value",
            "Zscore",
            "beta",
            "SE",
        ]
    )
    result = get_haplotypes(snp_bases, snp_file_dict)
    if result is None:
        print(f"{chrom} {strand}中没有符合条件的snp")
        return None
    else:
        ith = 0
        for i in result:
            haplotype_df.loc[ith] = (
                chrom,
                strand,
                i[6] + 1,
                i[4],
                i[0],
                i[1],
                i[5],
                i[2],
                i[3],
                i[7],
                i[8],
                i[9],
                i[10],
            )
            ith += 1
        return haplotype_df


def get_haplotypes(snp_bases, snp_file_dict):
    res_l = []
    # snp_bases [pos0]:{read1:"A",read2:"A",read3:"T"}
    for snp_pos, snp_base in snp_bases.items():
        snpID, A1, A2, eaf = snp_file_dict[snp_pos].split(";")[:4]
        snp_base_count = {A1: 0, A2: 0}
        for _, base in snp_base.items():
            if base == A1:
                snp_base_count[A1] += 1
            elif base == A2:
                snp_base_count[A2] += 1
        if snp_base_count == {A1: 0, A2: 0}:
            continue
        p_value, z_score, beta, se = ase_equivalence_bayes(
            snp_base_count[A1], snp_base_count[A2], float(eaf)
        )
        if p_value:
            res_l.append(
                [
                    A1,
                    A2,
                    snp_base_count[A1],
                    snp_base_count[A2],
                    snpID,
                    eaf,
                    snp_pos,
                    p_value,
                    z_score,
                    beta,
                    se,
                ]
            )
    if len(res_l) != 0:
        return res_l
    else:
        return None


def ase_equivalence_bayes(a1_count, a2_count, eaf, n_samples=10000000, random_state=42):
    total_count = a1_count + a2_count
    if a1_count < 5 or a2_count < 5 or total_count < 208:
        return None, None, None, None
    
    prior_alpha = eaf * total_count
    prior_beta = (1 - eaf) * total_count
    post_alpha = prior_alpha + a1_count
    post_beta = prior_beta + a2_count

    rng = np.random.default_rng(random_state)
    theta_samples = stats.beta.rvs(post_alpha, post_beta, size=n_samples, random_state=rng)

    mean_effect_size = np.mean(theta_samples - eaf)
    se_effect_size = np.std(theta_samples - eaf)
    if se_effect_size == 0:
        return None, None, None, None
    z_score = mean_effect_size / se_effect_size
    p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
    return p_value, z_score, mean_effect_size, se_effect_size


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Nanopore direct RNA data call ASE expression qtl."
    )
    parser.add_argument("-b", "--bam", type=str, help="bam file path")
    parser.add_argument("--snp_info", type=str, help="processed snp pkl file path")
    parser.add_argument("-o", "--outdirpre", type=str, help="outdir and pre")
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

    base_minQ = args.base_minQ - 1 if args.base_minQ != 0 else 0
    read_minQ = args.read_minQ - 1 if args.read_minQ != 0 else 0

    output_path = f"{args.outdirpre}_haplotype_{args.chrom}_{args.strand}_tmp.csv"

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
    end = geno_size_dict[args.chrom]
    haplotype_df = pd.DataFrame(
        columns=[
            "chrom",
            "strand",
            "snp_pos_1base",
            "rsID",
            "A1",
            "A2",
            "EAF",
            "A1_num",
            "A2_num",
            "p_value",
            "Zscore",
            "beta",
            "SE",
        ]
    )

    haplotype_df = pd.concat(
        [
            haplotype_df,
            count_haplotype(
                args.chrom,
                end,
                args.strand,
                args.bam,
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
        df["FDR"] = multipletests(df["p_value"], method="fdr_bh")[1]
        if df is None:
            print(
                f"{args.chrom} {args.strand}中的SNP没有足够的read覆盖度，无法进行统计"
            )
            exit()
        df = df.reset_index(drop=True)
        df.to_csv(output_path, index=None)
        print(f"{output_path}已保存")
