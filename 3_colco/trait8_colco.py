import pandas as pd
from subprocess import call
import itertools
import numpy as np
from rpy2 import robjects
from rpy2.robjects import pandas2ri
import argparse


def read_file(chrom, strand, full_summury, trait, trait_col=''):
    df = pd.read_csv(full_summury)
    df = df[(df['chrom'] == chrom)&(df['strand'] == strand)]
    df = df[df['p_value']<1]
    df = df.rename(columns={'p_value': f'p_value_{trait}'})
    df['snp_pos_1base'] = df['snp_pos_1base'].astype(str)
    common_col = ['rsID', 'chrom', 'snp_pos_1base', 'strand', 'EAF', f'p_value_{trait}']
    if trait_col != '':
        common_col.append(trait_col)
    df = df[common_col]
    return df


## 现在输入的是1个gene所有的share QTL，两两trait之间准备运行coloc
def gene_2trait_combiner(geneID, gene_df, min_qtl_count, p_threshold, coloc_threshold, p1, p2, p12):
    trait_columns = ['tss', 'promoter', 'polyA', 'isoform', 'stability', 'm6A', 'pseU', 'ASE']
    for col in trait_columns:
        gene_df[col] = gene_df[f'p_value_{col}'].notna()
    ## 下面用来生成所有两两trait组合的df
    sub_dfs = {}
    for trait1, trait2 in itertools.combinations(trait_columns, 2):
        sub_df = gene_df[(gene_df[trait1] != False) & (gene_df[trait2] != False)]
        if len(sub_df['snp_pos_1base']) >= min_qtl_count :
            sub_dfs[f"{trait1}_{trait2}"] = sub_df # 筛选2：两两trait组合的SNP数量大于阈值

    if len(sub_dfs) == 0:
        print(f'In {geneID},all trait combinations has less than {min_qtl_count} shared QTLs')
        return None
    else:
        gene_coloc_res = pd.DataFrame(columns=[
            "chrom","snp_pos_1base","rsID","strand","EAF",
            "p_value_1","p_value_2"
            "PPH4","trait1","trait2"
            ])
        for trait1_trait2, sub_df in sub_dfs.items():
            trait1, trait2 = trait1_trait2.split('_')
            sub_df = sub_df[["chrom" , "snp_pos_1base", "rsID", "strand", "EAF", f"p_value_{trait1}", f"p_value_{trait2}"]]
            sub_df = sub_df.drop_duplicates()
            if not {"m6A", "pseU"} & {trait1, trait2}:
                coloc_res = R_coloc_prepare(sub_df, trait1, trait2, p_threshold, coloc_threshold, p1, p2, p12)
            else:
                coloc_res = R_coloc_prepare(sub_df, trait1, trait2, p_threshold, coloc_threshold, p1, p2, p12)
            if coloc_res:
                gene_coloc_res = pd.concat([gene_coloc_res, coloc_res], ignore_index=True)
        if len(gene_coloc_res)==0:
            return None
        else:
            gene_coloc_res['geneID'] = geneID
            return gene_coloc_res # 返回chrom	snp_pos_1base	rsID	strand	EAF p_value_1	p_value_2 PPH4 trait1 trait2 geneID


## 输入的是某两个trait之间的share QTL，处理一下重复的QTL，再运行colco
def R_coloc_prepare(sub_df, trait1, trait2, p_threshold, coloc_threshold, p1, p2, p12):
    if (sub_df[f"p_value_{trait1}"] > p_threshold).all() and (sub_df[f"p_value_{trait2}"] > p_threshold).all():
        print(f"Both {trait1} and {trait2} are not significant in {geneID}")
        return None
    else:
        sub_df['idx'] = list(range(1,len(sub_df)+1))
        sub_df['snp_pos_1base'] = sub_df['snp_pos_1base'].astype(str)+"_"+sub_df['idx'].astype(str)
        sub_df['rsID'] = sub_df['rsID']+"_"+sub_df['idx'].astype(str)
        sub_df = sub_df.drop(columns=['idx'])
        new_column_names = {
            f"p_value_{trait1}": "p_value_1",
            f"p_value_{trait2}": "p_value_2"
        }
        sub_df = sub_df.rename(columns=new_column_names)
        coloc_res_df = R_coloc(sub_df, coloc_threshold, p1, p2, p12)

        if coloc_res_df:
            coloc_res_df['snp_pos_1base'] = coloc_res_df['snp_pos_1base'].apply(lambda x: x.split("_")[0])
            coloc_res_df['rsID'] = coloc_res_df['rsID'].apply(lambda x: x.split("_")[0])
            coloc_res_df['trait1'] = trait1
            coloc_res_df['trait2'] = trait2
            return coloc_res_df # 结果返回chrom	snp_pos_1base	rsID	strand	EAF p_value_1	p_value_2 PPH4 trait1 trait2


def R_coloc(df, coloc_threshold, p1=1e-4, p2=1e-4, p12=1e-5):
    pandas2ri.activate()
    robjects.r('library(coloc)')
    r_df = pandas2ri.py2rpy(df) # 将 df 从 Pandas 数据框转换为 R 数据框，并赋值给 r_df
    robjects.globalenv['df'] = r_df # 将 r_df 传递给 R 环境，并在 R 环境中将其命名为 df
    r_code = f"""
    result <- coloc.abf(
        dataset1 = list(snp = df$rsID, pvalues = df$p_value_1, type = "quant", N = 104),
        dataset2 = list(snp = df$rsID, pvalues = df$p_value_2, type = "quant", N = 104),
        MAF = df$EAF,
        p1 = {p1},
        p2 = {p2},
        p12 = {p12}
    )
    """
    robjects.r(r_code)
    coloc_all_result = robjects.r['result']
    coloc_summary = coloc_all_result.rx2('summary')
    print(coloc_summary[1])
    print("\n")
    coloc_results = coloc_all_result.rx2('results')
    # priors = coloc_result.rx2('priors')
    py_df = pd.DataFrame(robjects.conversion.rpy2py(coloc_results))
    coloc_res_df = py_df[py_df['SNP.PP.H4'] > coloc_threshold]
    if len(coloc_res_df) == 0:
        return None
    else:
        coloc_res_df = coloc_res_df.rename(columns={'snp': 'rsID', 'SNP.PP.H4': 'PPH4'})
        merged_df = pd.merge(df, coloc_res_df[['rsID', 'PPH4']], on='rsID', how='left')
        return merged_df # 返回的结果就是所有PPH4大于阈值的行，并且最后一列是PPH4值


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='8 Kinds of Nanopore Direct RNA QTL Colocalization.')
    parser.add_argument("--ASE", type=str, help="ASE full summury statistics file")
    parser.add_argument("--m6A", type=str, help="m6A full summury statistics file")
    parser.add_argument("--pseU", type=str, help="pseU full summury statistics file")
    parser.add_argument("--stability", type=str, help="stability full summury statistics file")
    parser.add_argument("--isoform", type=str, help="isoform full summury statistics file")
    parser.add_argument("--polyA", type=str, help="polyA full summury statistics file")
    parser.add_argument("--promoter", type=str, help="promoter full summury statistics file")
    parser.add_argument("--tss", type=str, help="tss full summury statistics file")
    parser.add_argument("-o","--outdir", type=str, help="output directory")
    parser.add_argument("-g","--gene_bed", type=str, help="gene region bed file")
    parser.add_argument("-c","--chrom", type=str, help="chromosome")
    parser.add_argument("-s","--strand", type=str, help="different strand processing")
    parser.add_argument("--min_qtl_count", type=float, default=25, help="Minimum number of shared SNPs for coloc analysis (default=25)")
    parser.add_argument("--coloc_threshold", type=float, default=0.75, help="coloc PPH4 threshold (default=0.75)")
    parser.add_argument("--p_threshold", type=float, default=0.05, help="Pvalue threshold (default=0.05)")
    parser.add_argument("--coloc_p1", type=float, default=1e-4, help="coloc prior1 values (default=1e-4)")
    parser.add_argument("--coloc_p2", type=float, default=1e-4, help="coloc prior2 values (default=1e-4)")
    parser.add_argument("--coloc_p12", type=float, default=1e-5, help="coloc prior12 values (default=1e-5)")
    args = parser.parse_args()

    ASE_df = read_file(args.chrom, args.strand, args.ASE, "ASE")
    m6A_df = read_file(args.chrom, args.strand, args.m6A, "m6A", "m6A_pos_1base")
    pseU_df = read_file(args.chrom, args.strand, args.pseU, "pseU", "pseU_pos_1base")
    stability_df = read_file(args.chrom, args.strand, args.stability, "stability")
    isoform_df = read_file(args.chrom, args.strand, args.isoform, "isoform")
    polyA_df = read_file(args.chrom, args.strand, args.polyA, "polyA")
    promoter_df = read_file(args.chrom, args.strand, args.promoter, "promoter")
    tss_df = read_file(args.chrom, args.strand, args.tss, "tss")

    ## 合并所有的 DataFrame
    dfs = [ASE_df, stability_df, isoform_df, polyA_df, promoter_df, tss_df, m6A_df, pseU_df]
    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on=['rsID', 'chrom', 'snp_pos_1base', 'strand', 'EAF'], how='outer')
    merged_df['non_null_count'] = merged_df.notnull().sum(axis=1) # 统计每行的非空值数量
    merged_df = merged_df[merged_df['non_null_count'] >= 10]  # 4列键 + 至少2种trait(6列)，保证这个SNP至少是两种QTL
    merged_df = merged_df.drop(columns=['non_null_count'])

    ## gene注释
    gene_df = merged_df[['chrom', 'snp_pos_1base', 'strand', 'rsID']]
    gene_df = gene_df.drop_duplicates()
    gene_df['snp_pos_0base'] = gene_df['snp_pos_1base'].astype(int) - 1
    gene_df['score'] = 0
    gene_df = gene_df[['chrom', 'snp_pos_0base', 'snp_pos_1base', 'rsID', 'score', 'strand']]
    gene_df.to_csv(f'{args.outdir}/rs_{args.chrom}_{args.strand}.bed', sep='\t', index=False, header=False)
    call(f'bedtools intersect -wa -wb -s -a {args.outdir}/rs_{args.chrom}_{args.strand}.bed -b {args.gene_bed} > {args.outdir}/gene_rs_{args.chrom}_{args.strand}.bed',shell=True)
    gene_df = pd.read_csv(f'{args.outdir}/gene_rs_{args.chrom}_{args.strand}.bed',header=None,sep='\t',usecols = [0,2,3,5,9],names=['chrom','snp_pos_1base','rsID','strand','geneID'])
    gene_df['snp_pos_1base'] = gene_df['snp_pos_1base'].astype(str)
    merged_df = pd.merge(merged_df, gene_df[['chrom', 'snp_pos_1base', 'rsID', 'strand', 'geneID']], on=['chrom', 'snp_pos_1base', 'rsID', 'strand'], how='left')
    call(f'rm {args.outdir}/rs_{args.chrom}_{args.strand}.bed {args.outdir}/gene_rs_{args.chrom}_{args.strand}.bed',shell=True)

    ## 分区进行colco
    gene_groups = merged_df.groupby('geneID')
    for gene_group in gene_groups:
        geneID = gene_group[0]
        gene_df = gene_group[1]
        genes_coloc_res = pd.DataFrame(columns=[
            "chrom","snp_pos_1base","rsID","strand","EAF",
            "p_value_1","p_value_2"
            "PPH4","trait1","trait2"
            "geneID"
            ])
        # 筛选1:首先对gene整体筛选 保证gene内至少有min_qtl_count个重叠QTL
        if len(set(gene_df['snp_pos_1base'])) >= args.min_qtl_count:
            gene_coloc_res = gene_2trait_combiner(geneID, gene_df, args.min_qtl_count, args.p_threshold, args.coloc_threshold, args.coloc_p1, args.coloc_p2, args.coloc_p12)
            if gene_coloc_res:
                genes_coloc_res = pd.concat([genes_coloc_res, gene_coloc_res], ignore_index=True)
            else:
                continue
        else:
            # print(f"{geneID} has less than {args.min_qtl_count} shared QTLs, skip")
            continue

    if len(genes_coloc_res) != 0 :
        m6A1 = genes_coloc_res[genes_coloc_res['trait1'] == "m6A"]
        m6A2 = genes_coloc_res[genes_coloc_res['trait2'] == "m6A"]
        pseU1 = genes_coloc_res[genes_coloc_res['trait1'] == "pseU"]
        pseU2 = genes_coloc_res[genes_coloc_res['trait2'] == "pseU"]
        others = genes_coloc_res[
            (genes_coloc_res['trait1'] != "m6A") &
            (genes_coloc_res['trait1'] != "pseU") &
            (genes_coloc_res['trait2'] != "m6A") &
            (genes_coloc_res['trait2'] != "pseU")
        ]
        others = others.assign(m6A_pos_1base='', pseU_pos_1base='')

        m6A_df = m6A_df.rename(columns={'p_value_m6A': 'p_value_1'})
        m6A1 = pd.merge(
            m6A1,
            m6A_df[['chrom', 'snp_pos_1base', 'rsID', 'strand', 'EAF', 'm6A_pos_1base', 'p_value_1']],
            on=['chrom', 'snp_pos_1base', 'rsID', 'strand', 'EAF', 'p_value_1'],
            how='left'
        )
        m6A1['pseU_pos_1base'] = ''

        m6A_df = m6A_df.rename(columns={'p_value_1': 'p_value_2'})
        m6A2 = pd.merge(
            m6A2,
            m6A_df[['chrom', 'snp_pos_1base', 'rsID', 'strand', 'EAF', 'm6A_pos_1base', 'p_value_2']],
            on=['chrom', 'snp_pos_1base', 'rsID', 'strand', 'EAF', 'p_value_2'],
            how='left'
        )
        m6A2['pseU_pos_1base'] = ''


        pseU_df = m6A_df.rename(columns={'p_value_pseU': 'p_value_1'})
        pseU1 = pd.merge(
            pseU1,
            pseU_df[['chrom', 'snp_pos_1base', 'rsID', 'strand', 'EAF', 'pseU_pos_1base', 'p_value_1']],
            on=['chrom', 'snp_pos_1base', 'rsID', 'strand', 'EAF', 'p_value_1'],
            how='left'
        )
        pseU1['m6A_pos_1base'] = ''

        pseU_df = pseU_df.rename(columns={'p_value_1': 'p_value_2'})
        pseU2 = pd.merge(
            pseU2,
            pseU_df[['chrom', 'snp_pos_1base', 'rsID', 'strand', 'EAF', 'pseU_pos_1base', 'p_value_2']],
            on=['chrom', 'snp_pos_1base', 'rsID', 'strand', 'EAF', 'p_value_2'],
            how='left'
        )
        pseU2['m6A_pos_1base'] = ''

        result_df = pd.concat([m6A1, m6A2, pseU1, pseU2, others], ignore_index=True)
        result_df.to_csv(f'{args.outdir}/coloc_{args.chrom}_{args.strand}.csv', index=False)
    else:
        print(f"{args.chrom} {args.strand} has no significant colocalization")

