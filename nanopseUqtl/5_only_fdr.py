import numpy as np
import os
import glob
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy as np
from subprocess import call
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Organize the results and calculate FDR.')
    parser.add_argument("-o","--outdirpre", type=str, help="outdir and pre")
    args = parser.parse_args()

    output_all_path = f"{args.outdirpre}_pseU_summary.csv"

    files = glob.glob(f'{args.outdirpre}_haplotype_chr*_tmp.csv')
    dfs = [pd.read_csv(file) for file in files]
    df = pd.concat(dfs, ignore_index=True)
    df = df.sort_values(by=['chrom', 'pseU_pos_1base', 'snp_pos_1base'])
    df = df.reset_index(drop=True)

    print("pvalue小于0.05的个数：", len(df[df['p_value']<0.05]))
    print("fdr小于0.1的个数：", len(df[df['FDR']<0.1]))
    print("pseUQTL同时也是pseU的个数： ", len(df[(df['ispseU?']=="Yes")&(df['p_value']<0.05)]))

    df.to_csv(output_all_path, index=False)
    call(f'rm {args.outdirpre}_haplotype_chr*_tmp.csv',shell=True)


