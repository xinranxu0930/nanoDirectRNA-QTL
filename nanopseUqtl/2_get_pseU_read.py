import pandas as pd
import pickle
import argparse
import numpy as np
from subprocess import run,call

def get_cov(chrom,strand,bam_file, df, pre,min_qscore):
    if "motif" in df.columns:
        df = df.groupby(['chrom', 'pos_0base', 'strand']).agg({
        'read_id': ';'.join,
        'motif': lambda x: x.value_counts().idxmax()
        }).reset_index()
    else:
        df = df.groupby(['chrom', 'pos_0base', 'strand']).agg({'read_id': ';'.join}).reset_index()
    df['mod_num'] = df['read_id'].str.count(";") + 1
    start_region = min(df['pos_0base'])-10
    end_region = max(df['pos_0base'])+10
    call(f'samtools view -b {bam_file} "{chrom}:{start_region}-{end_region}" > {pre}_{chrom}:{start_region}-{end_region}_{strand}.bam',shell=True)
    call(f'samtools mpileup {pre}_{chrom}:{start_region}-{end_region}_{strand}.bam -d 1000000 -Q {min_qscore} > {pre}_{chrom}:{start_region}-{end_region}_{strand}.txt',shell=True)
    cov_df = pd.read_csv(f'{pre}_{chrom}:{start_region}-{end_region}_{strand}.txt',header=None,sep="\t",usecols=[1,3],names=['pos_1base','cov'])
    cov_df['pos_0base'] = cov_df['pos_1base'] - 1
    df['cov'] = df['pos_0base'].map(cov_df.set_index('pos_0base')['cov']).astype(int)
    df['mod_rate'] = np.where(df['cov'] == 0, 0, (df['mod_num'] / df['cov'] * 100).round(2))
    call(f'rm {pre}_{chrom}:{start_region}-{end_region}_{strand}.bam {pre}_{chrom}:{start_region}-{end_region}_{strand}.txt',shell=True)
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", type=str, help="bam file path")
    parser.add_argument("-f", "--pileup_file", type=str, help="pileup result path")
    parser.add_argument("-o", "--output_prefix", type=str, help="outdir and pre")
    parser.add_argument("-c", "--chrom", type=str, help="different chromosome processing")
    parser.add_argument("-s", "--strand", type=str, help="different strand processing")
    parser.add_argument("-q", "--min_qscore", type=int, default=5, help="pseU min qscore(default=5)")
    args = parser.parse_args()

    min_qscore = args.min_qscore -1 if args.min_qscore > 0 else 0

    df_cols = ['chrom', 'pos_0base', 'strand', 'read_id']
    pileup_df_cols = ['chrom', 'pos_1base', 'strand', 'mod_num', 'cov', 'mod_rate']

    if "motif" in args.pileup_file:
        df_cols.append("motif")
        pileup_df_cols.append("motif")

    pileup_file_res = f"{args.output_prefix}_pseU_sites_{args.chrom}_{args.strand}_tmp.csv"
    pileup_file_dict = f"{args.output_prefix}_pseU_read_{args.chrom}_{args.strand}_tmp.pkl"
    # Using samtools to calculate the coverage of pseU loci in a specified chromosome
    command = f"awk '$1==\"{args.chrom}\"' {args.pileup_file} > {args.output_prefix}_{args.chrom}_{args.strand}_tmp.csv"
    run(command, shell=True, check=True)
    result = run(["wc", "-l", f'{args.output_prefix}_{args.chrom}_{args.strand}_tmp.csv'], text=True, capture_output=True)
    line_count = int(result.stdout.split()[0])
    if line_count == 0:
        print(f"{args.chrom} {args.strand}没有pseU位点")
        exit()
    df = pd.read_csv(f'{args.output_prefix}_{args.chrom}_{args.strand}_tmp.csv', header=None, sep="\t")
    df.columns = df_cols

    pileup_df = get_cov(args.chrom, args.strand, args.bam, df, args.output_prefix, min_qscore)
    pileup_df = pileup_df[pileup_df['mod_rate'] >= 10]
    pileup_df['pos_1base'] = pileup_df['pos_0base'] + 1

    pileup_df_res = pileup_df[pileup_df_cols]
    pileup_df_res.to_csv(pileup_file_res, index=None)
    print(f"{pileup_file_res}已保存")

    pileup_df['k'] = pileup_df['chrom'] + '_' + pileup_df['pos_0base'].astype(str) + '_' + pileup_df['strand']
    pileup_dict = {row['k']: row['read_id'] for index, row in pileup_df.iterrows()}

    with open(pileup_file_dict, 'wb') as file:
        pickle.dump(pileup_dict, file)
    print(f"{pileup_file_dict}已保存")

    call(f"rm {args.output_prefix}_{args.chrom}_{args.strand}_tmp.csv", shell=True)