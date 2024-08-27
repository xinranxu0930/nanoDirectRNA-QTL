import os
import glob
import pandas as pd
import pickle
from subprocess import call
import re
from itertools import groupby
import logomaker
import matplotlib.pyplot as plt
from Bio import SeqIO
import argparse

def fasta_iter(fasta_name):
    with open(fasta_name) as filename:
        faiter = (x[1] for x in groupby(filename, lambda line: line[0] == ">"))
        for header in faiter:
            header_str = next(header)[1:].strip()  # 使用next代替__next__，并修改变量名
            seq = "".join(s.strip() for s in next(faiter))  # 使用next代替__next__
            yield header_str, seq  # 使用yield代替return

def get_reverse_complementary_sequence(seq):
    seqreverse = seq[::-1]
    transtable = str.maketrans('ATGCatgc','TACGtacg')
    finalseq = seqreverse.translate(transtable)
    return finalseq

def get_motif(df,s,e,n,fa_file,pre_path):
    df['name'] = df['chrom'] + '_' + df['pos_1base'].astype(str)+ '_' + df['strand']
    df.loc[df['strand'] == '+', 'Start'] = df.loc[df['strand'] == '+', 'pos_1base'] - s - 1
    df.loc[df['strand'] == '+', 'End'] = df.loc[df['strand'] == '+', 'pos_1base'] + e
    df.loc[df['strand'] == '-', 'Start'] = df.loc[df['strand'] == '-', 'pos_1base'] - e - 1
    df.loc[df['strand'] == '-', 'End'] = df.loc[df['strand'] == '-', 'pos_1base'] + s
    df['Start'] = df['Start'].astype(int)
    df['End'] = df['End'].astype(int)
    get_fasta = df[['chrom','Start','End','name']]
    get_fasta.loc[:,'score'] = 0
    get_fasta.loc[:,'strand'] = df['strand']
    get_fasta.to_csv(f"{pre_path}_{n}.bed", header=None, sep="\t",index=False)
    call(f"bedtools getfasta -name -s -fi {fa_file} -bed {pre_path}_{n}.bed -fo {pre_path}_{n}.fa",shell=True)
    return df

def check_base(pre_path,df,n,pos):
    fasta_dict = {}
    with open(f"{pre_path}_{n}.fa", "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fasta_dict[record.id.split("::")[0]] = str(record.seq)
    df['motif_fa'] = df['name'].apply(lambda x: fasta_dict[x])
    df['base'] = df['motif_fa'].apply(lambda x: x[pos])
    df = df[(df['base']=="T")|(df['base']=="t")]
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_prefix", type=str, help="outdir and pre")
    parser.add_argument("--fa", type=str, help="fasta file path")
    args = parser.parse_args()

    res_path = f"{args.output_prefix}_pseU_sites.csv"
    pkl_f_path = f"{args.output_prefix}_pseU_read_tmp_f.pkl"
    pkl_r_path = f"{args.output_prefix}_pseU_read_tmp_r.pkl"

    # merge the results
    files_f = glob.glob(f'{args.output_prefix}_pseU_sites_chr*_+_tmp.csv')
    files_r = glob.glob(f'{args.output_prefix}_pseU_sites_chr*_-_tmp.csv')
    dfs_f = [pd.read_csv(file) for file in files_f]
    dfs_r = [pd.read_csv(file) for file in files_r]

    df_f = pd.concat(dfs_f, ignore_index=True)
    df_r = pd.concat(dfs_r, ignore_index=True)
    df = pd.concat([df_f, df_r], ignore_index=True)
    df = df[(df['motif'] == "AT") & (df['cov'] > 10) & (df['mod_rate'] > 10) | (df['motif'] != "AT")] # df1多的不太靠谱，所以筛选严格一点
    df = df.sort_values(by=['chrom', 'pos_1base'])

    files_r = glob.glob(f'{args.output_prefix}_pseU_read_chr*_-_tmp.pkl')
    files_f = glob.glob(f'{args.output_prefix}_pseU_read_chr*_+_tmp.pkl')
    dicts_r = [pickle.load(open(file, 'rb')) for file in files_r]
    dicts_f = [pickle.load(open(file, 'rb')) for file in files_f]
    merged_dict_r = {k: v for d in dicts_r for k, v in d.items()}
    merged_dict_f = {k: v for d in dicts_f for k, v in d.items()}

    df1 = df[df['motif'].apply(lambda x: len(x) == 3)]
    df7 = df[df['motif'].apply(lambda x: len(x) == 5)]
    df4 = df[df['motif'].apply(lambda x: len(x) == 6)]
    df1 = get_motif(df1,2,0,"pus1",args.fa,args.output_prefix)
    df4 = get_motif(df4,2,3,"pus4",args.fa,args.output_prefix)
    df7 = get_motif(df7,2,2,"pus7",args.fa,args.output_prefix)

    df1 = check_base(args.output_prefix,df1,"pus1",1)
    df4 = check_base(args.output_prefix,df4,"pus4",2)
    df7 = check_base(args.output_prefix,df7,"pus7",2)

    merged_df = pd.concat([df1, df4, df7], ignore_index=True)
    df_save = merged_df[['chrom','pos_1base','strand','mod_num', "cov", 'mod_rate', 'motif']]
    df_save.to_csv(res_path, index=False)

    df_save['pos_0base'] = df_save['pos_1base'] - 1
    df_save['k'] = df_save['chrom'] + '_' + df_save['pos_0base'].astype(str) + '_' + df_save['strand']
    readid_set = set(df_save['k'])
    filtered_merged_dict_r = {key: value for key, value in merged_dict_r.items() if key in readid_set}
    filtered_merged_dict_f = {key: value for key, value in merged_dict_f.items() if key in readid_set}
    with open(pkl_r_path, 'wb') as file:
        pickle.dump(filtered_merged_dict_r, file)
    with open(pkl_f_path, 'wb') as file:
        pickle.dump(filtered_merged_dict_f, file)

