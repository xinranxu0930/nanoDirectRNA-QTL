gtf_file = '/mnt/hpc/home/xuxinran/REF/hg19/gencode.v46lift37.annotation.gtf.gz'
res_dir = '/mnt/hpc/home/xuxinran/DirectSeq/data/zhaolin_240206/240201-zhaolin-RNA-merge/v0.7.2/stability'
full_bam2bed = "/mnt/hpc/home/xuxinran/DirectSeq/data/zhaolin_240206/240201-zhaolin-RNA-merge/v0.6.0/min_qscore_0/stability/nano_merge_calls_sort_map_bam2bed.bed"

import pandas as pd
from subprocess import call
import pickle
import glob
import os

gtf = pd.read_csv(gtf_file,sep="\t",comment="#",header=None,usecols=[0,2,3,4,6,8],names=['chr','type','start','end','strand','info'])

# 所有的gene转bed
gene_gtf = gtf[gtf['type']=='gene']
gene_gtf.loc[:, 'geneID'] = gene_gtf['info'].str.extract(r'gene_id "(.*?)";')
gene_gtf.loc[:, 'start'] = gene_gtf['start'] - 1
gene_gtf.loc[:, 'score'] = 0
gene_gtf = gene_gtf[['chr', 'start', 'end', 'geneID', 'score', 'strand']]
gene_gtf = gene_gtf[~gene_gtf['chr'].isin(['chrX', 'chrY', 'chrM'])]
gene_gtf.to_csv(f'{res_dir}/gene.bed',sep='\t',index=False,header=False)
call(f'bedtools merge -i {res_dir}/gene.bed -s -c 4,5,6 -o collapse > {res_dir}/gene_merge.bed',shell=True)
gene_gtf = pd.read_csv(f'{res_dir}/gene_merge.bed',sep="\t",header=None,names=['chr','start','end','geneID','score','strand'])
gene_gtf.loc[:, 'score'] = 0
gene_gtf.loc[:, 'strand'] = gene_gtf['strand'].astype(str).apply(lambda x: x.split(',')[0] if ',' in x else x)
gene_gtf.to_csv(f'{res_dir}/gene_merge.bed',sep='\t',index=False,header=False)
call(f'rm {res_dir}/gene.bed',shell=True)

# 所有的exon转bed
exon_gtf = gtf[gtf['type']=='exon']
exon_gtf.loc[:, 'geneID'] = exon_gtf['info'].str.extract(r'gene_id "(.*?)";')
exon_gtf.loc[:, 'start'] = exon_gtf['start'] - 1
exon_gtf.loc[:, 'score'] = 0
exon_gtf = exon_gtf[['chr', 'start', 'end', 'geneID', 'score', 'strand']]
exon_gtf = exon_gtf[~exon_gtf['chr'].isin(['chrX', 'chrY', 'chrM'])]
exon_gtf.to_csv(f'{res_dir}/exon.bed',sep='\t',index=False,header=False)
call(f'sort -k 1,1 -k2,2n {res_dir}/exon.bed > {res_dir}/exon_sorted.bed',shell=True)
call(f'bedtools merge -i {res_dir}/exon_sorted.bed -s -c 4,5,6 -o collapse > {res_dir}/exon_merge.bed',shell=True)
exon_merged_df = pd.read_csv(f'{res_dir}/exon_merge.bed', sep='\t', header=None,names=['chr', 'start', 'end', 'geneID', 'score', 'strand'])
exon_merged_df.loc[:, 'score'] = 0
exon_merged_df['strand'] = exon_merged_df['strand'].astype(str).apply(lambda x: x.split(',')[0] if ',' in x else x)
exon_merged_df.to_csv(f'{res_dir}/exon_merge.bed',sep='\t',index=False,header=False)
call(f'rm {res_dir}/exon.bed {res_dir}/exon_sorted.bed',shell=True)

# 获取intron
call(f'bedtools subtract -a {res_dir}/gene_merge.bed -b {res_dir}/exon_merge.bed -s > {res_dir}/intron.bed',shell=True)

# 获取没有intron的gene
call(f'bedtools intersect -a {res_dir}/gene_merge.bed -b {res_dir}/exon_merge.bed -s -wa -wb > {res_dir}/fully_overlapped_genes.bed',shell=True)
df = pd.read_csv(f"{res_dir}/fully_overlapped_genes.bed", sep='\t', header=None,usecols=[0,1,2,3,4,5,7,8],names=['chr','start','end','gene','score','strand','s_t','e_t'])
df_no_intron = df[(df['start']==df['s_t']) & (df['end']==df['e_t'])]
# df_intron = df[(df['start']!=df['s_t']) | (df['end']!=df['e_t'])]
# del df_no_intron['s_t'],df_no_intron['e_t'],df_intron['s_t'],df_intron['e_t']
del df_no_intron['s_t'],df_no_intron['e_t']
df_no_intron.to_csv(f"{res_dir}/genes_no_intron.bed", sep='\t', index=False, header=False)
# df_intron.to_csv(f"{res_dir}/genes_intron.bed", sep='\t', index=False, header=False)
call(f'rm {res_dir}/fully_overlapped_genes.bed',shell=True)

# 获取read和gene的比对情况（没比对上的、比对到no intron上的、有效read）,生成一个统计表
# 没比对上的 ！！！这里应该是对{bam2bed}做个统计 计算每个readID出现的次数，如果和{res_dir}/read_no_gene.bed次数相同证明是这部分没有比对到gene上，
call(f'bedtools subtract -a {full_bam2bed} -b {res_dir}/gene_merge.bed -A -s > {res_dir}/read_no_gene.bed',shell=True)
read_no_gene = pd.read_csv(f'{res_dir}/read_no_gene.bed',sep='\t',header=None,usecols = [3],names=['readID'])
read_no_gene_set = set(read_no_gene['readID'])
# 比对到no intron gene上的,直接对read ID计数就行
call(f'bedtools intersect -a {full_bam2bed} -b {res_dir}/genes_no_intron.bed -s -wa > {res_dir}/read_gene_no_intron.bed',shell=True)
read_gene_no_intron = pd.read_csv(f'{res_dir}/read_gene_no_intron.bed',sep='\t',header=None,usecols = [3],names=['readID'])
read_gene_no_intron_set = set(read_gene_no_intron['readID'])
print(f"read on an intronless gene:{len(read_gene_no_intron_set)}")
print(f"read not in gene:{len(read_no_gene_set)}")

valid_read = pd.read_csv(full_bam2bed,sep='\t',header=None,usecols = [3],names=['readID'])
valid_read_set = set(valid_read['readID'])
valid_read_set  = valid_read_set - read_no_gene_set - read_gene_no_intron_set
print(f"valid read:{len(valid_read_set)}")

# 保存valid_read_set为txt
with open(f'{res_dir}/valid_readID.txt','w') as f:
    for readID in valid_read_set:
        f.write(readID+'\n')

call(f'rm {full_bam2bed} {res_dir}/read_gene_no_intron.bed {res_dir}/read_no_gene.bed {res_dir}/genes_no_intron.bed',shell=True)