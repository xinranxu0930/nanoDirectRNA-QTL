input_file = '/mnt/hpc/home/xuxinran/REF/hg19/FANTOM5_CAGEpeak_TSS_human.bed.gz'
pre = '/mnt/hpc/home/xuxinran/DirectSeq/data/zhaolin_240206/240201-zhaolin-RNA-merge/v0.7.2/promoter/nano_merge'
res_dir = '/mnt/hpc/home/xuxinran/DirectSeq/data/zhaolin_240206/240201-zhaolin-RNA-merge/v0.7.2/promoter'
bam0 = "/mnt/hpc/home/xuxinran/DirectSeq/data/zhaolin_240206/240201-zhaolin-RNA-merge/v0.7.2/bam/m6A/nano_merge_map0.bam"
bam16 = "/mnt/hpc/home/xuxinran/DirectSeq/data/zhaolin_240206/240201-zhaolin-RNA-merge/v0.7.2/bam/m6A/nano_merge_map16.bam"
max_distance = 200


import pandas as pd
from subprocess import call
import pysam

def get_TSS(bamfile, strand):
    res_dict = set()
    samfile = pysam.AlignmentFile(bamfile, "rb")
    for read in samfile:
        chrom = read.reference_name
        if chrom in [f'chr{i}' for i in range(1, 23)]:
            if strand == '+':
                s = read.reference_start -1
                e = read.reference_start +1
            else:  # strand == '-'
                s = read.reference_end
                e = read.reference_end +2
            res_dict.add((chrom, s, e, read.query_name, 0, strand))
    return res_dict

# 获取所有gene的promoter
call(f"zcat {input_file} | sed '1d' > {res_dir}/cage.bed", shell=True)
call(f"sort -k1,1 -k2,2n {res_dir}/cage.bed -o {res_dir}/cage_sorted.bed",shell=True)
call(f'bedtools merge -i {res_dir}/cage_sorted.bed -d {max_distance} -s -c 4,5,6 -o collapse > {res_dir}/promoter_merge.bed',shell=True)
df = pd.read_csv(f'{res_dir}/promoter_merge.bed',sep="\t",header=None,names=['chrom','start','end','geneID','score','strand'])
autosomes = [f'chr{i}' for i in range(1, 23)]
df = df[df['chrom'].isin(autosomes)]
df = df.drop_duplicates().sort_values(by=['chrom', 'start'])
df['score'] = 0
df['strand'] = df['strand'].astype(str).apply(lambda x: x.split(',')[0] if ',' in x else x)
df_f = df[df['strand'] == '+'].copy()
df_r = df[df['strand'] == '-'].copy()
df_f['start'] = df_f['start'] - 1000
df_f['end'] = df_f['end'] + 100
df_r['start'] = df_r['start'] -100
df_r['end'] = df_r['end'] + 1000
res_df = pd.concat([df_f, df_r])
res_df[['start', 'end']] = res_df[['start', 'end']].astype(int)
res_df.to_csv(f'{res_dir}/promoter_merge.bed',sep='\t',index=False,header=False)

# 获取read的起始位置
res_dict_0 = get_TSS(bam0, "+")
res_dict_16 = get_TSS(bam16, "-")
merged_set = res_dict_0.union(res_dict_16)
# 将集合内容写入BED文件
with open(f'{pre}_read5.bed', 'w') as file:
    for item in merged_set:
        file.write('\t'.join(map(str, item)) + '\n')

call(f'bedtools intersect -s -a {pre}_read5.bed -b {res_dir}/promoter_merge.bed -wa -wb> {pre}_overlap.bed',shell=True)
call(f"awk '!seen[$4]++' {pre}_overlap.bed > {pre}_overlap_uniq.bed",shell=True)
call(f'rm {res_dir}/cage.bed {res_dir}/cage_sorted.bed {pre}_overlap.bed',shell=True)