import pandas as pd

res_pre = "/mnt/hpc/home/xuxinran/DirectSeq/data/zhaolin_240206/240201-zhaolin-RNA-merge/v0.7.2/stability/nano_merge"
bam2bed = "/mnt/hpc/home/xuxinran/DirectSeq/data/zhaolin_240206/240201-zhaolin-RNA-merge/v0.7.2/stability/nano_merge_calls_sort_map_valid_bam2bed.bed"

exon_df = pd.read_csv(f'{res_pre}_exon.bed', sep="\t",header=None,usecols=[1,2,3],names=["s","e","readID"])
exon_df.loc[:, 'exon_overlap'] = exon_df['e']-exon_df['s']
exon_df = exon_df.groupby('readID')['exon_overlap'].sum().reset_index(name='exon_len')

intron_df = pd.read_csv(f'{res_pre}_intron.bed', sep="\t",header=None,usecols=[1,2,3],names=["s","e","readID"])
intron_df.loc[:, 'intron_overlap'] = intron_df['e']-intron_df['s']
intron_df = intron_df.groupby('readID')['intron_overlap'].sum().reset_index(name='intron_len')

df = pd.read_csv(bam2bed, sep='\t', header=None, usecols=[3], names=['readID'])
df = df.drop_duplicates(keep='first')
df = df.merge(intron_df[['readID', 'intron_len']], on=['readID'], how='left')
df = df.merge(exon_df[['readID', 'exon_len']], on=['readID'], how='left')
df[['exon_len', 'intron_len']] = df[['exon_len', 'intron_len']].fillna(0)
df = df[df['exon'] != 0]
df.loc[:, 'len'] = df['intron_len']+df['exon_len']
df['stability'] = df.apply(lambda x: x['exon_len']/x['len'], axis=1)
mean_stability = df['stability'].mean()
std_stability = df['stability'].std()
df['stability_zscore'] = (df['stability'] - mean_stability) / std_stability # 计算标准化的RNA稳定性得分（Z-score）

df.to_csv(f'{res_pre}_read_stability.csv', index=False)