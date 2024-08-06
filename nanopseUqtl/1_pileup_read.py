import pysam
import pandas as pd
import re
import argparse

def find_position(cigar_tuples, read_start, target_index):
    # op:M0 I1 N3 D2 S4
    current_index = 0
    current_position = read_start
    for cigar in cigar_tuples:
        op, length = cigar
        if op in [0,1,4]:
            current_index += length
        if current_index <= target_index:
            if op in [0,2,3]:
                current_position += length
        else:
            if op == 0:
                return current_position + (target_index - (current_index-length))
            else:
                return None

def get_reverse_complementary_sequence(seq):
    seqreverse = seq[::-1]
    transtable = str.maketrans('ATGCatgc','TACGtacg')
    finalseq = seqreverse.translate(transtable).upper()
    return finalseq

def get_read_modify_base(map_bam_file,mod_threshold,min_qscore,threads,strand):
    mod_base = "T" if strand == "+" else "A"
    haplotype_set = set()
    with pysam.AlignmentFile(map_bam_file, "rb", threads=threads) as bamfile:
        for read in bamfile:
            read_id,read_chr,read_seq,read_qscore = read.query_name,read.reference_name,read.query_sequence,read.query_qualities
            if len(read.modified_bases) != 1:
                continue
            MM_u_tag = read.modified_bases[list(read.modified_bases.keys())[0]]
            if len(MM_u_tag) == 0:
                continue
            for tup in MM_u_tag:
                if read_qscore[tup[0]] < min_qscore:
                    continue
                if read_seq[tup[0]] != mod_base:
                    continue
                if tup[1] < (mod_threshold*256):
                    continue
                else:
                    target_index = tup[0]
                    geno_pos = find_position(read.cigartuples, read.reference_start, target_index)
                    if geno_pos is not None:
                        haplotype_set.add(f"{read_chr}\t{geno_pos}\t{strand}\t{read_id}")
    return haplotype_set

def get_read_modify_motif(map_bam_file,mod_threshold,min_qscore,threads,strand):
    pus1_motif_pattern = re.compile(r"AT")
    pus4_motif_pattern = re.compile(r"[AG]GTTC[ATCG]A[ATCG][CT]C[CT]")
    pus7_motif_pattern = re.compile(r"TGTA[AG]")
    haplotype_set = set()
    with pysam.AlignmentFile(map_bam_file, "rb", threads=threads) as bamfile:
        for read in bamfile:
            read_id,read_chr,read_seq,read_qscore = read.query_name,read.reference_name,read.query_sequence,read.query_qualities
            if len(read.modified_bases) != 1:
                continue
            MM_u_tag = read.modified_bases[list(read.modified_bases.keys())[0]]
            if len(MM_u_tag) == 0:
                continue
            for tup in MM_u_tag:
                if read_qscore[tup[0]] < min_qscore:
                    continue

                pus1_region = read_seq[tup[0]-1:tup[0]+1] if strand == "+" else get_reverse_complementary_sequence(read_seq[tup[0]:tup[0]+2])
                pus4_region = read_seq[tup[0]-3:tup[0]+8] if strand == "+" else get_reverse_complementary_sequence(read_seq[tup[0]-7:tup[0]+4])
                pus7_region = read_seq[tup[0]-2:tup[0]+3] if strand == "+" else get_reverse_complementary_sequence(read_seq[tup[0]-2:tup[0]+3])
                match1 = pus1_motif_pattern.match(pus1_region)
                match4 = pus4_motif_pattern.match(pus4_region)
                match7 = pus7_motif_pattern.match(pus7_region)
                if not match1 and not match4 and not match7:
                    continue
                elif match1:
                    motif = match1.group()
                elif match4:
                    motif = match4.group()
                elif match7:
                    motif = match7.group()

                if tup[1] < (mod_threshold*256):
                    continue
                else:
                    target_index = tup[0]
                    geno_pos = find_position(read.cigartuples, read.reference_start, target_index)
                    if geno_pos is not None:
                        haplotype_set.add(f"{read_chr}\t{geno_pos}\t{strand}\t{read_id}\t{motif}")
    return haplotype_set

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--dir_pre", type=str, help="dir pre file path")
    parser.add_argument("-b","--bam", type=str, help="bam file path")
    parser.add_argument("-m","--mod_threshold", type=float, default=0.7 ,help="pseU mod threshold(default=0.7)")
    parser.add_argument("-t","--threads", type=int, default=10, help="threads")
    parser.add_argument("-q","--min_qscore", type=int, default=5, help="base min qscore(default=5)")
    parser.add_argument("-s","--strand", type=str, help="strand")
    parser.add_argument("--motif", help="motif context", action="store_true")
    args = parser.parse_args()

    if args.motif:
        read_pseU_set = get_read_modify_motif(args.bam,args.mod_threshold,args.min_qscore,args.threads,args.strand)
        n = "_motif"
    else:
        read_pseU_set = get_read_modify_base(args.bam,args.mod_threshold,args.min_qscore,args.threads,args.strand)
        n = ""
    read_pseU_file = f"{args.dir_pre}_read_pseU_pos_{'f' if args.strand == '+' else 'r'}{n}_tmp.txt"

    with open(read_pseU_file, 'w') as file:
        # 遍历集合，将每个元素转换为字符串并写入文件
        for item in read_pseU_set:
            file.write(f"{item}\n")
