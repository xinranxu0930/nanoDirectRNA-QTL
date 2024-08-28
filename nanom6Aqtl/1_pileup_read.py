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
    methy_base = "A" if strand == "+" else "T"
    haplotype_set = set()
    with pysam.AlignmentFile(map_bam_file, "rb", threads=threads) as bamfile:
        for read in bamfile:
            read_id, read_chr, read_seq, read_qscore = (
                read.query_name,
                read.reference_name,
                read.query_sequence,
                read.query_qualities,
            )
            MM_a_tag = read.modified_bases[list(read.modified_bases.keys())[0]]
            if len(read.modified_bases) != 1 or len(MM_a_tag) == 0:
                continue
            for tup in MM_a_tag:
                if (
                    (read_qscore[tup[0]] >= min_qscore)
                    and (read_seq[tup[0]] == methy_base)
                    and (tup[1] >= (mod_threshold * 256))
                ):
                    target_index = tup[0]
                    geno_pos = find_position(
                        read.cigartuples, read.reference_start, target_index)
                    if geno_pos:
                        haplotype_set.add(
                            f"{read_chr}\t{geno_pos}\t{strand}\t{read_id}")
    return haplotype_set

def get_read_modify_motif(map_bam_file,mod_threshold,min_qscore,threads,strand):
    motif_pattern = re.compile(r"[GAT][GA]AC[ATC]")
    haplotype_set = set()
    with pysam.AlignmentFile(map_bam_file, "rb", threads=threads) as bamfile:
        for read in bamfile:
            read_id, read_chr, read_seq, read_qscore = (
                read.query_name,
                read.reference_name,
                read.query_sequence,
                read.query_qualities,
            )
            MM_a_tag = read.modified_bases.get(list(read.modified_bases.keys())[0], [])
            if len(read.modified_bases) != 1 or len(MM_a_tag) == 0:
                continue
            for tup in MM_a_tag:
                if read_qscore[tup[0]] >= min_qscore:
                    motif_region = read_seq[tup[0] - 2 : tup[0] + 3] if strand == "+" else get_reverse_complementary_sequence(read_seq[tup[0] - 2 : tup[0] + 3])

                    match = motif_pattern.match(motif_region)
                if not match:
                    continue
                motif = match.group()
                if tup[1] >= (mod_threshold*256):
                    target_index = tup[0]
                    geno_pos = find_position(read.cigartuples, read.reference_start, target_index)
                    if geno_pos:
                        haplotype_set.add(f"{read_chr}\t{geno_pos}\t{strand}\t{read_id}\t{motif}")
    return haplotype_set


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o","--output_prefix", type=str, help="output folder and file prefix")
    parser.add_argument("-b","--bam", type=str, help="bam file path")
    parser.add_argument("-m","--mod_threshold", type=float, default=0.8 ,help="m6A mod threshold(default=0.8)")
    parser.add_argument("-t","--threads", type=int, default=10, help="threads")
    parser.add_argument("-q","--min_qscore", type=int, default=5, help="base min qscore(default=5)")
    parser.add_argument("-s","--strand", type=str, help="strand")
    parser.add_argument("--motif", help="motif context", action="store_true")
    args = parser.parse_args()

    if args.motif:
        read_m6A_set = get_read_modify_motif(
            args.bam,args.mod_threshold,args.min_qscore,args.threads,args.strand)
        n = "_motif"
    else:
        read_m6A_set = get_read_modify_base(
            args.bam,args.mod_threshold,args.min_qscore,args.threads,args.strand)
        n = ""

    read_m6A_file = f"{args.output_prefix}_read_m6A_pos_{'f' if args.strand == '+' else 'r'}{n}_tmp.txt"

    with open(read_m6A_file, 'w') as file:
        # 遍历集合，将每个元素转换为字符串并写入文件
        for item in read_m6A_set:
            file.write(f"{item}\n")
