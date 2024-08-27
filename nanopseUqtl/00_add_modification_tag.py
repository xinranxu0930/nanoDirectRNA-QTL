import mappy as mp
import argparse
import pysam

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Add modification tag to bam file")
    parser.add_argument("-f", "--fastq", type=str, help="fastq file path")
    parser.add_argument("-b", "--bam", type=str, help="bam file without modification tag")
    parser.add_argument("-t", "--threads", type=int, default=2, help="threads")
    args = parser.parse_args()

    newbam = args.bam.replace(".bam","_mod_map.bam")

    read_flag = {}
    for i in mp.fastx_read(args.fastq, read_comment=True):
        read_name = i[0]
        MM_tag_value = i[3].split("\t")[0].split("T+17802.,")[1][:-1] if len(i[3].split("\t")[0].split("T+17802.,")) > 1 else None
        ML_tag_value = i[3].split("\t")[1].split("B:C,")[1] if len(i[3].split("\t")[1].split("B:C,")) > 1 else None
        if MM_tag_value is not None and ML_tag_value is not None:
            read_flag[read_name] = [MM_tag_value, ML_tag_value]

    with pysam.AlignmentFile(args.bam, "rb",threads=args.threads) as samfile:
        with pysam.AlignmentFile(newbam, "wb", header=samfile.header) as out_file:
            for read in samfile:
                read_name = read.query_name
                if (read_name in read_flag) and (read.flag == 0 or read.flag == 16):
                    MM_tag_value = read_flag[read_name][0]
                    ML_tag_value = read_flag[read_name][1]
                    ML_tag_value_list = list(map(int, ML_tag_value.split(",")))
                    read.set_tag("MM", f"T+17802?,{MM_tag_value};", replace=False)
                    read.set_tag("ML", ML_tag_value_list, replace=False)
                    out_file.write(read)
