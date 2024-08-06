import argparse
import pysam

def res_0(bam_file,threads):
    with pysam.AlignmentFile(bam_file, "rb",threads=threads) as samfile:
        with pysam.AlignmentFile(map_bam_file0, "wb", header=samfile.header) as out_file:
            for read in samfile:
                if read.flag == 0:
                    out_file.write(read)
def res_16(bam_file,threads):
    with pysam.AlignmentFile(bam_file, "rb",threads=threads) as samfile:
        with pysam.AlignmentFile(map_bam_file16, "wb", header=samfile.header) as out_file:
            for read in samfile:
                if read.flag == 16:
                    out_file.write(read)
def res_all(bam_file,threads):
    with pysam.AlignmentFile(bam_file, "rb",threads=threads) as samfile:
        with pysam.AlignmentFile(map_bam_file, "wb", header=samfile.header) as out_file:
            for read in samfile:
                if read.flag == 0 or read.flag == 16:
                    out_file.write(read)
                    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--basecall_bam", type=str, help="dorado basecall bam")
    parser.add_argument("-p", "--dir_pre", type=str, help="dir pre file path")
    parser.add_argument("-t", "--threads", type=int, default=2, help="threads")
    parser.add_argument("-m", "--mode", choices=['all', 'strand', 'map'], default='all', help="mode to run: all or strand or map")
    args = parser.parse_args()

    map_bam_file0 = f"{args.dir_pre}_calls_sorted_map0.bam"
    map_bam_file16 = f"{args.dir_pre}_calls_sorted_map16.bam"
    map_bam_file = f"{args.dir_pre}_calls_sorted_map.bam"

    if args.mode == 'all':
        res_16(args.basecall_bam,args.threads)
        res_0(args.basecall_bam,args.threads)
        res_all(args.basecall_bam,args.threads)
    elif args.mode == 'strand':
        res_16(args.basecall_bam,args.threads)
        res_0(args.basecall_bam,args.threads)
    else:
        res_all(args.basecall_bam,args.threads)

