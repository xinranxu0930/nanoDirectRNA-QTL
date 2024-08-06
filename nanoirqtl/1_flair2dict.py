import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--files', nargs='+', help='Input flair results file(s),*.combined.isoform.read.map.txt')
parser.add_argument("-p","--dir_pre", type=str, help="dir pre file path")
args = parser.parse_args()

all_dicts = {}
res_dict_path = f"{args.dir_pre}_isoform_read.pkl"

for file in args.files:
    isoform_dict = {}
    with open(file, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            isoformID = parts[0]
            reads_l = parts[1].split(",")
            for read in reads_l:
                isoform_dict[read] = isoformID

with open(res_dict_path, 'wb') as file:
    pickle.dump(isoform_dict, file)
print(f"{res_dict_path}已保存")