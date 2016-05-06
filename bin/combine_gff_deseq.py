import os
import sys
import csv
import argparse
from gff3 import Gff3Parser

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-g","--gene_file",help="input file")
parser.add_argument("-t","--type_",help="input file")
args = parser.parse_args()

def main():
    names = []
    if (args.type_ == "deseq") or (args.type_ == "fold_first"):
        gh = open(args.gene_file, "r")
        for row in csv.reader(gh, delimiter='\t'):
            if "Orientation" not in row[0]:
                if row[0] == "sense":
                    datas = row[9].split(";")
                    for data in datas:
                        if "Name" in data:
                            name = data.split("=")[-1]
                    names.append({"info": "\t".join(row[:10]), "start": int(row[4]), "end": int(row[5]), "strand": row[7], "name": name})
                    
    else:
        gh = open(args.gene_file, "r")
        for row in csv.reader(gh, delimiter='\t'):
            datas = row[9].split(";")
            for data in datas:
                if "Name" in data:
                    name = data.split("=")[-1]
            names.append({"info": "\t".join(row), "start": int(row[4]), "end": int(row[5]), "strand": row[7], "name": name})
    fh = open(args.input_file, "r")
    out = open("tmp_fold_change.csv", "w")
    for row in csv.reader(fh, delimiter='\t'):
#        if (row[0] != "baseMean") and ("sRNA" not in row[0]) and (
#            "_" in row[0]) and (("+" in row[0]) or ("-" in row[0])):
        if (row[0] != "baseMean"):
            info = row[0].split("_")
            for name in names:
                if (name["start"] == int(info[0])) and (
                    name["end"] == int(info[1])) and (
                    name["strand"] == info[2]):
                    if args.type_ == "deseq":
                        print(name["info"] + "\t" + "\t".join(row[1:]))
                    elif args.type_ == "fold_first":
                        out.write(name["info"] + "\t" + "0" + "\t")
                        out.write(row[2] + "\n")
                    elif args.type_ == "fold_middle":
                        out.write(name["info"] + "\t" + row[2] + "\n")
                    elif args.type_ == "fold_last":
                        out.write(name["info"] + "\t" + row[2] + "\n")
#        elif (row[0] != "baseMean"):
#            for name in names:
#                if (name["name"] in row[0]):
#                    if args.type_ == "deseq":
#                        print(name["info"] + "\t" + "\t".join(row[1:]))
#                    elif args.type_ == "fold_first":
#                        out.write(name["info"] + "\t" + "0" + "\t")
#                        out.write(row[2] + "\n")
#                    elif args.type_ == "fold_middle":
#                        out.write(name["info"] + "\t" + row[2] + "\n")
#                    elif args.type_ == "fold_last":
#                        out.write(name["info"] + "\t" + row[2] + "\n")

if __name__ == "__main__":
    main()
