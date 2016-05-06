import os
import sys
import csv
import argparse
import numpy as np
import math
import matplotlib.pyplot as plt
from gff3 import Gff3Parser

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-o","--output_file",help="output file")
parser.add_argument("-s","--srna_file",help="output file")
parser.add_argument("-t","--type_", help="output file")
args = parser.parse_args()

def get_srna(name):
    names = []
    fh = open(args.srna_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if row[2] == name:
            srnas = row[-4].split(";")
            for srna in srnas:
                names.append(srna.split("|")[-2])
            break
    return ",".join(names)

def get_name(row):
    datas = row[9].split(";")
    name = ""
    locus = ""
#    if row[3] == "CDS":
    for data in datas:
        name = str(row[4]) + "_" + str(row[5]) + "_" + row[7]
#    else:
#        for data in datas:
#            if data[:4] == "Name":
#                name = data[5:]
#            elif data[:8] == "sRNA_hit":
#                if (data[-2:] == "NA"):
#                    srna = name
#                else:
#                    srna = get_srna(name)
#        if name != srna:
#            name = name + ":" + srna
    return name

def check_hypo(row):
    detect = False
    if row[3] == "CDS":
        for infos in row[9].split(";"):
            if ("product" in infos) and ("hypothetical" not in infos):
                detect = True
    return detect

def main():
    tsbs = []
    pmems = []
    fh = open(args.input_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        detect = False
        if (not row[0].startswith("Orientation")) and (
            row[0] == "sense"):
            if (row[3] == "CDS") and ((args.type_ == "both") or (args.type_ == "CDS")):
                detect = check_hypo(row)
            if (row[3] == "sRNA") and ((args.type_ == "both") or (args.type_ == "sRNA")):
                detect = True
            if detect:
                gene_name = get_name(row)
                name = row
                feature = row[3]
                row = map(float, [row[11], row[13], row[15],
                      row[19], row[21], row[23],
                      row[17], row[25], row[27], row[29],
                      row[33], row[35], row[37], row[31]])
                if max(row) >= 10:
                    row = map(str, row)
                    print(gene_name + "\t" + "\t".join(row))

if __name__ == "__main__":
    main()
