#!/usr/bin/python

import os
import sys
import csv
import argparse
import numpy as np
import math
import matplotlib as pylab
import matplotlib.pyplot as plt

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-s","--rank_file",help="input file")
parser.add_argument("-r","--rpkm_file",help="output file")
parser.add_argument("-l","--log", default=False, action="store_true", help="output file")
parser.add_argument("-o","--out_file",help="output file")
args = parser.parse_args()

def get_gene(row):
    genes = []
    for rpkm in row[1:]:
        genes.append(math.log((1 + float(rpkm)), 2))
    return genes

def main():
    fh = open(args.rank_file, "r")
    highs = []
    for row in csv.reader(fh, delimiter='\t'):
        all_rows = map(float, row[1:-1])
        if ((int(row[-1]) <= 150) and (max(all_rows) <= 50)):
#        if (int(row[-1]) <= 200):
            highs.append(row[0])
    rh = open(args.rpkm_file, "r")
    plt.figure(figsize= (12.5, 8))
    x = np.arange(14)
    labels = ["TSB_OD_0.2", "TSB_OD_0.5", "TSB_OD_1", "TSB_t0", "TSB_t1", "TSB_t2", "TSB_ON",
              "pMEM_OD_0.2", "pMEM_OD_0.5", "pMEM_OD_1", "pMEM_t0", "pMEM_t1", "pMEM_t2", "pMEM_ON"]
    num = 1
    for row in csv.reader(rh, delimiter='\t'):
        if args.log:
            genes = get_gene(row)
        else:
            genes = map(float, row[1:])
        if row[0].strip() in highs:
            if ":" in row[0]:
                name = "ncRNA_" + str(num)
                num += 1
            else:
                name = row[0].split(",")[0]
            plt.plot(x,genes, label=name)
        else:
            plt.plot(x,genes,color='lightgrey')
    plt.legend(loc=4, fontsize=12)
    if args.log:
        plt.ylabel("log2(1 + RPKM)", fontsize=12)
    else:
        plt.ylabel("RPKM", fontsize=12)
    plt.xticks(x,labels,rotation=30, fontsize=12)
    plt.tight_layout()
    plt.savefig(args.out_file)
if __name__ == "__main__":
    main()
