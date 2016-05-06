#!/usr/bin/python

import os
import sys
import csv
import argparse
from scipy.stats.stats import pearsonr, spearmanr
import numpy
from scipy.spatial import distance

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-o","--output_file",help="output file")
args = parser.parse_args()

def main():
    exps = []
    fh = open(args.input_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if (row[0] != "Orientation") and (row[0] == "sense"):
            exps.append([float(row[10]), float(row[11]), float(row[12]),
                         float(row[13]), float(row[14]), float(row[15]), float(row[16]),
                         float(row[17]), float(row[18]), float(row[19]),
                         float(row[20]), float(row[21]), float(row[22]), float(row[23])])
    corrs = []
    for exp1 in exps:
        for exp2 in exps:
            if exp1 != exp2:
                print(spearmanr(exp1, exp2)[0])
#                print(distance.euclidean(exp1, exp2))

if __name__ == "__main__":
    main()
