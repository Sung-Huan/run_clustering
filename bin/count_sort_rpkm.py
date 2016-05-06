#!/usr/bin/python

import os
import sys
import csv
import argparse

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-o","--output_file",help="output file")
args = parser.parse_args()

def main():
    genes = {}
    fh = open(args.input_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if row[0] != "TSB_OD_0.2":
            genes[row[0]] = [-1, -1, -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1]
    fh.close()
    rank = 0
    fh = open(args.input_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if row[0] != "TSB_OD_0.2":
            rank += 1
            for index in range(len(row)):
                genes[row[index]][index] = rank
    for gene, ranks in genes.items():
        total = 0
        for rank in ranks:
            total = total + rank
        ranks = map(str, ranks)
        print("\t".join([gene] + ranks + [str(total)]))

if __name__ == "__main__":
    main()
