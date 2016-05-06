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
    pop = open("test_all", "w")
    data = open("test_data", "w")
    fh = open(args.input_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if row[0] != "strain":
            pop.write(row[4] + "\n")
            if len(row[-1]) != 0:
                data.write("\t".join([row[4], row[-1].replace(" ", "")]) + "\n")

if __name__ == "__main__":
    main()
