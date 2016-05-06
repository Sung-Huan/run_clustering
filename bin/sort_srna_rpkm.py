#!/usr/bin/python

import os
import sys
import csv
import argparse
import numpy as np
import math
import matplotlib.pyplot as plt

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-s","--srna_file",help="output file")
args = parser.parse_args()

def get_srna(data):
    names = []
    fh = open(args.srna_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if (data[1] == row[1]) and (data[4] == row[3]) and (
            data[5] == row[4]) and (data[7] == row[5]):
            srnas = row[-4].split(";")
            for srna in srnas:
                names.append(srna.split("|")[-2])
            break
    return ",".join(names)

def get_name(row):
    datas = row[9].split(";")
    name = ""
    for data in datas:
        if data[:4] == "Name":
            name = data[5:]
    return name

def main():
    fh = open(args.input_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if (not row[0].startswith("Orientation")) and (
            row[0] == "sense") and (row[3] != "CDS"):
            gene_name = get_name(row)
            row = [row[10], row[12], row[14],
                  row[18], row[20], row[22],
                  row[16], row[24], row[26], row[28],
                  row[32], row[34], row[36], row[30]]
            print("\t".join([gene_name] + row))

if __name__ == "__main__":
    main()
