import os
import sys
import csv
import argparse
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-o","--output_file",help="output file")
args = parser.parse_args()

def main():
    tables = []
    fh = open(args.input_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if (row[2] != "sRNA") and (row[0] == "group_28"):
            datas = row[3].split(";")
            gos = []
            for data in datas:
                gos.append(data.split("(")[0])
            tables.append([row[1], row[2]])
    plt.figure(figsize=(25, 10))
    columns = ["name", "number"]
    plt.table(cellText=tables,
              colLabels=columns,
              loc='bottom')
    plt.savefig("test.png")    

if __name__ == "__main__":
    main()
