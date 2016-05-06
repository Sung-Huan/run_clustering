#!/usr/bin/python

import os
import sys
import csv
import argparse

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-c","--class_file",help="output file")
args = parser.parse_args()

def main():
    fh = open(args.input_file, "r")
    gos = []
    pre_group = ""
    groups = []
    for row in csv.reader(fh, delimiter='\t'):
        if row[0] != "group":
            if pre_group != row[0]:
                ch = open(args.class_file, "r")
                for pos in row[5].split(";"):
                    for c in csv.reader(ch, delimiter='\t'):
                        if len(c) == 1:
                            num = c[0]
                        else:
                            start = pos.split("-")[0]
                            end = pos.split("-")[1].split("_")[0]
                            strand = pos.split("_")[1]
                            if (start == c[3]) and (
                                end == c[4]) and (strand == c[6]):
                                group = num
                                break
                ch.close()
                ch = open(args.class_file, "r")
                for c in csv.reader(ch, delimiter='\t'):
                    if len(c) == 1:
                        num = c[0]
                        srna_num = 1
                    elif group == num:
                        groups.append(num)
                        datas = c[8].split(";")
                        for data in datas:
                            if c[2] == "CDS":
                                if "Name=" in data:
                                    gene_name = data.split("=")[-1]
                                if "product=" in data:
                                    product = data.split("=")[-1]
                            elif c[2] == "sRNA":
                                if "Name=" in data:
                                    srna_name = data.split("=")[-1]
                                    if "sRNA_00" in srna_name:
                                        gene_name = "novel_" + str(srna_num)
                                        srna_num += 1
                                product = "sRNA"
                        if c[2] != "gene":
                            print("\t".join([row[0], product, gene_name,
                                     c[3] + "-" + c[4] + "_" + c[6]]))
                pre_group = row[0]
    ch.close()
    ch = open(args.class_file, "r")
    detect = False
    for c in csv.reader(ch, delimiter='\t'):
        if len(c) == 1:
            if c[0] not in groups:
                detect = True
            else:
                detect = False
        elif detect:
            datas = c[8].split(";")
            for data in datas:
                if c[2] == "CDS":
                    if "Name=" in data:
                        gene_name = data.split("=")[-1]
                    if "product=" in data:
                        product = data.split("=")[-1]
                elif c[2] == "sRNA":
                    if "Name=" in data:
                        srna_name = data.split("=")[-1]
                        if "sRNA_00" in srna_name:
                            gene_name = "novel_" + str(srna_num)
                            srna_num += 1
                    product = "sRNA"
            if c[2] != "gene":
                print("\t".join(["group_60", product, gene_name,
                         c[3] + "-" + c[4] + "_" + c[6]]))

if __name__ == "__main__":
    main()
