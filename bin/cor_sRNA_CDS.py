#!/usr/bin/python

import os
import sys
import csv
import argparse
from scipy.stats.stats import pearsonr, spearmanr
import numpy
from scipy.spatial import distance
from gff3 import Gff3Parser

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-g","--gff_file",help="output file")
parser.add_argument("-o","--go_file",help="output file")
parser.add_argument("-p","--pan_file",help="output file")
args = parser.parse_args()

def main():
    pans = []
    ph = open(args.pan_file, "r")
    for row in csv.reader(ph, delimiter='\t'):
        pans.append({"locus": row[0], "name": row[3]})
    gos = {}
    oh = open(args.go_file, "r")
    for row in csv.reader(oh, delimiter='\t'):
        gos[row[0]] = row[1]
    exps = []
    srnas = []
    fh = open(args.input_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if row[3] == "CDS":
            datas = row[9].split(";")
            for data in datas:
                if "Parent=" in data:
                    parent = data.split("=")[-1]
                elif "product=" in data:
                    product = data.split("=")[-1]
                elif "protein_id=" in data:
                    protein = data.split("=")[-1]
            gh = open(args.gff_file, "r")
            for entry in Gff3Parser().entries(gh):
                if (entry.feature == "gene") and (
                    entry.attributes["ID"] == parent):
                    name = entry.attributes["Name"]
                    for pan in pans:
                        if entry.attributes["Name"] == pan["locus"]:
                            if pan["name"] != "-":
                                name = pan["name"]
                            break
                    break
            gh.close()
            term = "NA"
            for id_, go in gos.items():
                if id_ == protein:
                    term = go
        elif row[3] == "sRNA":
            datas = row[9].split(";")
            for data in datas:
                if "Name=" in data:
                    name = data.split("=")[-1]
            product = "sRNA"
            protein = "NA"
            term = "NA"
        if (row[0] != "Orientation") and (row[0] == "sense"):
            exps.append({"info": "\t".join(row[1:10]), "pos": row[4] + "-" + row[5] + "_" + row[7],
                         "name": name, "product": product, "protein": protein, "go": term,
                         "quanti": [float(row[10]), float(row[11]), float(row[12]),
                          float(row[13]), float(row[14]), float(row[15]), float(row[16]),
                          float(row[17]), float(row[18]), float(row[19]),
                          float(row[20]), float(row[21]), float(row[22]), float(row[23])]})
        if (row[0] != "Orientation") and (row[0] == "sense") and (
            row[3] == "sRNA"):
            srnas.append({"info": "\t".join(row[1:10]), "pos": row[4] + "-" + row[5] + "_" + row[7],
                          "name": name, "product": product, "protein": protein, "go": term,
                          "quanti": [float(row[10]), float(row[11]), float(row[12]),
                           float(row[13]), float(row[14]), float(row[15]), float(row[16]),
                           float(row[17]), float(row[18]), float(row[19]),
                           float(row[20]), float(row[21]), float(row[22]), float(row[23])]})
    num = 1
    for srna in srnas:
        print(num)
        print("\t".join([srna["name"], srna["pos"]]))
        for exp in exps:
            if srna["info"] != exp["info"]:
                if (spearmanr(srna["quanti"], exp["quanti"])[0] >= 0.9):
                    print("\t".join(["plus_" + str(spearmanr(srna["quanti"], exp["quanti"])[0]), exp["name"], exp["product"], exp["protein"], exp["go"], exp["pos"]]))
                elif (spearmanr(srna["quanti"], exp["quanti"])[0] <= -0.9):
                    print("\t".join(["minus_" + str(spearmanr(srna["quanti"], exp["quanti"])[0]), exp["name"], exp["product"], exp["protein"], exp["go"], exp["pos"]]))
        num += 1
#                print(distance.euclidean(exp1, exp2))

if __name__ == "__main__":
    main()
