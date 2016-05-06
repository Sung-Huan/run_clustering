#!/usr/bin/python

import os
import sys
import csv
import argparse
from scipy.stats.stats import spearmanr
import numpy

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-t","--table_file",help="input file")
parser.add_argument("-g","--gene_file",help="output file")
parser.add_argument("-s","--srna_file",help="output file")
args = parser.parse_args()

def get_guanti(row, genes):
    quantis = []
    for pos in row[5].split(";"):
        strand = pos.split("_")[1]
        start = pos.split("-")[0]
        end = pos.split("_")[0].split("-")[-1]
        for gene in genes:
            if (gene["strand"] == strand) and (
                gene["start"] == start) and (
                gene["end"] == end):
                quantis.append(gene["quanti"])
                break
    return quantis

def main():
    genes = []
    srnas = []
    gh = open(args.gene_file, "r")
    for row in csv.reader(gh, delimiter='\t'):
        if row[3] == "sRNA":
            datas = row[9].split(";")
            for data in datas:
                if "Name=" in data:
                    name = data.split("=")[-1]
            if "sRNA_00" in name:
                 name = "novel_sRNA"
            sh = open(args.srna_file, "r")
            for srna in csv.reader(sh, delimiter='\t'):
                if (srna[3] == row[4]) and (
                    srna[4] == row[5]) and (srna[5] == row[7]):
                    srnas.append({"feature": row[3], "start": row[4], "cover": srna[10],
                          "end": row[5], "strand": row[7], "name": name,
                          "quanti": [float(row[10]), float(row[11]), float(row[12]),
                          float(row[13]), float(row[14]), float(row[15]), float(row[16]),
                          float(row[17]), float(row[18]), float(row[19]),
                          float(row[20]), float(row[21]), float(row[22]), float(row[23])]})
        genes.append({"feature": row[3], "start": row[4],
                      "end": row[5], "strand": row[7],
                      "quanti": [float(row[10]), float(row[11]), float(row[12]),
                      float(row[13]), float(row[14]), float(row[15]), float(row[16]),
                      float(row[17]), float(row[18]), float(row[19]),
                      float(row[20]), float(row[21]), float(row[22]), float(row[23])]})

    fh = open(args.table_file, "r")
    pre_group = None
    for row in csv.reader(fh, delimiter='\t'):
        print("\t".join(row))
        if row[0] != "group":
            if pre_group is None:
                pre_group = row[0]
                quantis = get_guanti(row, genes)
            else:
                if pre_group == row[0]:
                    quantis = quantis + get_guanti(row, genes)
                else:
                    antis = {}
                    for srna in srnas:
                        corrs = []
                        for quanti in quantis:
                            corrs.append(spearmanr(srna["quanti"], quanti)[0])
                        if numpy.median(corrs) <= -0.8:
                            if srna["name"] not in antis.keys():
                                antis[srna["name"]] = {"num": 1, "info": [srna]}
                            else:
                                antis[srna["name"]]["num"] += 1
                                antis[srna["name"]]["info"].append(srna)
                    for name, anti in antis.items():
                        covers = []
                        poss = []
                        for data in anti["info"]:
                            covers.append(data["start"] + "-" + data["end"] + "_" + data["strand"] + "_" + data["cover"])
                            poss.append(data["start"] + "-" + data["end"] + "_" + data["strand"])
                        print("\t".join([pre_group, name, str(anti["num"]), "Anti", ";".join(covers),
                                         ";".join(poss), "NA"]))
                    pre_group = row[0]
                    quantis = get_guanti(row, genes)
                

if __name__ == "__main__":
    main()
