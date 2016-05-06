import os
import sys
import csv
import argparse
from gff3 import Gff3Parser
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import colors
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
import matplotlib.patches as patches
import six
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx


__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-t","--table_file",help="output file")
args = parser.parse_args()

def main():
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
                    break            gh.close()
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
    th = open(args.table_file, "r")
    pre_group = None
    for row in csv.reader(th, delimiter='\t'):
        if row[0] != "group":
            if pre_group is None:
                pre_group = row[0]
                gene = row
                quantis = []
            else:
                if row[0] != pre_group:
                    plot(row, quantis)
                    pre_group = row[0]
                    gene = row
                    quantis = []
                if row[3] == "Anti":
                    for poss in row[5].split(";"):
                        strand = pos.split("_")[1]
                        start = pos.split("-")[0]
                        end = pos.split("_")[0].split("-")[-1]
                        for srna in srnas:
                            if (srna["strand"] == strand) and (
                                srna["start"] == start) and (
                                srna["end"] == end):
                                quantis.append({"name": srna["name"], "quanti": srna["quanti"]})
                                break
    plot(row, quantis)
                        
if __name__ == "__main__":
    main()
