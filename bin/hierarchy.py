import os
import sys
import csv
import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import colors
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import inconsistent
from scipy.cluster.hierarchy import fcluster, fclusterdata
import numpy as np
import math
import random
import six
from ete3 import Tree
from gff3 import Gff3Parser
from numpy import percentile

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-s","--srna_file",help="output file")
parser.add_argument("-g","--gff_file",help="output file")
parser.add_argument("-t","--type_", default=None, help="output file")
parser.add_argument("-d","--max_d", type=float, help="output file")
args = parser.parse_args()

def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)
    ddata = dendrogram(*args, **kwargs)
    colors = []
    color_elements = ["FF", "00", "33", "66", "99", "AA", "DD"]
    color = ["00", "00", "00"]
    num = 0
    if not kwargs.get('no_plot', False):
        plt.title('Dendrogram of hierarchical clustering', fontsize=28)
        plt.tick_params(labelsize=22)
        plt.xlabel('Genes', fontsize=24, color="black")
        plt.ylabel('Distance', fontsize=24, color="black")
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            color_index = random.randint(0, 6)
            if num % 3 == 0:
                color = color[:2] + [color_elements[color_index]]
            elif num % 3 == 1:
                color = [color_elements[color_index]] + color[1:]
            elif num % 3 == 2:
                color = [color[0]] + [color_elements[color_index]] + [color[-1]]
            if ("".join(color) != "FFFFFF") and \
               (color not in colors):
                c = '#' + "".join(color)
                colors.append(color)
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, color=c)
            num += 1
        if max_d:
            plt.axhline(y=max_d, linewidth=3, color='orange')
            
    return ddata

def get_srna(name):
    names = []
    fh = open(args.srna_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if row[2] == name:
            srnas = row[-4].split(";")
            for srna in srnas:
                names.append(srna.split("|")[-2])
            break
    return ",".join(names)

def get_name(row):
    datas = row[9].split(";")
    name = ""
    locus = ""
    if row[3] == "CDS":
        for data in datas:
#            if data[:4] == "Name":
#                name = data[5:]
            if data[:5] == "locus":
                name = data[10:]
#        if name != locus:
#            name = locus + ":" + name
    else:
        for data in datas:
            if data[:4] == "Name":
                name = "|".join([data[5:], row[4], row[5], row[7]])
#            elif data[:8] == "sRNA_hit":
#                if (data[-2:] == "NA"):
#                    srna = name
#                else:
#                    srna = get_srna(name)
#        if name != srna:
#            name = name + ":" + srna
    return name

def main():
    all_rpkms = {"names": [], "rpkms": []}
    srna_rpkms = {"names": [], "rpkms": []}
    gene_rpkms = {"names": [], "rpkms": []}
    gff_f = open(args.gff_file, "r")
    genes = []
    for entry in Gff3Parser().entries(gff_f):
        if entry.feature != "source":
            genes.append(entry)
    libs =  {"TSB_OD_0.2": [], "TSB_OD_0.5": [], "TSB_OD_1": [], "TSB_t0": [], "TSB_t1": [], "TSB_t2": [], "TSB_ON": [],
             "pMEM_OD_0.2": [], "pMEM_OD_0.5": [], "pMEM_OD_1": [], "pMEM_t0": [], "pMEM_t1": [], "pMEM_t2": [], "pMEM_ON": []}
    fh = open(args.input_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if (not row[0].startswith("Orientation")) and (
            row[0] == "sense"):
            gene_name = get_name(row)
            rpkm_row = [float(row[10]), float(row[11]), float(row[12]),
                        float(row[13]), float(row[14]), float(row[15]),
                        float(row[16]), float(row[17]), float(row[18]),
                        float(row[19]), float(row[20]), float(row[21]),
                        float(row[22]), float(row[23])]
            if row[3] == "CDS":
                all_rpkms["names"].append(gene_name)
                gene_rpkms["names"].append(gene_name)
                all_rpkms["rpkms"].append(rpkm_row)
                gene_rpkms["rpkms"].append(rpkm_row)
            elif row[3] == "sRNA":
                all_rpkms["names"].append(gene_name)
                srna_rpkms["names"].append(gene_name)
                all_rpkms["rpkms"].append(rpkm_row)
                srna_rpkms["rpkms"].append(rpkm_row)
    data = np.array(all_rpkms["rpkms"])
    Z = linkage(data, method='ward', metric='euclidean')
    c, coph_dists = cophenet(Z, pdist(data))
    clusters = fcluster(Z, args.max_d, criterion='distance')
    nums = {}
    names = {}
    c_genes = {}
    index = 0
    for c in clusters:
        if c not in nums.keys():
            nums[c] = 1
            names[c] = [all_rpkms["names"][index]]
            c_genes[c] = [all_rpkms["rpkms"][index]]
        else:
            nums[c] += 1
            names[c].append(all_rpkms["names"][index])
            c_genes[c].append(all_rpkms["rpkms"][index])
        index += 1
    print(nums)
#    x = np.arange(14)
#    labels = ["TSB_OD_0.2", "TSB_OD_0.5", "TSB_OD_1", "TSB_t0", "TSB_t1", "TSB_t2", "TSB_ON",
#              "pMEM_OD_0.2", "pMEM_OD_0.5", "pMEM_OD_1", "pMEM_t0", "pMEM_t1", "pMEM_t2", "pMEM_ON"]
#    color_list = list(six.iteritems(colors.cnames))
#    for index, gene_list in c_genes.items():
#        plt.figure(figsize=(12.5, 8))
#        srna_detect = False
#        srna_num = 1
#        color_num = 0
#        for i in range(len(gene_list)):
#            if "sRNA" in names[index][i]:
#                srna_detect = True
#                if ":" in names[index][i]:
#                    srna_name = names[index][i].split(":")[-1]
#                else:
#                    srna_name = "novel_" + str(srna_num)
#                    srna_num += 1
#                if ("grey" not in color_list[color_num][0]) and (
#                    "gray" not in color_list[color_num][0]) and (
#                    "white" not in color_list[color_num][0]) and (
#                    "snow" not in color_list[color_num][0]) and (
#                    color_list[color_num][0] != "w"):
#                    plt.plot(x,gene_list[i], color=color_list[color_num][0], label=srna_name)
#                    color_num += 1
#            else:
#                plt.plot(x,gene_list[i],color='lightgrey')
##        plt.axhline(y=0, linewidth=2, color='red')
#        plt.ylabel("log2 fold change", fontsize=10)
#        plt.xticks(x,labels,rotation=45, fontsize=8)
#        if srna_detect:
#            plt.legend(loc=9, bbox_to_anchor=(1.065, 1), fontsize=8)
#        plt.savefig("test_" + str(index) + ".png")
    for index, gene_names in names.items():
        print(index)
        for name in gene_names:
            for gene in genes:
                if ("locus_tag" in gene.attributes.keys()):
                    if name == gene.attributes["locus_tag"]:
                        print(gene.info)
                elif ("sRNA_hit" in gene.attributes.keys()):
                    infos = name.split("|")
                    if (infos[0] == gene.attributes["Name"]) and (
                        infos[1] == str(gene.start)) and (
                        infos[2] == str(gene.end)) and (
                        infos[3] == gene.strand):
                        print(gene.info)
    plt.style.use('ggplot')
    plt.figure(figsize=(25, 10))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('Genes')
    plt.ylabel('distance')
    fancy_dendrogram(
        Z,
#        truncate_mode='lastp',
#        p=12,
        leaf_rotation=90.,
        leaf_font_size=12.,
#        show_contracted=True,
#        annotate_above=10,
        no_labels=True,
        show_leaf_counts=False,
        max_d=args.max_d,  # plot a horizontal cut-off line
    )
    plt.savefig("hierarchical_tree.png")

if __name__ == "__main__":
    main()

