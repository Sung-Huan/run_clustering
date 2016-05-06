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
parser.add_argument("-g","--gff_file",help="output file")
parser.add_argument("-n","--gene_file",help="output file")
parser.add_argument("-c","--class_file",help="output file")
parser.add_argument("-l","--go_level", type=int, help="output file")
parser.add_argument("-p","--pan_file", help="output file")
args = parser.parse_args()

def get_class(info):
    start = info.split("-")[0]
    end = info.split("-")[1].split("_")[0]
    strand = info.split("_")[1]
    groups = []
    fh = open(args.class_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if len(row) == 1:
            group = row[0]
        if len(row) != 1:
            if (start == row[3]) and (end == row[4]) and (strand == row[6]):
                break
    fh.close()
    fh = open(args.class_file, "r")
    detect = False
    for row in csv.reader(fh, delimiter='\t'):
        if (len(row) == 1) and (row[0] == group):
            detect = True
        elif (len(row) == 1) and (row[0] != group):
            detect = False
        if (len(row) != 1) and (detect):
            groups.append({"start": row[3], "end": row[4], "strand": row[6], "feature": row[2]})
    fh.close()
    return groups

def plot(datas, x, type_, color_list, color_num, srna_num, product):
    if type_ == "CDS":
        groups = get_class(datas["info"][0])
        for info in groups:
            fh = open(args.gene_file, "r")
            for row in csv.reader(fh, delimiter='\t'):
                if (not row[0].startswith("Orientation")) and (
                    row[0] == "sense") and (row[3] == "CDS"):
                    if (info["start"] == row[4]) and (info["end"] == row[5]) and (
                        info["strand"] == row[7]):
                        rpkm_row = [float(row[10]), float(row[11]), float(row[12]),
                                    float(row[13]), float(row[14]), float(row[15]),
                                    float(row[16]), float(row[17]), float(row[18]),
                                    float(row[19]), float(row[20]), float(row[21]),
                                    float(row[22]), float(row[23])]
                        if (type_ == "CDS") and (info["feature"] == "CDS"):
                            plt.plot(x,rpkm_row,color='silver')
                            break
            fh.close()
    elif (type_ == "sRNA"):
        for info in datas["info"]:
            start = info.split("-")[0]
            end = info.split("-")[1].split("_")[0]
            strand = info.split("_")[1]
            fh = open(args.gene_file, "r")
            for row in csv.reader(fh, delimiter='\t'):
                if (not row[0].startswith("Orientation")) and (
                    row[0] == "sense") and (row[3] == "sRNA"):
                    if (start == row[4]) and (end == row[5]) and (
                        strand == row[7]):
                        rpkm_row = [float(row[10]), float(row[11]), float(row[12]),
                                    float(row[13]), float(row[14]), float(row[15]),
                                    float(row[16]), float(row[17]), float(row[18]),
                                    float(row[19]), float(row[20]), float(row[21]),
                                    float(row[22]), float(row[23])]
                        while(1):
                            if ("grey" not in color_list[color_num][0]) and (
                                "gray" not in color_list[color_num][0]) and (
                                "white" not in color_list[color_num][0]) and (
                                "snow" not in color_list[color_num][0]) and (
                                color_list[color_num][0] != "w"):
                                break
                            else:
                                color_num += 1
                        if "novel" in product:
                            plt.plot(x,rpkm_row, color=color_list[color_num][0], label="ncRNA_" + str(srna_num))
                            srna_num += 1
                        else:
                            plt.plot(x,rpkm_row, color=color_list[color_num][0], label=product)
                        color_num += 1
            fh.close()
    return color_num, srna_num

def get_gene_name(infos, genes, pans):
    names = []
    for info in infos.split(", "):
        start = int(info.split("-")[0])
        end = int(info.split("-")[1].split("_")[0])
        strand = info.split("_")[1]
        name = "NA"
        for gene in genes:
            if gene.strand == strand:
                if (gene.start == start) and (
                    gene.end == end):
                    name = gene.attributes["Name"]
                    break
                elif (gene.start <= start) and (
                    gene.end >= end):
                    name = gene.attributes["Name"]
        if "SAOUHSC" in name:
            for pan in pans:
                if name == pan["locus"]:
                    if pan["name"] != "-":
                        name = pan["name"]
                    break
        names.append(name)
    return names

def get_pan(name, pans):
    for pan in pans:
        if name == pan["locus"]:
            if pan["name"] != "-":
                name = pan["name"]
            break
    return name

def main():
    pans = []
    ph = open(args.pan_file, "r")
    for row in csv.reader(ph, delimiter='\t'):
        pans.append({"locus": row[0], "name": row[3]})
    gff_f = open(args.gff_file, "r")
    genes = []
    names = []
    for entry in Gff3Parser().entries(gff_f):
        if (entry.feature == "CDS"):
            genes.append(entry)
        if (entry.feature == "gene"):
            names.append(entry)
    pros = []
    fh = open(args.input_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if row[3] != "sRNA":
            if ("." in row[1]):
                level = 0
                for dot in row[1]:
                    if dot == ".":
                        level += 1
                gene_names = get_gene_name(row[-1], names, pans)
                pros.append({"group": row[0], "go": row[1].replace(".", "") + "(" + \
                             row[3] + ", " + row[4] + ")", "go_name": row[3], "term": row[1].replace(".", ""),
                             "product": row[-2].split(", "), "num": 0, "pos": row[-1].split(", "), "lv": level,
                             "gene_name": gene_names})
        else:
            pros.append({"group": row[0], "go": row[1].replace(".", ""), "num": row[6],
                         "go_name": row[3], "product": row[-1].split(", "), "pos": row[-1].split(", ")})
    pre_group = None
    finals = {}
    for pro in pros:
        if pro["group"] != pre_group:
            finals[pro["group"]] = {}
            pre_group = pro["group"]
        if pro["go_name"] != "sRNA": 
            for product in pro["product"]:
                if product not in finals[pro["group"]].keys():
                    finals[pro["group"]][product] = {"go": [], "num": 0, "term": [], "lv": []}
                    num = 0
                    infos = []
                    gene_names = []
                    poss = []
                    for pos in pro["pos"]:
                        start = pos.split("-")[0]
                        end = pos.split("-")[1].split("_")[0]
                        strand = pos.split("_")[1]
                        if "-".join([start, end, strand]) not in poss:
                            poss.append("-".join([start, end, strand]))
                            for gene in genes:
                                if (gene.attributes["product"] == product) and (
                                    gene.start == int(start)) and (
                                    gene.end == int(end)) and (gene.strand == strand):
                                    num += 1
                                    parent = gene.attributes["Parent"]
                                    infos.append(str(gene.start) + "-" + str(gene.end) + "_" + gene.strand)
                                    for name in names:
                                        if parent == name.attributes["ID"]:
                                            gene_names.append(get_pan(name.attributes["Name"], pans))
                                            break
                    finals[pro["group"]][product]["num"] = str(num)
                    finals[pro["group"]][product]["info"] = infos
                    finals[pro["group"]][product]["gene_names"] = gene_names
                finals[pro["group"]][product]["go"].append(pro["go"])
                finals[pro["group"]][product]["term"].append(pro["term"])
                finals[pro["group"]][product]["lv"].append(pro["lv"])
        else:
            finals[pro["group"]][pro["go"]] = {"go": ["sRNA"], "num": pro["num"], "info": pro["product"]}
    printed = []
    group_num = 1
    labels = ["TSB_OD_0.2", "TSB_OD_0.5", "TSB_OD_1", "TSB_t0", "TSB_t1", "TSB_t2", "TSB_ON",
              "pMEM_OD_0.2", "pMEM_OD_0.5", "pMEM_OD_1", "pMEM_t0", "pMEM_t1", "pMEM_t2", "pMEM_ON"]
    color_list = list(six.iteritems(colors.cnames))
    x = np.arange(14)
    while True:
        first = True
        for group, products in finals.items():
            if group not in printed:
                if first:
                    less = len(finals[group])
                    gos = {"group": group, "product": products}
                    first = False
                else:
                    if less > len(finals[group]):
                        less = len(finals[group])
                        gos = {"group": group, "product": products}
        color_num = 0
        srna_num = 1
        srna_detect = False
        plot_cds = True
        length = 0
        if not first:
            fig = plt.figure(figsize=(14, 8))
#            gs = gridspec.GridSpec(2, 1,height_ratios=[1,2])
#            ax1 = plt.subplot(gs[0])
            tables = []
            terms = []
            num_rows = 0
            for product, datas in gos["product"].items():
                if datas["go"] != ["sRNA"]:
                    if int(datas["num"]) != 0:
                        top = "NA"
                        for goi in range(len(datas["lv"])):
                            level = datas["lv"][goi]
                            if top == "NA":
                                top = level
                                go_cut = datas["go"][goi]
                            else:
                                if top >= args.go_level:
                                    if (level < top) and (level >= args.go_level):
                                        top = level
                                        go_cut = datas["go"][goi]
                                    elif level == top:
                                        go_cut = go_cut + "; " + datas["go"][goi]
                                else:
                                    if level > top:
                                        top = level
                                        go_cut = datas["go"][goi]
                                    elif level == top:
                                        go_cut = go_cut + "; " + datas["go"][goi]
                        print("\t".join(["group_" + str(group_num), product, datas["num"], go_cut,
                                         ";".join(datas["go"]), ";".join(datas["info"]), ";".join(datas["gene_names"])]))
                        num_rows += 1
#                        tables.append([product, datas["num"], ""])
#                        terms.append([product, datas["num"], datas["term"]])
                    if plot_cds:
                        color_num, srna_num = plot(datas, x, "CDS", color_list, color_num, srna_num, product)
                        plot_cds = False
            for product, datas in gos["product"].items():
                if datas["go"] == ["sRNA"]:
                    srna_detect = True
                    color_num, srna_num = plot(datas, x, "sRNA", color_list, color_num, srna_num, product)
                    print("\t".join(["group_" + str(group_num), product, datas["num"],
                                     ";".join(datas["go"]), ";".join(datas["info"]), "NA"]))
            printed.append(gos["group"])
            plt.ylabel("log2 fold change", fontsize=16)
            plt.xticks(x,labels,rotation=30, fontsize=14)
            plt.yticks(fontsize=14)
            plt.xlim([0, 13])
            if srna_detect:
#                plt.legend(loc=9, bbox_to_anchor=(1.065, 1), fontsize=14)
                plt.legend(loc=0, fontsize=12, ncol=2)
#            columns = ["name", "number", "GO term"]
#            ax2 = plt.subplot(gs[1])
#            ax2.axis('tight')
#            ax2.axis('off')
#            the_table = ax2.table(cellText=tables,
#                        colLabels=columns,
#                        loc='top', cellLoc="left", bbox=[0, 1 - 0.0175*num_rows, 1, 0.0175*num_rows])
#            cellDict=the_table.get_celld()
#            for index in range(num_rows + 1):
#                cellDict[(index,1)].set_width(0.05)
#                cellDict[(index,2)].set_width(0.6)
#                cellDict[(index,0)].set_height(0.0175)
#                cellDict[(index,1)].set_height(0.0175)
#                cellDict[(index,2)].set_height(0.0175)
#            print(num_rows)
#            print((-0.0105, 0.0190972*num_rows), 0.5- 0.0175*num_rows)
#            table_gos = []
#            for table1 in terms:
#                print(table1)
#                num_gos = {}
#                for table2 in terms:
#                    if table1 != table2:
#                        for term1 in table1[2]:
#                            if term1 not in num_gos.keys():
#                                num_gos[term1] = 1
#                            for term2 in table2[2]:
#                                if term1 == term2:
#                                    num_gos[term1] += 1
#                table_gos.append(sorted(num_gos.values(), reverse=True))
#            cm = plt.get_cmap('RdBu_r') 
#            cNorm  = colors.Normalize(vmin=1, vmax=num_rows)
#            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
#            for y in range(len(table_gos)):
#                x = 0
#                for go in table_gos[y]:
#                    verts = [
#                             (-0.0105 + x*0.0005, 0.05 - 0.00175*(float(num_rows)/float(num_rows+1)) - y*(0.00175*float(num_rows)/float(num_rows+1))),
#                             (-0.011 + x*0.0005, 0.05 - 0.00175*(float(num_rows)/float(num_rows+1)) - y*(0.00175*float(num_rows)/float(num_rows+1))),
#                             (-0.011 + x*0.0005, 0.05 - 2*0.00175*(float(num_rows)/float(num_rows+1)) - y*(0.00175*float(num_rows)/float(num_rows+1))),
#                             (-0.0105 + x*0.0005, 0.05 - 2*0.00175*(float(num_rows)/float(num_rows+1)) - y*(0.00175*float(num_rows)/float(num_rows+1))),
#                             (-0.0105 + x*0.0005, 0.05 - 0.00175*(float(num_rows)/float(num_rows+1)) - y*(0.00175*float(num_rows)/float(num_rows+1)))
#                             ]
#                    codes = [Path.MOVETO,
#                             Path.LINETO,
#                             Path.LINETO,
#                             Path.LINETO,
#                             Path.CLOSEPOLY,
#                             ]
#                    path = Path(verts, codes)
#                    patch = patches.PathPatch(path, facecolor=scalarMap.to_rgba(go), lw=1)
#                    ax2.add_patch(patch)
#                    x += 1
#            norm = colors.Normalize(vmin=0, vmax=42)
#            cax = fig.add_axes([0.1, 0.05, 0.8, 0.02])
#            cb1 = matplotlib.colorbar.ColorbarBase(cax, cmap="RdBu_r",
#                                norm=norm, spacing='proportional',
#                                orientation='horizontal')
#            fig.tight_layout()
            plt.tight_layout()
            plt.savefig("group_" + str(group_num) + ".png")
            group_num += 1
        if first:
            break

if __name__ == "__main__":
     main()
