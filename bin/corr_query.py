#!/usr/bin/python

import os
import sys
import csv
import shutil
import argparse
from subprocess import call
from scipy.stats.stats import pearsonr, spearmanr
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.spatial import distance

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--gene_quanti_file",help="gene quantification file")
parser.add_argument("-g","--gene_name_file", default=None, help="gene name file")
parser.add_argument("-p","--position", help="gene position and strand (ex. 121434:122441:+)")
parser.add_argument("-pc","--pos_cut", type=float, default=0.8, help="the cutoff of positive correlation, default is 0.8")
parser.add_argument("-nc","--neg_cut", type=float, default=-0.8, help="the cutoff of negative correlation, default is -0.8")
parser.add_argument("-k","--known_srna_only", action="store_true", help="only analyze known sRNA")
parser.add_argument("-kc","--known_srna_color", default="silver", help="the color for known sRNA")
parser.add_argument("-dc","--novel_srna_color", default="silver", help="the color for novel sRNA")
parser.add_argument("-cc","--cds_color", default="silver", help="the color for CDS")
parser.add_argument("-go","--go_file", help="association go file")
parser.add_argument("-gp","--goatools_path", help="path of find_enrichment.py in goatools")
parser.add_argument("-gb","--obo_file", help="path fo go.obo")
parser.add_argument("-po","--population_file", help="Go population file")
parser.add_argument("-ga","--go_association", help="Go association file")
parser.add_argument("-s","--strain", help="Strain name")

args = parser.parse_args()

def get_product(info):
    attrs = "_".join(info.split("_")[3:])
    datas = attrs.split(";")
    for data in datas:
        if data.startswith("product"):
            product = data.split("=")[-1]
            return product
    for data in datas:
        if data.startswith("Name"):
            product = data.split("=")[-1]
            return product

def get_pro_id(info):
    attrs = "_".join(info.split("_")[3:])
    datas = attrs.split(";")
    for data in datas:
        if data.startswith("protein_id"):
            product = data.split("=")[-1]
            return product
    return None

def get_gene_name(genes, info):
    attrs = "_".join(info.split("_")[3:])
    datas = attrs.split(";")
    print(info)
    if (("ID=" + args.strain) not in info):
        for data in datas:
            if data.startswith("locus_tag"):
                locus = data.split("=")[-1]
        for g_locus, name in genes.items():
            if g_locus == locus:
                return name
    else:
        for data in datas:
            if data.startswith("Name"):
                return(data.split("=")[-1])
    return "-"

def gen_input_goatools(filename, infos, exps, cutoff, q_info, q_exp, gos):
    out = open(filename, "w")
    datas = []
    ids = []
    pro_gos = []
    ccs = []
    for info, exp in zip(infos, exps):
        if info != q_info:
            if "positive" in filename:
                if (float(spearmanr(exp, q_exp)[0]) >= cutoff):
                    if get_pro_id(info) is not None:
                        out.write(get_pro_id(info) + "\n")
                    detect = False
                    for pro_id, go_list in gos.items():
                        if pro_id in info:
                            pro_gos.append(go_list)
                            detect = True
                            break
                    if not detect:
                        pro_gos.append("NA")
                    datas.append(exp)
                    ids.append(info)
                    ccs.append("{0:.5f}".format(float(spearmanr(exp, q_exp)[0])))
            elif "negative" in filename:
                if (float(spearmanr(exp, q_exp)[0]) <= cutoff):
                    detect = False
                    if get_pro_id(info) is not None:
                        out.write(get_pro_id(info) + "\n")
                    for pro_id, go_list in gos.items():
                        if pro_id in info:
                            pro_gos.append(go_list)
                            detect = True
                            break
                    if not detect:
                        pro_gos.append("NA")
                    datas.append(exp)
                    ids.append(info)
                    ccs.append("{0:.5f}".format(float(spearmanr(exp, q_exp)[0])))
    out.close()
    out_go = open(filename + "_go", "w")
    call(["python3", args.goatools_path, "--pval=0.05", "--indent",
          "--obo=" + args.obo_file, filename, args.population_file,
          args.go_association], stdout=out_go)
    out_go.close()
    fh = open(filename + "_go", "r")
    start = False
    enrichs = []
    for row in csv.reader(fh, delimiter='\t'):
        if start:
            if row[2] == "e":
                enrichs.append(row[0].replace(".", ""))
        if row[0] == "GO":
            start = True
    fh.close()
    os.remove(filename)
    os.remove(filename + "_go")
    return datas, ids, pro_gos, enrichs, ccs

def plot(filename, q_info, q_exp, pngname, infos, exps, cutoff, genes, gos):
    tags = ["TSB_OD_0.2", "TSB_OD_0.5", "TSB_OD_1", "TSB_t0", "TSB_t1", "TSB_t2", "TSB_ON",
            "pMEM_OD_0.2", "pMEM_OD_0.5", "pMEM_OD_1", "pMEM_t0", "pMEM_t1", "pMEM_t2", "pMEM_ON"]
    datas, ids, pro_gos, enrichs, ccs = gen_input_goatools(
         filename, infos, exps, cutoff, q_info, q_exp, gos)
    f_datas = [q_exp]
    f_ids = [q_info]
    out = open(filename, "w")
    q_items = q_info.split("_")
    product = get_product(q_info)
    gene_name = get_gene_name(genes, q_info)
    out.write("Query gene: " + "_".join(q_items[:3] + [product]) + "\n")
    out.write("start\tend\tstrand\tgene\tgene_name\tGO\tC.C.\n")
    out.write("\t".join(q_items[:3] + [product, gene_name]) + "\n")
    for info, exp, go_list, cc in zip(ids, datas, pro_gos, ccs):
        product = get_product(info)
        gene_name = get_gene_name(genes, info)
        items = info.split("_")
        if go_list != "NA":
            for go in go_list:
                if go in enrichs:
                    f_datas.append(exp)
                    f_ids.append(info)
                    out.write("\t".join(items[:3] + [product, gene_name,
                                        ";".join(go_list), str(cc)]) + "\n")
                    break
        else:
            f_datas.append(exp)
            f_ids.append(info)
            out.write("\t".join(items[:3] + [product, gene_name, "-", str(cc)]) + "\n")
    out.close()
    fig = plt.figure(figsize=(14, 10))
    for data, id_ in zip(f_datas, f_ids):
        if "ID=srna" in id_:
            pass
        else:
            plt.plot(data, color=args.cds_color)
    for data, id_ in zip(f_datas, f_ids):
        if "ID=srna" in id_:
            if "Name=sRNA" in id_:
                plt.plot(data, color=args.novel_srna_color)
            else:
                plt.plot(data, color=args.known_srna_color)
    plt.plot(f_datas[0], color="red")
    plt.ylabel("log2 fold change", fontsize=20)
    x = np.arange(14)
    plt.xticks(x,tags,rotation=30, fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlim([0, 13])
    plt.savefig(pngname)

def main():
    gos = {}
    goh = open(args.go_file, "r")
    for row in csv.reader(goh, delimiter='\t'):
        gos[row[0]] = row[1].split(";")
    start = args.position.split(":")[0]
    end = args.position.split(":")[1]
    strand = args.position.split(":")[2]
    exps = []
    infos = []
    genes = {}
    if args.gene_name_file is not None:
        gh = open(args.gene_name_file, "r")
        for row in csv.reader(gh, delimiter="\t"):
            genes[row[0]] = row[3]
    fh = open(args.gene_quanti_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if (row[0] != "Orientation") and (row[0] == "sense"):
            if (not args.known_srna_only) or (
                   (args.known_srna_only) and ("Name=sRNA" not in row[9])):
                exps.append([float(row[10]), float(row[11]), float(row[12]),
                             float(row[13]), float(row[14]), float(row[15]), float(row[16]),
                             float(row[17]), float(row[18]), float(row[19]),
                             float(row[20]), float(row[21]), float(row[22]), float(row[23])])
            infos.append("_".join([row[4], row[5], row[7], row[9]]))
            if (row[4] == start) and (row[5] == end) and (row[7] == strand):
                q_exp = [float(row[10]), float(row[11]), float(row[12]),
                         float(row[13]), float(row[14]), float(row[15]), float(row[16]),
                         float(row[17]), float(row[18]), float(row[19]),
                         float(row[20]), float(row[21]), float(row[22]), float(row[23])]
                q_info = "_".join([row[4], row[5], row[7], row[9]])
    name = get_product(q_info)
    if "/" in name:
        name = "_".join([name.replace("/", "_or_"), start, end, strand])
    else:
        name = "_".join([name, start, end, strand])
    plot(name + "_positive", q_info, q_exp,
         name + "_positive.png", infos, exps, args.pos_cut, genes, gos)
    plot(name + "_negative", q_info, q_exp,
         name + "_negative.png", infos, exps, args.neg_cut, genes, gos)
#                print(distance.euclidean(exp1, exp2))

if __name__ == "__main__":
    main()
