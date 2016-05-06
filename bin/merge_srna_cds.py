#!/usr/bin/python

import os
import sys
import csv
import argparse

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-c","--class_file",help="input file")
parser.add_argument("-g","--group_file",help="output file")
parser.add_argument("-o","--go_file",help="output file")
parser.add_argument("-s","--srna_file",help="output file")
args = parser.parse_args()

def main():
    print("\t".join(["group_id", "GO_term(dots for level)", "GO_type", "GO_class_name",
                     "p-value", "total_number(GO_mapped)", "counts_this_GO", "protein_products"]))
    gos = []
    oh = open(args.go_file, "r")
    for row in csv.reader(oh, delimiter='\t'):
        gos.append({"protein": row[0], "product": row[1], "go": row[2]})
    gh = open(args.group_file, "r")
    line = 0
    start = False
    group = "group_" + args.group_file.split("_")[-1]
    for row in csv.reader(gh, delimiter='\t'):
        line += 1
        if line == 4:
            total = row[0].strip().split(" ")[0]
        if row[0] == "GO":
            start = True
        else:
            if start:
                term = row[0].split(".")[-1]
                finals = []
                poss = []
                for go in gos:
                    if term in go["go"]:
                        detect = False
                        ch = open(args.class_file, "r")
                        for cl in csv.reader(ch, delimiter='\t'):
                            if (len(cl) == 1) and (cl[0] == group.split("_")[-1]):
                                detect = True
                            elif (len(cl) == 1) and (cl[0] != group.split("_")[-1]):
                                detect = False
                            if (len(cl) != 1) and (detect):
                                datas = cl[8].split(";")
                                for data in datas:
                                    if ("product" in data) and (data.split("=")[-1] == go["product"]):
                                        poss.append(cl[3] + "-" + cl[4] + "_" + cl[6])
                        finals.append(go["product"])
                        ch.close()
                if (row[2] == "e") and (len(poss) != 0):
                    print("\t".join([group, row[0], row[1], row[3], row[6],
                                     total, str(len(finals)), ", ".join(finals), ", ".join(poss)]))
    if start:
        ch = open(args.class_file, "r")
        detect = False
        stats = {}
        srna_total = 0
        for row in csv.reader(ch, delimiter='\t'):
            if len(row) == 1:
                if args.group_file.split("_")[-1] == row[0]:
                    detect = True
                else:
                    detect = False
            else:
                if detect:
                    if row[2] == "sRNA":
                        srna_total += 1
                        sh = open(args.srna_file, "r")
                        for r in csv.reader(sh, delimiter='\t'):
                            if (r[3] == row[3]) and (
                                r[4] == row[4]) and (
                                r[5] == row[6]):
                                srnas = r[-4].split(";")
                                if srnas[0] == "NA":
                                    name = "novel_sRNA"
                                else:
                                    name = r[2]
#                                    for srna in srnas:
#                                        name = srna.split("|")[-2]
                                break
                        if name not in stats.keys():
                            stats[name] = {"member": [r[3] + "-" + r[4] + "_" + r[5]], "num": 1}
                        else:
                            stats[name]["member"].append(r[3] + "-" + r[4] + "_" + r[5])
                            stats[name]["num"] += 1
        for name, datas in stats.items():
            print("\t".join([group, name, "NA", "sRNA", "NA",
                             str(srna_total), str(datas["num"]), ", ".join(datas["member"]), ", ".join(datas["member"])]))    

if __name__ == "__main__":
    main()
