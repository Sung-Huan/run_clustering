#!/usr/bin/python

import os
import sys
import csv
import argparse
from gff3 import Gff3Parser

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-g","--go_file",help="output file")
parser.add_argument("-o","--go_obo",help="output file")
parser.add_argument("-s","--goslim_obo",help="output file")
parser.add_argument("-f","--gff_file",help="output file")
parser.add_argument("-r","--srna_file",help="output file")
parser.add_argument("-t","--go_type",help="output file")
parser.add_argument("-c","--go_goslim_cds_table",help="output file")
parser.add_argument("-pg","--per_go",type=float, help="output file")
parser.add_argument("-pc","--per_cds",type=float, help="output file")
args = parser.parse_args()

def import_obo(filename):
    obos = []
    start = False
    with open(filename, "r") as o_h:
        for line in o_h:
            line = line.strip()
            if line == "[Term]":
                obo = {}
                start = True
            elif start:
                if len(line) == 0:
                    obos.append(obo.copy())
                    start = False
                else:
                    datas = line.split(": ")
                    if datas[0] == "is_a":
                        if "is_a" not in obo.keys():
                            obo["is_a"] = []
                        obo["is_a"].append(datas[1].strip())
                    else:
                        obo[datas[0]] = datas[1].strip()
    return obos

def compare_go_slim(go, term_obos, slim_obos, stats, genes, row, compares):
    target_terms = [go]
    detect = False
    for target_term in target_terms:
        for term_obo in term_obos:
            if target_term == term_obo["id"]:
                if "is_a" in term_obo.keys():
                    for is_a in term_obo["is_a"]:
                        go_a = is_a.split(" ! ")
                        if (go_a[1] != "biological_process") and (
                            go_a[1] != "cellular_component") and (
                            go_a[1] != "molecular_function"):
                            target_terms.append(go_a[0])
                elif ("is_obsolete" in term_obo.keys()):
                    if term_obo["is_obsolete"] == "true":
                        break
                for slim_obo in slim_obos:
                    for target_term in target_terms:
                        if target_term == slim_obo["id"]:
                            detect = True
                            go_slim = slim_obo["id"]
                break
        if detect:
            detect = False
            detect_slim = False
            if go_slim not in stats.keys():
                stats[go_slim] = {"go_num": 1, "products": [], "name": ""}
                compares.append(go_slim)
                detect_slim = True
            else:
                if go_slim not in compares:
                    stats[go_slim]["go_num"] += 1
                    compares.append(go_slim)
                    detect_slim = True
            if detect_slim:
                for obo in slim_obos:
                    if obo["id"] == go_slim:
                        stats[go_slim]["name"] = obo["name"]
                for gene in genes:
                    if gene.feature == "CDS":
                        if (int(row[3]) == gene.start) and (
                            int(row[4]) == gene.end) and (
                            row[6] == gene.strand):
                            stats[go_slim]["products"].append(gene.attributes["product"])
            break

def main():
    go_nums = {}
    ch = open(args.go_goslim_cds_table, "r")
    rh = open(args.srna_file, "r")
    srnas = []
    for row in csv.reader(rh, delimiter='\t'):
        infos = row[-4].split(";")
        for info in infos:
            if info != "NA":
                srnas.append({"name": row[2], "homo": info.split("|")[-2],
                              "start": row[3], "end": row[4], "strand": row[5]})
            else:
                srnas.append({"name": row[2], "homo": "noval",
                              "start": row[3], "end": row[4], "strand": row[5]})
            break
    for row in csv.reader(ch, delimiter='\t'):
        go_nums[row[0]] = len(row[1].split(";"))
    gff_f = open(args.gff_file, "r")
    genes = []
    for entry in Gff3Parser().entries(gff_f):
        genes.append(entry)
    term_obos = import_obo(args.go_obo)
    slim_obos = import_obo(args.goslim_obo)
    gos = []
    gh = open(args.go_file, "r")
    for row in csv.reader(gh, delimiter='\t'):
        if row[0] != "strain":
            gos.append({"strain": row[0], "strand": row[1],
                        "start": row[2], "end": row[3], "go": row[5]})
    fh = open(args.input_file, "r")
    first = True
    group = 1
    srna_names = {}
    for row in csv.reader(fh, delimiter='\t'):
        if len(row) == 1:
            index = row[0]
            if first:
                num = 0
                srna_nums = 0
                first = False
            else:
                proteins = {}
                for go_term, datas in stats.items():
                    if datas["go_num"] >= (go_nums[go_term]*args.per_go) or (
                       datas["go_num"] >= (num*args.per_cds)):
                        print("\t".join(["group_" + str(group), 
                                     str(num), go_term, datas["name"], str(datas["go_num"]),
                                     ";".join([str(float(datas["go_num"])/float(go_nums[go_term])),
                                               str(float(datas["go_num"])/float(num))]),
                                     ", ".join(datas["products"])]))
                for srna, homos in srna_names.items():
                    print("\t".join(["group_" + str(group),
                                     str(srna_nums), "sRNA", srna, str(homos["num"]),
                                     str(float(homos["num"])/float(srna_nums)),
                                     ", ".join(homos["name"])]))
                group += 1
                num = 0
                srna_nums = 0
            stats = {}
            stats_goslim = {}
            srna_names = {}
        if len(row) > 1:
            if row[2] == "CDS":
                go_terms = []
                for go in gos:
                    if (row[0] == go["strain"]) and (
                        row[3] == go["start"]) and (
                        row[4] == go["end"]) and (
                        row[6] == go["strand"]):
                        go_terms = go["go"].split("; ")
                if len(go_terms) != 0:
                    num += 1
                    compares = []
                    for go_term in go_terms:
                        if "GO" in go_term:
                            if args.go_type == "goslim":
                                compare_go_slim(go_term, term_obos, slim_obos,
                                                stats, genes, row, compares)
                            if args.go_type == "go":
                                if go_term not in stats.keys():
                                    stats[go_term] = {"go_num": 1, "products": [], "name": ""}
                                else:
                                    stats[go_term]["go_num"] += 1
                                for obo in term_obos:
                                    if obo["id"] == go_term:
                                        stats[go_term]["name"] = obo["name"]
                                for gene in genes:
                                    if gene.feature == "CDS":
                                        if (int(row[3]) == gene.start) and (
                                            int(row[4]) == gene.end) and (
                                            row[6] == gene.strand):
                                            stats[go_term]["products"].append(gene.attributes["product"])
            elif row[2] == "sRNA":
                srna_nums += 1
                for srna in srnas:
                    if (row[3] == srna["start"]) and (
                        row[4] == srna["end"]) and (
                        row[6] == srna["strand"]):
                        if (srna["homo"] not in srna_names.keys()):
                            srna_names[srna["homo"]] = {"name": [row[3] + "-" + row[4] + "_" + row[6]], "num": 1}
                        else:
                            srna_names[srna["homo"]]["num"] += 1
                            srna_names[srna["homo"]]["name"].append(row[3] + "-" + row[4] + "_" + row[6])
                
    proteins = {}
    for go_term, datas in stats.items():
        if datas["go_num"] >= (go_nums[go_term]*args.per_go) or (
           datas["go_num"] >= (num*args.per_cds)):
            print("\t".join(["group_" + str(group),
                         str(num), go_term, datas["name"], str(datas["go_num"]),
                         ";".join([str(float(datas["go_num"])/float(go_nums[go_term])),
                                   str(float(datas["go_num"])/float(num))]),
                         ", ".join(datas["products"])]))

    for srna, homos in srna_names.items():
        print("\t".join(["group_" + str(group),
                         str(srna_nums), "sRNA", srna, str(homos["num"]),
                         str(float(homos["num"])/float(srna_nums)),
                         ", ".join(homos["name"])]))

if __name__ == "__main__":
    main()
