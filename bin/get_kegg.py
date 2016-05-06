#!/usr/bin/python

import os
import sys
import csv
import argparse

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-k","--kegg_file",help="output file")
args = parser.parse_args()


def main():
    gos = []
    gh = open(args.kegg_file, "r")
    for row in csv.reader(gh, delimiter='\t'):
        gos.append({"kegg": row[0], "locus_tag": row[1]})
    fh = open(args.input_file, "r")
    first = True
    for row in csv.reader(fh, delimiter='\t'):
        if len(row) == 1:
            if first:
                first = False
            else:
                print(int(row[0]) - 1)
                for key, value in stats.items():
                    print(key + "=" + str(value))
            stats = {}
        if len(row) > 1:
            if row[2] == "gene":
                infos = row[8].split(";")
                for info in infos:
                    if "locus_tag" in info:
                        locus = info.split("=")[-1]
                for kegg in gos:
                    if kegg["locus_tag"] == locus:
                        if kegg["kegg"] not in stats.keys():
                            stats[kegg["kegg"]] = 1
                        else:
                            stats[kegg["kegg"]] += 1
#                go_slims = compare_go_slim(go_terms, term_obos, slim_obos)
#                for go_slim in go_slims:
#                    if go_slim not in stats.keys():
#                        stats[go_slim] = 1
#                    else:
#                        stats[go_slim] += 1
#                print("; ".join(go_slims))
    print(row)
    for key, value in stats.items():
        print(key + "=" + str(value))

if __name__ == "__main__":
    main()
