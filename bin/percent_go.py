import os
import sys
import csv
import argparse

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-s","--srna_file",help="output file")
args = parser.parse_args()

def main():
    proteins = []
    fh = open(args.input_file, "r")
    print("\t".join(["group", "product", "number", "GO_term_level_4(GO_name, p-value, percentage of CDS)",
                     "GO_term(for sRNA is \"start-end_strand_coverage(avg)\")", "CDS_positions", "gene_names"]))
    for row in csv.reader(fh, delimiter='\t'):
        if row[3] != "sRNA":
            proteins.append({"group": row[0], "name": row[1], "num": row[2],
                             "cut_go": row[3].replace("; ", ";"), "go": row[4],
                             "pos": row[5], "gene_name": row[-1]})
        else:
            coverages = []
            for info in row[4].split(";"):
                strand = info.split("_")[-1]
                start = info.split("-")[0]
                end = info.split("_")[0].split("-")[-1]
                sh = open(args.srna_file, "r")
                for srna in csv.reader(sh, delimiter='\t'):
                    if (srna[3] == start) and (
                        srna[4] == end) and (srna[5] == strand):
                        coverages.append(start + "-" + end + "_" + strand + "_" + srna[10])
                        break
                sh.close()
            proteins.append({"group": row[0], "name": row[1], "num": row[2],
                             "cut_go": "NA", "go": ";".join(coverages),
                             "pos": row[4], "gene_name": row[-1]})
    pre_group = None
    stats = {}
    for protein in proteins:
        if protein["cut_go"] != "NA":
            if pre_group is None:
                pre_group = protein["group"]
                stats[protein["group"]] = {}
                stats[protein["group"]]["total"] = 0
            else:
                if pre_group != protein["group"]:
                    stats[protein["group"]] = {}
                    stats[protein["group"]]["total"] = 0
                    pre_group = protein["group"]
            stats[protein["group"]]["total"] += 1
            for go in protein["cut_go"].split(";"):
                if go not in stats[protein["group"]].keys():
                    stats[protein["group"]][go] = 1
                else:
                    stats[protein["group"]][go] += 1
    for protein in proteins:
        if protein["cut_go"] != "NA":
            cut_gos = []
            for go in protein["cut_go"].split(";"):
                for stat, num in stats[protein["group"]].items():
                    if stat == go:
                        per = "{0:.0f}%".format(float(num)/float(stats[protein["group"]]["total"]) * 100)
                        cut_gos.append(go.replace(")", ", " + str(per) + ")"))
            print("\t".join([protein["group"], protein["name"], protein["num"],
                             ";".join(cut_gos), protein["go"], protein["pos"], protein["gene_name"]]))
        else:
            print("\t".join([protein["group"], protein["name"], protein["num"],
                             protein["cut_go"], protein["go"], protein["pos"], protein["gene_name"]]))

if __name__ == "__main__":
    main()
