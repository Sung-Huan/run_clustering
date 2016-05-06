import os
import sys
import csv
import argparse
import gff3

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-g","--go_file",help="input file")
parser.add_argument("-d","--go_obo",help="input file")
args = parser.parse_args()

def trace_back(go):
    roots = []
    detect = False
    with open(args.go_obo) as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                if detect:
                    break
            if line == "id: " + go:
                detect = True
            if detect:
                if line.startswith("is_a:"):
                    roots.append(line.split(" ")[1])
    return roots

def main():
    gh = open(args.go_file, "r")
    gos = []
    for row in csv.reader(gh, delimiter='\t'):
        if row[0] != "strain":
            gos.append({"protein": row[-2], "go": row[-1]})
    fh = open(args.input_file, "r")
    first = True
    num = 1
    for row in csv.reader(fh, delimiter='\t'):
        if len(row) == 1:
            if first:
                first = False
            else:
                out.close()
            out = open("test_" + str(num), "w")
            gofile = open("go_" + str(num), "w")
            num += 1
        if len(row) > 1:
            if row[2] == "CDS":
                infos = row[8].split(";")
                detect = False
                for info in infos:
                    if "protein_id" in info:
                        out.write(info.split("=")[-1] + "\n")
                        for go in gos:
                            if info.split("=")[-1] == go["protein"]:
                                detect = True
                                final_go = go
                    if "product" in info:
                        product = info.split("=")[-1]
                if detect:
                    all_gos = final_go["go"].split("; ")
                    run_gos = []
                    while (len(all_gos) != 0):
                        tmp_gos = []
#                        print(all_gos)
#                        print("VVV")
                        for go in all_gos:
#                            print(go)
                            if go not in run_gos:
                                run_gos.append(go)
                                roots = trace_back(go)
#                                print(roots)
                                for root in roots:
                                    if root not in all_gos:
                                        tmp_gos.append(root)
                        if len(tmp_gos) == 0:
                            break
                        else:
                            for tmp_go in tmp_gos:
                                if tmp_go not in all_gos:
                                    all_gos.append(tmp_go)
                    gofile.write("\t".join([final_go["protein"], product, ";".join(all_gos)]) + "\n")
        

if __name__ == "__main__":
    main()
