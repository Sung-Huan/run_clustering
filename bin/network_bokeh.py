#!/usr/bin/python

import os
import sys
import csv
import argparse
from scipy.stats.stats import pearsonr, spearmanr
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.spatial import distance
import networkx as nx
import pandas as pd
import bokeh as bk
from bokeh.plotting import figure, output_file, save
from bokeh.models import (
    Plot, 
    ColumnDataSource, 
    HoverTool, 
    TapTool, 
    Callback, 
    Circle, 
    MultiLine, 
    DataRange1d,
    BoxZoomTool,
    WheelZoomTool,
    ResetTool,
    ResizeTool,
    Text,
)

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
args = parser.parse_args()

def get_product(info):
    datas = info.split(";")
    for data in datas:
        if data.startswith("product"):
            product = data.split("=")[-1]
            return product
    for data in datas:
        if data.startswith("Name"):
            product = data.split("=")[-1]
            return product

def get_pro_id(info):
    datas = info.split(";")
    for data in datas:
        if data.startswith("protein_id"):
            product = data.split("=")[-1]
            return product
    return None

def plot(filename, q_info, q_exp, pngname, infos, exps, cutoff):
    G = nx.Graph()

    datas = []
    ids = []
    out.write(get_product(q_info) + "\n")
    datas.append(q_exp)
    ids.append(q_info)
    for info, exp in zip(infos, exps):
        if info != q_info:
            if "positive" in filename:
                if (float(spearmanr(exp, q_exp)[0]) >= cutoff):
                    out.write(get_product(info) + "\n")
                    datas.append(exp)
                    ids.append(info)
            elif "negative" in filename:
                if (float(spearmanr(exp, q_exp)[0]) <= cutoff):
                    out.write(get_product(info) + "\n")
                    datas.append(exp)
                    ids.append(info)
    out.close()
    fig = plt.figure()
    for data, id_ in zip(datas, ids):
        if "ID=srna" in id_:
            pass
        else:
            plt.plot(data, color="silver")
    for data, id_ in zip(datas, ids):
        if "ID=srna" in id_:
            if "Name=sRNA" in id_:
                plt.plot(data, color="lightblue", alpha=0.5)
            else:
                plt.plot(data, color="pink", alpha=0.5)
    plt.plot(datas[0], color="red")
    plt.savefig(pngname)

def graph_draw(g, pos, colors):
    labels = [ str(v) for v in g.nodes() ]
    colors = []
    tips = {"start": [], "end": [], "strand": [], "product": [], "link": []}
    for label in labels:
        if ("ID=srna" in label):
            colors.append("pink")
        else:
            colors.append("lightblue")
        product = get_product("_".join(label.split("_")[3:]))
        tips["start"].append(label.split("_")[0])
        tips["end"].append(label.split("_")[1])
        tips["strand"].append(label.split("_")[2])
        tips["product"].append(product)
        tips["link"].append("-")
    vx, vy = zip(*[ pos[v] for v in g.nodes() ])
    xs, ys = [], []
    lines = {"start": [], "end": [], "strand": [], "product": [], "link": []}
    for (a, b) in g.edges():
        x0, y0 = pos[a]
        x1, y1 = pos[b]
        p1 = get_product("_".join(a.split("_")[3:]))
        p2 = get_product("_".join(b.split("_")[3:]))
        lines["link"].append("<-->".join([p1, p2]))
        lines["start"].append("-")
        lines["end"].append("-")
        lines["product"].append("-")
        lines["strand"].append("-")
        xs.append([x0, x1])
        ys.append([y0, y1])
    hover = HoverTool(
        tooltips=[
            ("start", "@start"),
            ("end", "@end"),
            ("strand", "@strand"),
            ("product", "@product"),
            ("link", "@link")
        ]
    )
    f = figure(width=1000, height=1000, x_axis_type=None, y_axis_type=None,
               outline_line_color=None,
               tools=["pan,wheel_zoom,box_zoom,reset", hover])
    f.multi_line(xs, ys, line_color="black", source=lines)
    f.circle(vx, vy, size=8, line_color="black", fill_color=colors, source=tips, alpha=0.75)
    save(f, "network.html")

def main():
#    plt.figure(figsize=(25, 25))
#    plt.axis('off')
    genes = {}
    G = nx.Graph()
    exps = []
    infos = []
    fh = open(args.input_file, "r")
    for row in csv.reader(fh, delimiter='\t'):
        if (row[0] != "Orientation") and (row[0] == "sense"):
            exps.append([float(row[10]), float(row[11]), float(row[12]),
                         float(row[13]), float(row[14]), float(row[15]), float(row[16]),
                         float(row[17]), float(row[18]), float(row[19]),
                         float(row[20]), float(row[21]), float(row[22]), float(row[23])])
            infos.append("_".join([row[4], row[5], row[7], row[9]]))
    all_n = []
    all_e = []
    real_n = []
    real_colors = []
    datas = []
    all_nodes = []
    for info1, exp1 in zip(infos, exps):
        nodes = []
        edges = []
        colors = []
        num = 0
        product = get_product("_".join(info1.split("_")[3:]))
        name = "_".join(info1.split("_")[:3] + [product])
        genes[name] = 0
        for info2, exp2 in zip(infos, exps):
            if info1 != info2:
                if (float(spearmanr(exp1, exp2)[0]) >= 0.77):
                     num += 1
                     genes[name] += 1
                     if (info1 not in nodes):
                         if "ID=srna" in info1:
                             colors.append("lightblue")
                         else:
                             colors.append("red")
                         nodes.append(info1)
                         all_n.append(info1)
                     if (info2 not in nodes):
                         if "ID=srna" in info2:
                             colors.append("lightblue")
                         else:
                             colors.append("red")
                         nodes.append(info2)
                         all_n.append(info2)
                     edges.append([info1, info2])
                     all_e.append([info1, info2])
        if num >= 10:
#        real_n = real_n + nodes
#        real_colors = real_colors + colors
#        for node, color in zip(nodes, colors):
#            if node not in all_nodes:
#                product = get_product("_".join(node.split("_")[3:]))
#                datas.append({"gene": product, "start": node.split("_")[0],
#                              "end": node.split("_")[1], "color": color,
#                              "strand": node.split("_")[2]})
#                all_nodes.append(node)
            G.add_nodes_from(nodes)
            G.add_edges_from(edges)
    pos=nx.spring_layout(G)
    p = Plot(
    plot_height=600, 
    plot_width=800,
    )
    graph_draw(G, pos, datas)
    sort_genes = sorted(genes.items(), reverse=True, key=lambda value: value[1])
    out = open("sort_gene.csv", "w")
    for gene in sort_genes:
        out.write(str(gene) + "\n")
#
#
#
#    nodes = nx.draw_networkx_nodes(G,pos,node_size=40,node_shape='o', nodelist=real_n, node_color=real_colors, linewidths=0.5)
#    nx.draw_networkx_labels(G,pos,labels,font_size=60,font_color='r')
#    nodes.set_edgecolor('black')
#    nx.draw_networkx_edges(G,pos, width=0.5, len=10, edge_color='black')
#    plt.savefig("network.pdf")
#    plt.clf()

if __name__ == "__main__":
    main()
