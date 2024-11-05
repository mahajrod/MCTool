#!/usr/bin/env python
__author__ = 'tomarovsky'

import argparse
from collections import OrderedDict

import matplotlib.pyplot as plt
import pandas as pd

from MACE.Routines import Plotter

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, dest="input", type=lambda s: s.split(","),
                    help="Input file (or Comma-separated files) with two columns containing label in the first one and filename in the second."
                         "Stripped histograms will be drawn in the same order as labels")
parser.add_argument("--reference", action="store", dest="reference", type=lambda s: s.split(","),
                    help="Comma-separated list of species labels of output figure (for input1 and input2, respectively).")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=300,
                    help="Dpi of figure")

parser.add_argument("--colors", action="store", dest="colors", type=lambda s: s.split(","),
                    default=("#E04B4B", "#6094C3"), help="Color of double stripped histograms. Default: red, blue")
parser.add_argument("--palette", action="store", dest="palette", type=str,
                    default="turbo", help="Color palette of stripped histograms. Default: turbo")
parser.add_argument("--figure_height", action="store", dest="figure_height", type=float, default=6.0,
                    help="Height of figure in inches. Default: 6")
parser.add_argument("--figure_width_per_sample", action="store", dest="figure_width_per_sample", type=float, default=0.5,
                    help="Per sample width of figure in inches. Default: 1")
parser.add_argument("--figure_grid", action="store_true", default=False,
                    help="Add grid lines to the figure. Default: False")
parser.add_argument("--font-size", action="store", dest="font_size", type=float, default=16,
                    help="Font size. Default: 16")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=("png", "svg"), help="Comma-separated list of formats (supported by matplotlib) of "
                         "output figure.Default: png,svg")

parser.add_argument("-l", "--title", action="store", dest="title", default="",
                    help="Suptitle of figure. Default: ''")
parser.add_argument("--ylabel", action="store", dest="ylabel", default="Heterozygous SNPs/kbp",
                    help="y-axis label. Default: 'Heterozygous SNPs/kbp'")
parser.add_argument("-w", "--window_size", action="store", dest="window_size", required=True, type=float,
                    help="Size of the windows use for counts.")
parser.add_argument("-b", "--bin_width", action="store", dest="bin_width", required=True, default=0.1, type=float,
                    help="Width of bins in SNPs/1kbp")
parser.add_argument("-m", "--multiplicator", action="store", dest="multiplicator", default=1000, type=float,
                    help="Multiplicator for variant counts. "
                         "Default: 1000, i.e variant counts will be scaled to per 1 kbp ")
parser.add_argument("--ymin", action="store", dest="ymin", type=float, default=-0.1,
                    help="Minimum limit for Y axis . Default: -0.1")
parser.add_argument("--ymax", action="store", dest="ymax", type=float, default=None,
                    help="Maximum limit for Y axis. Default: not set")
parser.add_argument("--yticklist", action="store", dest="yticklist", type=lambda s: list(map(float, s.split(","))),
                    default=[0, 1, 2, 3, 4, 5],
                    help="Comma-separated tick list for Y axis. "
                         "Default: 0, 1, 2, 3, 4, 5")
parser.add_argument("--rotation", action="store", dest="rotation", type=float, default=45,
                    help="Rotation angle for X labels. Default: 45")
parser.add_argument("--horizontal_lines", action="store", dest="horizontal_lines",
                    type=lambda s: list(map(float, s.split(","))),
                    help="Comma-separated list of y-coordinates to draw horizontal lines. "
                         "Default: not set")

parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float,
                    help="Adjust left border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_top", action="store", dest="subplots_adjust_top", type=float,
                    help="Adjust top border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_right", action="store", dest="subplots_adjust_right", type=float,
                    help="Adjust right border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_bottom", action="store", dest="subplots_adjust_bottom", type=float,
                    help="Adjust bottom border of subplots on the figure. Default: matplotlib defaults")

parser.add_argument("--no_x", action="store_true", dest="no_x", default=False,
                    help="Do not use counts from X chromosome. Default: False")
parser.add_argument("--only_count", action="store_true", dest="only_count", default=False,
                    help="Only count and save variants to CSV, do not draw them. Default: False")

args = parser.parse_args()


def load_and_prepare_data(file_dict, label, args):
    data_list = []
    for entry in file_dict:
        df = pd.read_csv(file_dict[entry], sep="\t")
        if args.no_x:
            df = df[~df["CHROM"].str.contains("chrX")]
        df.set_index("CHROM", inplace=True)
        density = df.iloc[:, -1] / args.window_size * args.multiplicator
        data_list.extend([{'density': d, 'id': entry, 'Reference': label} for d in density])
    return data_list


file_dict = []
for file_name in args.input:
    with open(file_name, "r") as in_fd:
        file_dict.append(OrderedDict([line.strip().split("\t") for line in in_fd]))

if len(args.input) == 1:
    file_dict = file_dict[0]
    data = pd.DataFrame(load_and_prepare_data(file_dict, None, args))
else: # elif len(args.input) == 2
    data1 = load_and_prepare_data(file_dict[0], args.reference[0], args)
    data2 = load_and_prepare_data(file_dict[1], args.reference[1], args)
    data = pd.DataFrame(data1 + data2)

if args.only_count:
    data.to_csv(f"{args.input}.counts.csv", index=False)

fig, ax = plt.subplots(figsize=(args.figure_width_per_sample * len(data['id'].unique()), args.figure_height),
                       dpi=args.dpi)

if len(args.input) == 1:
    Plotter.draw_stripped_histograms(
        ax=ax,
        data=data,
        ymin=args.ymin,
        ymax=args.ymax,
        yticklist=args.yticklist,
        title=args.title,
        ylabel=args.ylabel,
        bin_width=args.bin_width,
        palette=args.palette,
        rotation=args.rotation,
        horizontal_lines=args.horizontal_lines,
        figure_grid=args.figure_grid,
        font_size=args.font_size,
    )
else: # elif len(args.input) == 2
    Plotter.draw_double_stripped_histograms(
        ax=ax,
        data=data,
        ymin=args.ymin,
        ymax=args.ymax,
        yticklist=args.yticklist,
        title=args.title,
        ylabel=args.ylabel,
        bin_width=args.bin_width,
        references=args.reference,
        colors=args.colors,
        rotation=args.rotation,
        horizontal_lines=args.horizontal_lines,
        figure_grid=args.figure_grid,
        font_size=args.font_size,
    )

plt.subplots_adjust(top=args.subplots_adjust_top, bottom=args.subplots_adjust_bottom,
                    left=args.subplots_adjust_left, right=args.subplots_adjust_right)

plt.tight_layout()

for ext in "png", "svg":
    plt.savefig("{0}.{1}".format(args.output_prefix, ext), transparent=False)



