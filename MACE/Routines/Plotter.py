#!/usr/bin/env python
__author__ = "tomarovsky"

import distinctipy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Rectangle


class Plotter:
    def __init__(self):
        pass

    def histogram_scaling(self, histogram, scale=0.45):
        scaling_factor = scale / max(histogram, default=1)
        histogram *= scaling_factor
        return histogram

    def histogram_with_double_first_column_and_bin_centers(self, df, bins):
        hist, bin_edges = np.histogram(df, bins=bins, density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_centers = np.insert(bin_centers, 0, 0)
        hist = np.insert(hist, 0, hist[0])
        bin_centers = np.insert(bin_centers, 0, 0)
        hist = np.insert(hist, 0, 0)
        bin_centers = np.append(bin_centers, bin_centers[-1])
        hist = np.append(hist, 0)
        hist = self.histogram_scaling(hist)
        return hist, bin_centers

    def draw_stripped_histograms(
        self,
        ax,
        data,
        ymin,
        ymax,
        yticklist,
        title="",
        ylabel="Heterozygous SNPs/kbp",
        bin_width=0.1,
        palette="turbo",
        rotation=45,
        horizontal_lines=[],
        figure_grid=True,
        font_size=None,
    ):

        ax.spines[["right", "top"]].set_visible(False)
        ax.set_axisbelow(True)

        if font_size is not None:
            plt.rcParams.update({"font.size": font_size})

        unique_ids = data["id"].unique()
        colors = sns.color_palette(palette, len(unique_ids))

        for i, unique_id in enumerate(unique_ids):
            df = data[data["id"] == unique_id]["density"]

            # bins
            bins = np.arange(df.min(), df.max() + bin_width, bin_width)

            # draw stripped histograms
            hist, bin_centers = self.histogram_with_double_first_column_and_bin_centers(df, bins)
            ax.fill_betweenx(bin_centers, i, i - hist, color=colors[i], edgecolor=colors[i])
            ax.fill_betweenx(bin_centers, i, i + hist, color=colors[i], edgecolor=colors[i])

            # boxplot
            ax.boxplot(
                df,
                vert=True,
                positions=[i],
                widths=0.15,
                whis=0.1,
                showfliers=True,
                showmeans=True,
                patch_artist=False,
                meanprops=dict(marker="o", markerfacecolor="white", markeredgecolor="white", markersize=0.5),
                flierprops=dict(marker="o", markerfacecolor="#575757", markeredgecolor="#575757", markersize=0.5),
                medianprops=dict(color="#575757", solid_capstyle="butt"),
                boxprops=dict(color="#575757", linewidth=0.5),
                whiskerprops=dict(color="#575757", linewidth=0.5),
                capprops=dict(color="#575757", linewidth=0.5),
            )

        ax.set_xticks(range(len(unique_ids)))
        ax.set_xticklabels(unique_ids, ha="right", rotation=rotation)
        ax.set_yticks(yticklist)
        if horizontal_lines:
            for ycoord in horizontal_lines:
                ax.axhline(y=ycoord, color="red", linestyle="--", linewidth=1)
        ax.set_ylim(ymin=ymin, ymax=ymax)
        ax.set_xlim(xmin=-0.5, xmax=len(unique_ids) - 0.5)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if figure_grid:
            ax.grid(axis="y", linestyle="--", alpha=0.5)

    def draw_double_stripped_histograms(
        self,
        ax,
        data,
        ymin,
        ymax,
        yticklist,
        title="",
        ylabel="Heterozygous SNPs/kbp",
        bin_width=0.1,
        references=["Species 1", "Species 2"],
        colors=["#E04B4B", "#6094C3"],
        rotation=45,
        horizontal_lines=[],
        figure_grid=True,
        font_size=None,
    ):

        ax.spines[["right", "top"]].set_visible(False)
        ax.set_axisbelow(True)

        if font_size is not None:
            plt.rcParams.update({"font.size": font_size})

        unique_ids = data["id"].unique()

        for i, unique_id in enumerate(unique_ids):
            df_1 = data[(data["Reference"] == references[0]) & (data["id"] == unique_id)]["density"]
            df_2 = data[(data["Reference"] == references[1]) & (data["id"] == unique_id)]["density"]

            # bins
            bins = np.arange(min([df_1.min(), df_2.min()]), max([df_1.max(), df_2.max()]) + bin_width, bin_width)

            # stripped histograms
            hist_1, bin_centers_1 = self.histogram_with_double_first_column_and_bin_centers(df_1, bins)
            hist_2, bin_centers_2 = self.histogram_with_double_first_column_and_bin_centers(df_2, bins)
            ax.fill_betweenx(
                bin_centers_1, i, i - hist_1, color=colors[0], edgecolor=colors[0], linewidth=0.5, label=r"$\mathit{" + references[0] + "}$"
            )
            ax.fill_betweenx(
                bin_centers_2, i, i + hist_2, color=colors[1], edgecolor=colors[1], linewidth=0.5, label=r"$\mathit{" + references[1] + "}$"
            )
            ax.plot([i, i], [0, max([len(hist_1) * bin_width, len(hist_2) * bin_width])], linewidth=0.6, color="#575757")

            # boxplot
            for d, pos in zip([df_1, df_2], [i - 0.02, i + 0.02]):
                ax.boxplot(
                    d,
                    vert=True,
                    positions=[pos],
                    widths=0.04,
                    whis=0,
                    showfliers=False,
                    showmeans=False,
                    patch_artist=True,
                    boxprops=dict(facecolor="#575757", color="#575757", linewidth=0),
                    medianprops=dict(color="white", solid_capstyle="butt"),
                    whiskerprops=dict(color="#575757", linewidth=0),
                    capprops=dict(color="#575757", linewidth=0),
                )

            ax.legend(title="Reference", loc="upper left", ncols=2) if i == 0 else None

        ax.set_xticks(range(len(unique_ids)))
        ax.set_xticklabels(unique_ids, ha="right", rotation=rotation)
        ax.set_yticks(yticklist)
        if horizontal_lines:
            for ycoord in horizontal_lines:
                ax.axhline(y=ycoord, color="red", linestyle="--", linewidth=1)
        ax.set_ylim(ymax=ymax, ymin=ymin)
        ax.set_xlim(-0.5, len(unique_ids) - 0.5)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if figure_grid:
            ax.grid(axis="y", linestyle="--", alpha=0.5)

    def half_patch(self, color1, color2):
        class SplitPatchHandler(HandlerPatch):
            def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
                p1 = Rectangle((xdescent, ydescent), width / 2, height, facecolor=color1, transform=trans)
                p2 = Rectangle((xdescent + width / 2, ydescent), width / 2, height, facecolor=color2, transform=trans)
                return [p1, p2]

        return Rectangle((0, 0), 1, 1), SplitPatchHandler()

    def draw_double_admixture_barplot(
        self,
        ax,
        df,
        global_columns,
        local_columns,
        references,
        local_colors=["#e02828", "#4286c3"],
        global_colors=["#e06c6c", "#7da2c3"],
        rotation=45,
        yticks=[0, 25, 50, 75, 100],
        xlabel="",
        legend_position=(1.006, 1.2),
        font_size=None,
    ):

        ax.spines[["right", "top"]].set_visible(False)
        ax.set_axisbelow(True)

        if font_size is not None:
            plt.rcParams.update({"font.size": font_size})

        index_column = df.columns[0]  # first column as index ('Sample_ID')

        df.set_index(index_column, inplace=True)

        stacked_data, stacked_data2 = df[global_columns], df[local_columns]

        stacked_data2.plot(kind="bar", stacked=True, width=0.46, color=local_colors, ax=ax, position=0, legend=False)
        stacked_data.plot(kind="bar", stacked=True, width=0.46, color=global_colors, ax=ax, position=1, legend=False)

        # Custom legend
        species_1 = self.half_patch(global_colors[0], local_colors[0])
        species_2 = self.half_patch(global_colors[1], local_colors[1])
        ax.legend(
            [species_1[0], species_2[0]],
            [f"Global and local ADMIXTURE from $\\mathit{{{references[0]}}}$", f"Global and local ADMIXTURE from $\\mathit{{{references[1]}}}$"],
            bbox_to_anchor=legend_position,
            handler_map={species_1[0]: species_1[1], species_2[0]: species_2[1]},
        )

        ax.set_xlim(right=len(stacked_data) - 0.5)
        ax.set_xticklabels(df.index, ha="right", rotation=rotation)
        ax.set_xlabel(xlabel)
        ax.set_yticks(yticks)

    def draw_roh_cumsum_plot(
        self,
        ax,
        roh_files,
        genome_length,
        colors=None,
        xlim=(10e4, 150e6),
        ylim=(0, 0.35),
        xlabel="Homozygous tract length in base-pairs",
        ylabel="Cumulative genome fraction in ROHs",
        legend_title="",
        legend_loc="upper left",
        legend_ncol=2,
        linewidth=1.5,
        figure_grid=True,
        font_size=None,
    ):

        ax.spines[["right", "top"]].set_visible(False)
        ax.set_axisbelow(True)

        if font_size is not None:
            plt.rcParams.update({"font.size": font_size})

        if colors is None:
            colors = distinctipy.get_colors(len(roh_files))
        else:
            if type(colors) == str:
                colors = sns.color_palette(colors, len(roh_files))

        for count, f in enumerate(roh_files):
            df = pd.read_csv(f, sep="\t")
            df.sort_values(by=["length"], inplace=True, ignore_index=True)

            df["cumsum"] = np.cumsum(df["length"]) / genome_length
            df.loc[len(df.index)] = [None, None, None, xlim[1], df.iloc[-1]["cumsum"]]

            ax.step(df["length"], df["cumsum"], label=f.split("/")[-1].split(".")[0], linewidth=linewidth, c=colors[count], where="post")

        ax.set_xscale("log")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(title=r"$\mathit{" + legend_title + "}$", loc=legend_loc, ncol=legend_ncol)

        if figure_grid:
            ax.grid(True, linestyle="--", alpha=0.5)

    def draw_busco_summary_plot(
        self,
        ax,
        data,
        colors=["#23b4e8", "#008dbf", "#fbbc04", "#ea4335"],
        legend_loc=(-0.005, 0.97),
        legend_ncol=4,
        xticks=[0, 25, 50, 75, 100],
    ):

        ax.spines[["right", "top"]].set_visible(False)
        ax.set_axisbelow(True)

        position = range(len(data))

        ax.barh(position, data["S"], height=0.8, label="Complete and single-copy BUSCOs (S)", color=colors[0])
        ax.barh(position, data["D"], height=0.8, left=data["S"], label="Complete and duplicated BUSCOs (D)", color=colors[1])
        ax.barh(position, data["F"], height=0.8, left=data["S"] + data["D"], label="Fragmented BUSCOs (F)", color=colors[2])
        ax.barh(position, data["M"], height=0.8, left=data["S"] + data["D"] + data["F"], label="Missing BUSCOs (M)", color=colors[3])

        ax.legend(ncol=legend_ncol, loc=legend_loc, handlelength=0.8, frameon=False)

        ax.set_yticks(range(len(data["Species"].to_list())))
        ax.set_yticklabels(["" for _ in range(len(data["Species"].to_list()))])
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{i}%" for i in xticks])
        ax.set_xlim(xticks[0], xticks[-1])

        for i, row in data.iterrows():
            ax.text(xticks[0] - 1, i, f'{row["Species"]}', va="center", ha="right", fontweight="medium", style="italic")
            ax.text(
                xticks[0] + 1,
                i,
                f"S: {row['S']}%   |   D: {row['D']}%   |   F: {row['F']}%   |   M: {row['M']}%",
                va="center",
                ha="left",
                fontweight="bold",
                color="white",
            )



