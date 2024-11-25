#!/usr/bin/env python
__author__ = "tomarovsky"

import os
import re
from copy import deepcopy
from functools import partial
from pathlib import Path

import distinctipy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Rectangle
from RouToolPa.Collections.General import SynDict
from RouToolPa.Parsers.BED import CollectionBED
from RouToolPa.Parsers.BLAST import CollectionBLAST
from RouToolPa.Parsers.GFF import CollectionGFF
from RouToolPa.Parsers.STR import CollectionSTR
from RouToolPa.Parsers.VCF import CollectionVCF

from MACE.Routines import StatsVCF, Visualization

# TODO: add PCA plot from PRINK
# TODO: add simple plot (PAR)
# TODO: add function to remove empty plot
# TODO: add ete3 plots


class Plotter:
    def __init__(self):
        pass

    def set_paperticks_style(self, font_scale):
        """
        Configures a "ticks" style and "paper" context.
        """
        custom_params = {
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.grid": True,
            "axes.axisbelow": True,
            "grid.color": "#dfdfdf",
            "grid.linestyle": "--",
        }
        sns.set_theme(style="ticks", rc=custom_params)
        sns.set_context("paper", font_scale=font_scale)

    def set_figure_fontsize(self, font_scale):
        """
        plt.rcParams.update({'font.size': font_scale})
        """
        plt.rcParams.update({"font.size": font_scale})

    def annotate_subplot(self, ax, annotation, offset=(-0.1, 1.1), fontsize=12, fontweight="bold", color="black"):
        """
        Add an annotation to a specific subplot (Axes) in a figure.

        Parameters:
        - ax : matplotlib.axes.Axes
            The specific subplot (Axes) to label.
        - annotation : str
            The annotation to add (e.g., "A", "B", etc.).
        - offset : tuple of float, optional
            The x and y offsets for the label relative to the subplot in Axes coordinates.
        - fontsize : int, optional
            Font size for the label.
        - fontweight : str or int, optional
            Font weight for the label (e.g., "bold", "normal", or a numeric value).
        - color : str, optional
            The color used for the sign. Defaults to "black".
        """
        ax.text(
            offset[0], offset[1], annotation, transform=ax.transAxes, fontsize=fontsize, fontweight=fontweight, color=color, va="center", ha="center"
        )

    def scaled_histogram_with_extended_bins(self, df, bins, scale=0.45):
        hist, bin_edges = np.histogram(df, bins=bins, density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_centers = np.insert(bin_centers, 0, 0)
        hist = np.insert(hist, 0, hist[0])
        bin_centers = np.insert(bin_centers, 0, 0)
        hist = np.insert(hist, 0, 0)
        bin_centers = np.append(bin_centers, bin_centers[-1])
        hist = np.append(hist, 0)
        scaling_factor = scale / max(hist, default=1)
        hist *= scaling_factor
        return hist, bin_centers

    def add_boxplot_classic(self, ax, df, position, width=0.15):
        ax.boxplot(
            df,
            vert=True,
            positions=[position],
            widths=width,
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

    def add_boxplot_empty(self, ax, df, position, width=0.15):
        ax.boxplot(
            df,
            vert=True,
            positions=[position],
            widths=width,
            showfliers=False,
            showmeans=True,
            showcaps=False,
            patch_artist=True,
            boxprops=dict(alpha=0),  # color="#575757", facecolor="#575757",
            flierprops=dict(marker="o", markerfacecolor="#575757", markeredgecolor="#575757", markersize=0.5),
            whiskerprops=dict(color="#575757", linewidth=1),
            capprops=dict(color="#575757", linewidth=1),
            meanprops=dict(marker="o", markerfacecolor="w", markeredgecolor="w", markersize=2),
            medianprops=dict(color="w", linewidth=1),
        )
        q25, q75 = df.quantile([0.25, 0.75])
        ax.hlines(y=q25, xmin=position - 0.1, xmax=position + 0.1, colors="#575757", linewidth=1)
        ax.hlines(y=q75, xmin=position - 0.1, xmax=position + 0.1, colors="#575757", linewidth=1)

    def add_double_boxplot(self, ax, df_1, df_2, position, width=0.04):
        for d, pos in zip([df_1, df_2], [position - 0.02, position + 0.02]):
            ax.boxplot(
                d,
                vert=True,
                positions=[pos],
                widths=width,
                whis=0,
                showfliers=False,
                showmeans=False,
                patch_artist=True,
                boxprops=dict(facecolor="#575757", color="#575757", linewidth=0),
                medianprops=dict(color="white", solid_capstyle="butt", linewidth=2),
                whiskerprops=dict(color="#575757", linewidth=0),
                capprops=dict(color="#575757", linewidth=0),
            )

    def process_variant_counts(self, file_paths, removed_chrX=None, reference=None, window_size=1, multiplicator=1):
        data_list = []
        for count in file_paths:
            df = pd.read_csv(count, sep="\t")
            file_path = Path(count)
            file_name = file_path.stem
            id = file_name.split(".")[0]

            if removed_chrX:
                # If removed_chrX is a list, remove the specified chromosomes; otherwise, remove a single string
                if isinstance(removed_chrX, list):
                    for chrX in removed_chrX:
                        if chrX is not None:
                            df = df[~df["CHROM"].str.contains(chrX)]
                else:
                    df = df[~df["CHROM"].str.contains(removed_chrX)]

            df.set_index("CHROM", inplace=True)
            density = df.iloc[:, -1] / window_size * multiplicator

            if reference:
                data_list.extend([{"density": d, "id": id, "Reference": reference} for d in density])
            else:
                data_list.extend([{"density": d, "id": id} for d in density])

        return data_list

    def draw_stripped_histograms(
        self,
        ax,
        data,
        ymin,
        ymax,
        yticklist,
        window_size=1000000,
        multiplicator=1000,
        removed_chrX="",
        title="",
        ylabel="Heterozygous SNPs/kbp",
        bin_width=0.1,
        boxplot_type="empty",
        boxplot_width=0.15,
        palette="turbo",
        rotation=45,
        statuses=dict(),
        sample_ids=[],  # same sorting as in data
        horizontal_lines=[],
        figure_grid=True,
        font_size=None,
        sort_by=None,  # 'mean' or 'median'
    ):
        """
        Draws stripped histograms and boxplots to visualize heterozygous SNP density for multiple samples.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the histograms.

        data : list of str
            A list of file paths to tab-separated files containing SNP density data.

        ymin : float
            The minimum value for the y-axis.

        ymax : float
            The maximum value for the y-axis.

        yticklist : list of float
            List of y-axis tick positions.

        window_size : int, optional
            The size of the window used for calculating density. Defaults to 1,000,000.

        multiplicator : float, optional
            A multiplier applied to the density values for scaling. Defaults to 1000.

        removed_chrX : str, optional
            Chromosome name to be removed from the analysis (e.g., "chrX"). If not provided, no chromosomes are removed.
            Defaults to an empty string.

        title : str, optional
            The title of the plot. Defaults to an empty string.

        ylabel : str, optional
            The label for the y-axis. Defaults to "Heterozygous SNPs/kbp".

        bin_width : float, optional
            The width of the bins for the histograms. Defaults to 0.1.

        boxplot_type : str, optional
            Type of boxplot to draw. Options are 'classic' or 'empty'. Defaults to 'empty'.

        boxplot_width : float, optional
            The width of the boxplot elements. Defaults to 0.15.

        palette : str, optional
            The color palette used for the histograms. Defaults to "turbo".

        rotation : int or float, optional
            The rotation angle for x-axis labels. Defaults to 45.

        statuses : dict, optional
            A dictionary mapping sample IDs to status labels. Status labels are displayed above the histograms.
            Defaults to an empty dictionary.

        sample_ids : list of str, optional
            A list of sample IDs for labeling the x-axis. If not provided, the unique IDs from the data are used.
            Defaults to an empty list.

        horizontal_lines : list of float, optional
            List of y-coordinates at which to draw horizontal lines across the plot. Defaults to an empty list.

        figure_grid : bool, optional
            Whether to display a grid on the y-axis. Defaults to True.

        font_size : int or float, optional
            Font size for the plot text. If not provided, the default font size is used.

        sort_by : str, optional
            Specifies how to sort the sample IDs. Options are 'mean' or 'median'. If not provided, no sorting is applied.
            Defaults to None.

        Notes:
        ------
        - The function reads density data from each file, optionally filters out a specified chromosome, and calculates
          SNP density. It then creates stripped histograms and optional boxplots for each sample.
        - Status labels are displayed above the histograms if provided in the `statuses` dictionary.
        """

        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)

        ax.set_axisbelow(True)

        if font_size is not None:
            plt.rcParams.update({"font.size": font_size})

        data_list = self.process_variant_counts(data, removed_chrX=removed_chrX, window_size=window_size, multiplicator=multiplicator)
        data = pd.DataFrame(data_list)

        # Calculate mean or median values for sorting if needed
        if sort_by in ["mean", "median"]:
            sort_values = data.groupby("id")["density"].agg(sort_by).sort_values()
            unique_ids = sort_values.index
        else:
            unique_ids = data["id"].unique()

        colors = sns.color_palette(palette, len(unique_ids))

        for i, unique_id in enumerate(unique_ids):
            df = data[data["id"] == unique_id]["density"]

            # bins
            bins = np.arange(df.min(), df.max() + bin_width, bin_width)

            # draw stripped histograms
            hist, bin_centers = self.scaled_histogram_with_extended_bins(df, bins)
            ax.fill_betweenx(bin_centers, i, i - hist, color=colors[i], edgecolor=colors[i])
            ax.fill_betweenx(bin_centers, i, i + hist, color=colors[i], edgecolor=colors[i])

            # boxplot
            if boxplot_type == "classic":
                self.add_boxplot_classic(ax, df, i, width=boxplot_width)
            elif boxplot_type == "empty":
                self.add_boxplot_empty(ax, df, i, width=boxplot_width)

            # add statuses
            if statuses and unique_id in statuses:
                last_bin_center = bin_centers[-1]
                ax.text(i, last_bin_center + 0.05, statuses[unique_id], ha="center", va="bottom", color=colors[i])

        ax.set_xticks(range(len(unique_ids)))
        if sample_ids:
            ax.set_xticklabels(sample_ids, ha="right", rotation=rotation)
        else:
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
        left_data,
        right_data,
        ymin,
        ymax,
        yticklist,
        window_size=1000000,
        multiplicator=1000,
        removed_chrX=[],
        sort_by=None,
        title="",
        ylabel="Heterozygous SNPs/kbp",
        bin_width=0.1,
        boxplot_width=0.04,
        references=["Species 1", "Species 2"],
        colors=["#E04B4B", "#6094C3"],
        rotation=45,
        horizontal_lines=[],
        figure_grid=True,
        font_size=None,
    ):
        """
        Draw double stripped histograms comparing the density of data between two species across various IDs.

        This function plots histograms and boxplots comparing heterozygous SNP density from two different data sets,
        represented as `left_data` and `right_data`. The histograms are stripped and filled to visually compare
        the density distributions for each ID. The function also includes the ability to sort IDs by mean or
        median density values and customize various plot elements.

        Parameters
        ----------
        ax : matplotlib.axes._axes.Axes
            The axes object where the histograms and boxplots will be drawn.
        left_data : list of str
            List of file paths to data sets representing the left side (first species).
        right_data : list of str
            List of file paths to data sets representing the right side (second species).
        ymin : float
            Minimum value for the y-axis.
        ymax : float
            Maximum value for the y-axis.
        yticklist : list of float
            List of y-tick values to use on the y-axis.
        window_size : int, optional, default=1000000
            The size of the window used to normalize density values.
        multiplicator : int, optional, default=1000
            A multiplier applied to the density values for scaling purposes.
        removed_chrX : str, optional, default=""
            String pattern used to remove certain chromosomes from the data if needed.
        title : str, optional, default=""
            Title of the plot.
        ylabel : str, optional, default="Heterozygous SNPs/kbp"
            Label for the y-axis.
        bin_width : float, optional, default=0.1
            Width of the histogram bins.
        boxplot_width : float, optional, default=0.04
            Width of the boxplots.
        references : list of str, optional, default=["Species 1", "Species 2"]
            List of labels for the two data sets being compared.
        colors : list of str, optional, default=["#E04B4B", "#6094C3"]
            List of colors to use for the histograms of the two data sets.
        rotation : int or float, optional, default=45
            Angle of rotation for the x-axis labels.
        horizontal_lines : list of float, optional, default=[]
            List of y-coordinates for horizontal lines to be drawn across the plot.
        figure_grid : bool, optional, default=True
            Whether to display a grid on the y-axis.
        font_size : int, optional, default=None
            Font size for the plot text. If None, the default font size is used.
        sort_by : str, optional, default="mean"
            Method to sort the IDs, either 'mean' or 'median'.
            - 'mean': Sort IDs by the mean density value of the left data.
            - 'median': Sort IDs by the median density value of the left data.

        Notes
        -----
        - The function reads data from the provided file paths, calculates density values, and removes any data
          matching `removed_chrX` if specified.
        - It uses the `scaled_histogram_with_extended_bins` method to generate histograms and the
          `add_double_boxplot` method to add boxplots to the plot.
        - The x-axis labels are the unique IDs, and the y-axis limits, grid, and additional horizontal lines
          can be customized.
        """

        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)

        ax.set_axisbelow(True)

        if font_size is not None:
            plt.rcParams.update({"font.size": font_size})

        left_data_list = self.process_variant_counts(
            left_data,
            removed_chrX=[removed_chrX[0]] if removed_chrX else None,
            reference=references[0],
            window_size=window_size,
            multiplicator=multiplicator,
        )
        right_data_list = self.process_variant_counts(
            right_data,
            removed_chrX=[removed_chrX[1]] if removed_chrX else None,
            reference=references[1],
            window_size=window_size,
            multiplicator=multiplicator,
        )

        data = pd.DataFrame(left_data_list + right_data_list)

        # sorting
        if sort_by == "mean":
            sorting_values = data[data["Reference"] == references[0]].groupby("id")["density"].mean()
            unique_ids = sorting_values.sort_values().index
        elif sort_by == "median":
            sorting_values = data[data["Reference"] == references[0]].groupby("id")["density"].median()
            unique_ids = sorting_values.sort_values().index
        else:
            unique_ids = data["id"].unique()

        for i, unique_id in enumerate(unique_ids):
            df_1 = data[(data["Reference"] == references[0]) & (data["id"] == unique_id)]["density"]
            df_2 = data[(data["Reference"] == references[1]) & (data["id"] == unique_id)]["density"]

            # bins
            bins = np.arange(min([df_1.min(), df_2.min()]), max([df_1.max(), df_2.max()]) + bin_width, bin_width)

            # stripped histograms
            hist_1, bin_centers_1 = self.scaled_histogram_with_extended_bins(df_1, bins)
            hist_2, bin_centers_2 = self.scaled_histogram_with_extended_bins(df_2, bins)
            ax.fill_betweenx(
                bin_centers_1, i, i - hist_1, color=colors[0], edgecolor=colors[0], linewidth=0.5, label=r"$\mathit{" + references[0] + "}$"
            )
            ax.fill_betweenx(
                bin_centers_2, i, i + hist_2, color=colors[1], edgecolor=colors[1], linewidth=0.5, label=r"$\mathit{" + references[1] + "}$"
            )
            ax.plot([i, i], [0, max([len(hist_1) * bin_width, len(hist_2) * bin_width])], linewidth=0.6, color="#575757")

            # boxplot
            self.add_double_boxplot(ax, df_1, df_2, i, width=boxplot_width)

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

    def legend_half_patch(self, color1, color2):
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
        """
        Draw a double barplot to visualize global and local ADMIXTURE proportions for multiple samples.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the plot.

        df : pandas.DataFrame
            A DataFrame containing ADMIXTURE data. The first column must be `Sample_ID`,
            followed by columns with global and local ADMIXTURE proportions.

        global_columns : list of str
            Column names in `df` corresponding to global ADMIXTURE proportions.

        local_columns : list of str
            Column names in `df` corresponding to local ADMIXTURE proportions.

        references : list of str
            Names of the reference populations used in ADMIXTURE analysis. This should be a list of two items.

        local_colors : list of str, optional
            Colors for the local ADMIXTURE bars. Defaults to ["#e02828", "#4286c3"].

        global_colors : list of str, optional
            Colors for the global ADMIXTURE bars. Defaults to ["#e06c6c", "#7da2c3"].

        rotation : int, optional
            Rotation angle for the x-axis labels. Defaults to 45 degrees.

        yticks : list of int, optional
            Y-axis tick values. Defaults to [0, 25, 50, 75, 100].

        xlabel : str, optional
            Label for the x-axis. Defaults to an empty string.

        legend_position : tuple of float, optional
            Position of the legend box. Defaults to (1.006, 1.2).

        font_size : int, optional
            Font size for the plot text. If `None`, the default font size is used.

        Notes:
        ------
        - The first column of `df` is used as the index (sample identifiers) for the plot.
        - Two stacked bar plots are created side-by-side: one for global ADMIXTURE proportions
          and the other for local ADMIXTURE proportions.
        - A custom legend is generated to distinguish between global and local ADMIXTURE for
          each reference population.
        """

        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)
        ax.set_axisbelow(True)

        if font_size is not None:
            plt.rcParams.update({"font.size": font_size})

        index_column = df.columns[0]  # first column as index ('Sample_ID')

        df.set_index(index_column, inplace=True)

        stacked_data, stacked_data2 = df[global_columns], df[local_columns]

        stacked_data2.plot(kind="bar", stacked=True, width=0.46, color=local_colors, ax=ax, position=0, legend=False)
        stacked_data.plot(kind="bar", stacked=True, width=0.46, color=global_colors, ax=ax, position=1, legend=False)

        # Custom legend
        species_1 = self.legend_half_patch(global_colors[0], local_colors[0])
        species_2 = self.legend_half_patch(global_colors[1], local_colors[1])
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
        data,
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
    ):
        """
        Draws a cumulative sum plot of runs of homozygosity (ROH) for multiple samples.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the plot.

        data : list of str
            A list of file paths to tab-separated files containing ROH data. Each file should have a column "length"
            indicating the length of homozygous tracts.

        genome_length : float
            The total length of the genome, used to normalize cumulative ROH values.

        colors : list of str or str, optional
            The colors used for the plot lines. If a list of colors is provided, it should match the number of samples.
            If a string is provided, it specifies a color palette name from seaborn. If None, distinct colors are
            automatically generated. Defaults to None.

        xlim : tuple of float, optional
            Limits for the x-axis (log scale), representing the range of tract lengths to display. Defaults to (10e4, 150e6).

        ylim : tuple of float, optional
            Limits for the y-axis, representing the range of cumulative genome fractions. Defaults to (0, 0.35).

        xlabel : str, optional
            The label for the x-axis. Defaults to "Homozygous tract length in base-pairs".

        ylabel : str, optional
            The label for the y-axis. Defaults to "Cumulative genome fraction in ROHs".

        legend_title : str, optional
            The title for the legend. Defaults to an empty string.

        legend_loc : str, optional
            The location of the legend box within the plot. Defaults to "upper left".

        legend_ncol : int, optional
            Number of columns in the legend. Defaults to 2.

        linewidth : float, optional
            The width of the plot lines. Defaults to 1.5.

        figure_grid : bool, optional
            Whether to display a grid on the plot. Defaults to True.

        Notes:
        ------
        - The function reads ROH data from each file, sorts it by tract length, and calculates the cumulative sum of
          tract lengths normalized by the genome length.
        - The plot uses a logarithmic x-axis to display a wide range of homozygous tract lengths.
        - A step plot is created for each sample, and a legend is added to differentiate between samples.
        """
        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)
        ax.set_axisbelow(True)

        if colors is None:
            colors = distinctipy.get_colors(len(data))
        else:
            if type(colors) == str:
                colors = sns.color_palette(colors, len(data))

        for count, f in enumerate(data):
            df = pd.read_csv(f, sep="\t", header=None, names=["scaffold", "start", "end", "length"])
            df.sort_values(by=["length"], inplace=True, ignore_index=True)

            df["cumsum"] = np.cumsum(df["length"]) / genome_length

            label = f.split("/")[-1].split(".")[0]
            print(f"{label}\t{len(df.index)}\t{df['length'].sum()}\t{df['length'].sum() / genome_length * 100}")

            df.loc[len(df.index)] = [None, None, None, xlim[1], df.iloc[-1]["cumsum"]]

            ax.step(df["length"], df["cumsum"], label=label, linewidth=linewidth, c=colors[count], where="post")

        ax.set_xscale("log")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(title=r"$\mathit{" + legend_title.replace(" ", r"\,") + "}$", loc=legend_loc, ncol=legend_ncol)

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
        bold_species_indices=None,
    ):
        """
        Draws a horizontal bar plot to visualize BUSCO (Benchmarking Universal Single-Copy Orthologs) results
        for multiple species.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the plot.

        data : str or list of str
            Path to a tab-separated file with columns [file_path, species_name, species_id]
            or a list of paths to BUSCO short summary text files (short_summary_{species}.txt).

        colors : list of str, optional
            A list of four color hex codes used for the bar segments:
            - colors[0]: Color for "Complete and single-copy BUSCOs (S)"
            - colors[1]: Color for "Complete and duplicated BUSCOs (D)"
            - colors[2]: Color for "Fragmented BUSCOs (F)"
            - colors[3]: Color for "Missing BUSCOs (M)"
            Default is ["#23b4e8", "#008dbf", "#fbbc04", "#ea4335"].

        legend_loc : tuple of float, optional
            Position of the legend box within the plot. Defaults to (-0.005, 0.97).

        legend_ncol : int, optional
            Number of columns in the legend. Defaults to 4.

        xticks : list of int, optional
            X-axis tick positions, representing percentages. Defaults to [0, 25, 50, 75, 100].

        bold_species_indices : list of int, optional
            A list of indices (0-based) of species that should be labeled in bold. If specified,
            these species will have their labels printed in bold font. Indices are counted from the
            top to bottom, but can also be provided in reverse order (from bottom to top).
            For example, if `bold_species_indices=[1, 2, 3]`, the species at positions 1, 2, and 3 from the bottom
            will be labeled in bold.
        """
        # Check if data is a TSV file or a list of file paths
        if isinstance(data, str):
            df_input = pd.read_csv(data, sep="\t", header=None, names=["file_path", "Species", "ID"])
            df_input["Species"] = df_input["Species"].str.replace("_", " ")
            file_paths = df_input["file_path"].tolist()
        else:
            file_paths = data
            df_input = pd.DataFrame(
                {
                    "file_path": file_paths,
                    "Species": [file_path.split("/")[-1][14:-4].replace("_", " ") for file_path in file_paths],
                    "ID": [""] * len(file_paths),  # Placeholder if IDs are not provided
                }
            )

        # Extract data from each BUSCO summary file
        data_list = []
        for file_path, species, species_id in zip(df_input["file_path"], df_input["Species"], df_input["ID"]):
            with open(file_path, "r") as file:
                lines = file.readlines()
                target_line = lines[8]
                matches = re.findall(r"\d+\.\d+|\d+", target_line)
                numbers = [float(match) if "." in match else int(match) for match in matches]
                line = [species, species_id, numbers[1], numbers[2], numbers[3], numbers[4], numbers[5]]
                data_list.append(line)

        # Create DataFrame for plotting
        df = pd.DataFrame(data_list, columns=["Species", "ID", "S", "D", "F", "M", "N"]).iloc[::-1].reset_index(drop=True)

        # Customize plot
        # ax.spines[["left", "right", "top"]].set_visible(False)
        for spine in ["left", "right", "top"]:
            ax.spines[spine].set_visible(False)

        ax.set_axisbelow(True)
        position = range(len(file_paths))

        # Plot bars
        ax.barh(position, df["S"], height=0.9, label="Complete and single-copy BUSCOs (S)", color=colors[0])
        ax.barh(position, df["D"], height=0.9, left=df["S"], label="Complete and duplicated BUSCOs (D)", color=colors[1])
        ax.barh(position, df["F"], height=0.9, left=df["S"] + df["D"], label="Fragmented BUSCOs (F)", color=colors[2])
        ax.barh(position, df["M"], height=0.9, left=df["S"] + df["D"] + df["F"], label="Missing BUSCOs (M)", color=colors[3])

        # Add legend
        ax.legend(ncol=legend_ncol, loc=legend_loc, handlelength=0.8, frameon=False)

        # Set ticks and labels
        # ax.yaxis.set_ticks_position('none')
        ax.set_yticks(range(len(df["Species"].tolist())))
        ax.set_yticklabels(["" for _ in range(len(df["Species"].tolist()))])
        ax.set_xlim(0, 100)
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{i}%" for i in xticks])

        # Invert bold_species_indices
        if bold_species_indices:
            bold_species_indices = [len(df) - 1 - idx for idx in bold_species_indices]

        # Annotate species, IDs, and BUSCO scores
        ax.text(-1, len(df.index), "Species", va="center", ha="right", fontweight="semibold")
        ax.text(1, len(df.index), "ID", va="center", ha="left", fontweight="semibold") if not df_input["ID"].eq("").all() else None
        for i, row in df.iterrows():
            fontweight = "bold" if bold_species_indices and i in bold_species_indices else "medium"
            ax.text(xticks[0] - 1, i, f'{row["Species"]}', va="center", ha="right", fontweight=fontweight, style="italic")
            ax.text(1, i, f'{row["ID"]}', va="center", ha="left", fontweight="bold", color="white")
            ax.text(
                df["S"].min() - 2,
                i,
                f"S: {row['S']}%   |   D: {row['D']}%   |   F: {row['F']}%   |   M: {row['M']}%",
                va="center",
                ha="right",
                fontweight="bold",
                color="white",
            )

    def draw_histogram_with_stats(self, ax, data, bins, show_legend=True, legend_loc="upper right", legend_ncol=1):
        """
        Draws a histogram with key statistics.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the histogram.

        data : array-like
            The dataset to plot as a histogram.

        bins : int or sequence of scalars or str
            Number of bins or bin edges for the histogram.

        show_legend : bool, optional
            Whether to display a legend. Default is True.

        legend_loc : str, optional
            Location of the legend on the plot. Default is "upper right".

        legend_ncol : int, optional
            Number of columns in the legend. Default is 1.

        Additional Details:
        -------------------
        - Adds vertical lines for:
            * Mean (red dashed line).
            * Median (blue dashed line).
            * 5th and 95th percentiles (purple dashed-dotted lines).
        - Shades regions within ±1, ±2, and ±3 standard deviations of the mean:
            * ±1σ in orange.
            * ±2σ in yellow.
            * ±3σ in green.
        """
        # Customize plot
        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)

        # Calculating statistics
        mean = np.mean(data)
        median = np.median(data)
        std_dev = np.std(data)
        percentile_5 = np.percentile(data, 5)
        percentile_95 = np.percentile(data, 95)

        # Mean line (red dashed)
        ax.axvline(mean, color="red", linestyle="--", label=f"Mean: {mean:.0f}")

        # Median line (blue dashed)
        ax.axvline(median, color="blue", linestyle="--", label=f"Median: {median:.0f}")

        # 5th and 95th percentiles lines
        ax.axvline(percentile_5, color="purple", linestyle="-.", label=f"5th percentile: {percentile_5:.0f}")
        ax.axvline(percentile_95, color="purple", linestyle="-.", label=f"95th percentile: {percentile_95:.0f}")

        # Shaded areas for mean ± std_dev * 1, 2, 3
        for i in range(1, 4):
            color = "orange" if i == 1 else "yellow" if i == 2 else "green"
            ax.axvspan(mean - i * std_dev, mean + i * std_dev, color=color, alpha=0.1)
            left_value = mean - i * std_dev
            right_value = mean + i * std_dev
            ax.axvline(left_value, color=color, linestyle="--", label=f"Mean - {i}σ: {left_value:.0f}")
            ax.axvline(right_value, color=color, linestyle="--", label=f"Mean + {i}σ: {right_value:.0f}")

        # Plotting histogram
        ax.hist(data, bins=bins, edgecolor="black")

        # Legend and labels
        if show_legend:
            ax.legend(ncol=legend_ncol, loc=legend_loc)

    def classify_roh(self, length):
        if length < 1_000_000:
            return "S"
        elif length >= 10_000_000:
            return "UL"
        else:
            return "L"

    def draw_roh_barplot(
        self,
        ax,
        data,
        genome_length,
        colors={"N": "#23b4e8", "S": "#008dbf", "L": "#fbbc04", "UL": "#ea4335"},
        xticks=[50, 60, 70, 80, 90, 100],
        xlim=(45, 100),
        sorting=True,
        groups=None,
        show_annotation=False,
        show_legend=True,
        legend_loc=(0.25, 0.97),
        legend_ncol=4,
    ):
        """
        Visualizes the distribution of ROHs (Runs of Homozygosity) across genome categories
        using a horizontal stacked bar plot.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the bar plot.

        data : list of str
            A list of file paths to BED files, where each file contains ROH data for a single sample.
            Each file should have tab-delimited columns: `scaffold`, `start`, `end`, `length`.

        genome_length : int or dict
            The total length of the genome, used to calculate percentages.
            If groups, dictionary with groups (keys) and genome sizes (values).

        colors : dict, optional
            A dictionary mapping ROH categories to their respective colors.
            Default is `{"N": "#23b4e8", "S": "#008dbf", "L": "#fbbc04", "UL": "#ea4335"}`.

        xticks : list of int, optional
            X-axis tick positions, representing percentages. Defaults to [0, 25, 50, 75, 100].

        xlim : tuple of float, optional
            Limits for the x-axis (log scale), representing the range of tract lengths to display. Defaults to (45, 100).

        sorting : bool, optional
            If True, bars are sorted by the percentage of Ultra Long ROHs (UL) in descending order.
            If False, bars follow the reverse order of `data`. Default is False.

        groups : dict, optional
            A dictionary where keys are group names and values are lists of sample names.
            This allows grouping samples on the y-axis with the group name followed by sample names.
            If None, no grouping is applied. Default is None.

        show_annotation : bool, optional
            Whether to display an annotation. Default is False.

        show_legend : bool, optional
            Whether to display a legend. Default is True.

        legend_loc : str, optional
            Location of the legend on the plot. Default is (0.25, 0.97).

        legend_ncol : int, optional
            Number of columns in the legend. Default is 4.

        Additional Details:
        --------------------
        - Each BED file represents a sample and contains ROH data with the following format:
            * `scaffold`: Chromosome or scaffold identifier.
            * `start`: Start position of the ROH.
            * `end`: End position of the ROH.
            * `length`: Length of the ROH (in base pairs).
        - Categories:
            * `"N"`: Non-ROHs (portion of the genome not covered by ROHs).
            * `"S"`: Short ROHs (<1,000,000 bp).
            * `"L"`: Long ROHs (1,000,000–10,000,000 bp).
            * `"UL"`: Ultra Long ROHs (≥10,000,000 bp).
        - The function calculates the percentage of the genome occupied by each ROH category for each sample.
        - The stacked bar plot ensures each bar represents 100% of the genome for a given sample.
        """
        # Customize plot
        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["left", "right", "top"]:
            ax.spines[spine].set_visible(False)

        # Load data for each file
        sample_data = {}
        sample_names = []
        for file_path in data:
            sample_name = os.path.basename(file_path).split(".")[0]
            sample_names.append(sample_name)
            df = pd.read_csv(file_path, sep="\t", header=None, names=["scaffold", "start", "end", "length"])
            df["length"] = pd.to_numeric(df["length"])
            df["classification"] = df["length"].apply(self.classify_roh)

            total_lengths = df.groupby("classification")["length"].sum().to_dict()

            for category in ["S", "L", "UL"]:
                if not groups:
                    total_lengths[category] = (total_lengths.get(category, 0) / genome_length) * 100
                else:
                    sample_group = next((key for key, values in groups.items() if sample_name in values), None)
                    total_lengths[category] = (total_lengths.get(category, 0) / genome_length[sample_group]) * 100

            # Add Non-ROHs
            total_roh_percentage = sum(total_lengths.values())
            total_lengths["N"] = max(0, 100 - total_roh_percentage)

            sample_data[sample_name] = total_lengths

        result_df = pd.DataFrame.from_dict(sample_data, orient="index").fillna(0)
        result_df = result_df[["N", "S", "L", "UL"]]

        # Sort the DataFrame if sorting=True
        if sorting:
            result_df = result_df.sort_values(by=["UL", "L", "S", "N"], ascending=[False, False, False, False])

        # If groups are not specified, all samples are considered as one group
        if not groups:
            groups = {"All": result_df.index.tolist()}

        # Create an ordered list of indices based on grouping
        grouped_samples = []
        for group, samples in groups.items():
            # Include only samples present in the DataFrame
            group_samples = [sample for sample in samples if sample in result_df.index]
            if sorting:
                # Sort within the group
                group_df = result_df.loc[group_samples]
                group_samples = group_df.sort_values(by=["UL", "L", "S", "N"], ascending=[False, False, False, False]).index.tolist()
            grouped_samples.extend(group_samples)

        # Reorder the DataFrame indices based on groups
        if groups and sorting:
            result_df = result_df.loc[grouped_samples]
        elif groups and not sorting:
            result_df = result_df.loc[sample_names]

        # Plot a cumulative bar plot (replace vertical bars with horizontal ones)
        category_labels = {"N": "Non-RoHs (N)", "S": "Short RoHs (S)", "L": "Long RoHs (L)", "UL": "Ultra Long RoHs (UL)"}
        bottom = None
        for category, color in colors.items():
            ax.barh(
                result_df.index,
                result_df[category],
                left=bottom,
                color=color,
                label=category_labels[category],
                height=0.9,
            )

            bottom = result_df[category] if bottom is None else bottom + result_df[category]

        # Configure initial area shading
        if xticks[0] != 0:
            ax.axvspan(xlim[0], xticks[0], color="white")

        # Configure axes and legend for the horizontal plot
        ax.yaxis.set_ticks_position("none")

        # Remove default Y-axis ticks
        ax.set_yticks([])
        ax.set_yticklabels([])

        # Add custom yticklabels
        for i, label in enumerate(result_df.index):
            ax.text(xticks[0] - 0.5, i, label, va="center", ha="right")

        ax.set_xlim(xlim)
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{i}%" for i in xticks])
        ax.spines["bottom"].set_bounds(xticks[0], xlim[1])
        if show_legend:
            ax.legend(loc=legend_loc, ncol=legend_ncol, handlelength=0.8, frameon=False)

        # Add vertical lines to indicate groups
        if groups and "All" not in groups:
            for group, samples in groups.items():
                group_samples = [sample for sample in result_df.index if sample in samples]
                if group_samples:
                    # Indices of the first and last samples in the group
                    start_idx = result_df.index.get_loc(group_samples[0])
                    end_idx = result_df.index.get_loc(group_samples[-1])

                    # Draw a vertical line
                    ax.vlines(
                        x=xlim[0],
                        ymin=start_idx - 0.25,
                        ymax=end_idx + 0.25,
                        colors="black",
                        linewidth=1,
                    )

                    # Add group name
                    y_pos = (start_idx + end_idx) / 2  # Average value for correct positioning
                    ax.text(xlim[0] - 1, y_pos, group, rotation=0, va="center", ha="right", fontstyle="italic")

        # Add annotation on each bar
        if show_annotation:
            for i, sample in enumerate(result_df.index):
                values = result_df.loc[sample, ["N", "S", "L", "UL"]]
                text = f'  N: {values["N"]:.1f}%   |   S: {values["S"]:.1f}%   |   L: {values["L"]:.1f}%   |   UL: {values["UL"]:.1f}%'
                print(f"{sample}\t{values['N']:.1f}%\t{values['S']:.1f}%\t{values['L']:.1f}%\t{values['UL']:.1f}%")
                ax.text(
                    xticks[0],
                    i,
                    text,
                    va="center",
                    ha="left",
                    fontsize=8,
                    color="white",
                    fontweight="bold",
                )

    def read_series(self, s):
        return pd.read_csv(s, header=None).squeeze("columns") if os.path.exists(s) else pd.Series(s.split(","))

    def rgb_tuple_to_hex(self, rgb_tuple):
        color_code = "#"
        for i in [0, 1, 2]:
            color_code += "{:02X}".format(int(255 * rgb_tuple[i]))
        return color_code

    def draw_variant_window_densities(
        self,
        ax,
        input_file,
        input_type="vcf",
        output_prefix=None,
        output_formats=[],
        title=False,
        title_fontsize=20,
        window_size=1000000,
        window_step=None,
        density_multiplier=1000,
        scaffold_white_list=pd.Series(dtype=str),
        scaffold_ordered_list=pd.Series(dtype=str),
        scaffold_length_file=[],
        scaffold_syn_file=None,
        syn_file_key_column=0,
        syn_file_value_column=1,
        figure_width=10,
        figure_height_per_scaffold=0.5,
        figure_header_height=0,
        colormap="jet",
        custom_color_list=None,
        density_thresholds=(0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5),
        density_thresholds_expression_type="left_open",
        skip_top_interval=False,
        skip_bottom_interval=False,
        test_colormaps=False,
        hide_track_label=True,
        subplots_adjust_left=None,
        subplots_adjust_top=None,
        subplots_adjust_right=None,
        subplots_adjust_bottom=None,
        only_count=False,
        x_tick_fontsize=None,
        stranded=False,
        rounded=True,
        middle_break=False,
        stranded_end=False,
        feature_name="SNPs",
        centromere_bed=None,
        highlight_file=None,
    ):
        # coverage=None,
        # scaffold_column_name="scaffold",
        # window_column_name="window",
        # coverage_column_name="median",
        # mean_coverage=None,
        # max_coverage_threshold=2.5,
        # min_coverage_threshold=0.5,
        # scaffold_black_list=pd.Series(dtype=str),
        # sort_scaffolds=False,
        # masking_gff_list=None,
        # masking_threshold=0.5,

        # if not output_prefix
        #     raise ValueError("Output prefix is required.")

        scaffold_white_list = self.read_series(scaffold_white_list)
        scaffold_ordered_list = self.read_series(scaffold_ordered_list)

        scaffold_ordered_list = scaffold_ordered_list[::-1]

        if isinstance(scaffold_ordered_list, list):
            if not scaffold_ordered_list:
                scaffold_ordered_list = scaffold_white_list
        else:
            if scaffold_ordered_list.empty:
                scaffold_ordered_list = scaffold_white_list

        if custom_color_list is not None:
            if len(custom_color_list) != len(density_thresholds):
                raise ValueError(
                    "ERROR!!! Custom color list is set, but the number of colors ({0}) in the list is not equal to the number of the thresholds (1)!".format(
                        len(custom_color_list), len(density_thresholds)
                    )
                )

        variants = CollectionVCF(input_file, parsing_mode="only_coordinates")

        chr_len_df = (
            pd.read_csv(scaffold_length_file, sep="\t", header=None, index_col=0) if scaffold_length_file else deepcopy(variants.scaffold_length)
        )
        chr_len_df.index = pd.Index(map(str, chr_len_df.index))
        chr_len_df.index.name = "scaffold"
        chr_len_df.columns = ["length"]

        chr_syn_dict = SynDict(filename=scaffold_syn_file, key_index=syn_file_key_column, value_index=syn_file_value_column)

        if centromere_bed:
            centromere_df = pd.read_csv(centromere_bed, usecols=(0, 1, 2), index_col=0, header=None, sep="\t", names=["scaffold_id", "start", "end"])
            centromere_df.rename(index=chr_syn_dict, inplace=True)
        else:
            centromere_df = None

        if input_type == "vcf":
            count_df = StatsVCF.count_variants_in_windows(
                variants,
                window_size,
                window_step,
                reference_scaffold_lengths=chr_len_df,
                ignore_scaffolds_shorter_than_window=True,
                # output_prefix=output_prefix,
                skip_empty_windows=False,
                expression=None,
                per_sample_output=False,
                scaffold_white_list=scaffold_white_list,
                scaffold_syn_dict=chr_syn_dict,
            )
            feature_df, track_df = StatsVCF.convert_variant_count_to_feature_df(count_df, window_size, window_step)
            # feature_df.to_csv("{}.features.counts".format(output_prefix), sep="\t", header=True, index=True)
            feature_df[feature_df.columns[-1]] = feature_df[feature_df.columns[-1]] * float(density_multiplier) / float(window_size)

            # feature_df.to_csv("{}.features.bed".format(output_prefix), sep="\t", header=True, index=True)

        elif input_type == "bedgraph":
            track_df = pd.read_csv(
                input,
                sep="\t",
                names=["scaffold", "start", "end", "value"],
                header=None,
                index_col=0,
                na_values=".",
                dtype={"scaffold": str, "start": int, "end": int, "value": float},
            )
            if scaffold_syn_file:
                track_df.rename(index=chr_syn_dict, inplace=True)
            track_df["value"] = track_df["value"].astype(float)

        if scaffold_syn_file:
            chr_len_df.rename(index=chr_syn_dict, inplace=True)

        # scale counts
        track_df[track_df.columns[-1]] = track_df[track_df.columns[-1]] * float(density_multiplier) / float(window_size)

        if track_df.index.nlevels > 1:
            # drop second level of index if it was added by groupby
            track_df = (
                track_df.groupby("scaffold").apply(lambda df: df[df["end"] <= chr_len_df.loc[df.index[0], "length"]]).reset_index(level=1, drop=True)
            )

        if not only_count:
            if custom_color_list is not None:
                cmap_list = ["custom_list"]
            else:
                cmap_list = Visualization.colormap_list if test_colormaps else [colormap]

            for colormap in cmap_list:
                if colormap == "custom_list":
                    colors = custom_color_list
                else:
                    cmap = plt.get_cmap(colormap, len(density_thresholds))
                    colors = [self.rgb_tuple_to_hex(cmap(i)[:3]) for i in range(0, len(density_thresholds))]

                color_expression = partial(
                    Visualization.color_threshold_expression,
                    thresholds=density_thresholds,
                    colors=colors,
                    background="white",
                    interval_type=density_thresholds_expression_type,
                    skip_top_interval=skip_top_interval,
                    skip_bottom_interval=skip_bottom_interval,
                )

                track_with_colors_df = Visualization.add_color_to_track_df(
                    track_df, color_expression, value_column_index=-1  # TODO fix it, add support for multiple tracks in the file
                )

                # track_with_colors_df.to_csv("{}.{}.track.bed".format(output_prefix, colormap), sep="\t", header=True, index=True)
                # print(feature_with_colors_df)
                # print(scaffold_ordered_list)
                Visualization.draw_features(
                    {"TR": track_with_colors_df},
                    chr_len_df,
                    scaffold_ordered_list,
                    output_prefix,
                    legend=Visualization.density_legend(
                        colors,
                        density_thresholds,
                        feature_name=feature_name,
                        interval_type=density_thresholds_expression_type,
                        skip_top_interval=skip_top_interval,
                        skip_bottom_interval=skip_bottom_interval,
                    ),
                    centromere_df=centromere_df,
                    highlight_df=highlight_file,
                    figure_width=figure_width,
                    figure_height_per_scaffold=figure_height_per_scaffold,
                    figure_header_height=figure_header_height,
                    dpi=300,
                    default_color="red",
                    title=title,
                    extensions=output_formats,
                    feature_start_column_id="start",
                    feature_end_column_id="end",
                    feature_color_column_id="color",
                    feature_length_column_id="length",
                    subplots_adjust_left=subplots_adjust_left,
                    subplots_adjust_bottom=subplots_adjust_bottom,
                    subplots_adjust_right=subplots_adjust_right,
                    subplots_adjust_top=subplots_adjust_top,
                    show_track_label=not hide_track_label,
                    show_trackgroup_label=True,
                    close_figure=False,
                    subplot_scale=False,
                    track_group_scale=False,
                    # track_group_distance=2,
                    # xmax_multiplier=1.3,
                    # ymax_multiplier=1.00,
                    stranded_tracks=stranded,
                    rounded_tracks=rounded,
                    middle_break=middle_break,
                    stranded_end_tracks=stranded_end,
                    xtick_fontsize=x_tick_fontsize,
                    subplot_title_fontsize=title_fontsize,
                    subplot_title_fontweight="bold",
                    axes=ax,
                )

    def draw_features(
        self,
        ax,
        input_file,
        input_type="str",
        header=None,
        legend=None,
        output_prefix=None,
        output_formats=[],
        title=False,
        start_column_name=None,
        end_column_name=None,
        color_column_name=None,
        default_color="tab:blue",
        scaffold_white_list=pd.Series(dtype=str),
        scaffold_ordered_list=pd.Series(dtype=str),
        scaffold_length_file=None,
        scaffold_syn_file=None,
        syn_file_key_column=0,
        syn_file_value_column=1,
        colormap="jet",
        hide_track_label=True,
        x_tick_type="nucleotide",
        feature_shape="rectangle",
        subplots_adjust_left=None,
        subplots_adjust_top=None,
        subplots_adjust_right=None,
        subplots_adjust_bottom=None,
        figure_width=10,
        figure_height_per_scaffold=0.5,
        figure_header_height=0,
        verbose=False,
        subplot_scale=False,
        track_group_scale=False,
        x_tick_fontsize=None,
        stranded=False,
        rounded=True,
        stranded_end=False,
        fill_empty_tracks=False,
        empty_color="lightgrey",
        centromere_bed=None,
        highlight_file=None,
        title_fontsize=20,
    ):
        # scaffold_column_name=None,
        # scaffold_black_list=pd.Series(dtype=str),
        # sort_scaffolds=False,
        # figure_height_per_scaffold=0.5,
        # print(scaffold_ordered_list)
        scaffold_white_list = self.read_series(scaffold_white_list)
        scaffold_ordered_list = self.read_series(scaffold_ordered_list)
        scaffold_ordered_list = scaffold_ordered_list[::-1]

        chr_syn_dict = SynDict(filename=scaffold_syn_file, key_index=syn_file_key_column, value_index=syn_file_value_column)

        # print(scaffold_ordered_list)
        if isinstance(scaffold_ordered_list, list):
            if not scaffold_ordered_list:
                scaffold_ordered_list = deepcopy(scaffold_white_list)
                scaffold_ordered_list.replace(chr_syn_dict, inplace=True)
        else:
            # print("AAAA")
            if scaffold_ordered_list.empty:
                scaffold_ordered_list = deepcopy(scaffold_white_list)
                scaffold_ordered_list.replace(chr_syn_dict, inplace=True)

        if centromere_bed:
            centromere_df = pd.read_csv(centromere_bed, usecols=(0, 1, 2), index_col=0, header=None, sep="\t", names=["scaffold_id", "start", "end"])
            centromere_df.rename(index=chr_syn_dict, inplace=True)
        else:
            centromere_df = None
        try:
            if input_type == "str":
                feature_df = CollectionSTR(
                    in_file=input_file, records=None, format="filtered_str", parsing_mode="all", black_list=(), white_list=(), syn_dict=chr_syn_dict
                )

                feature_df.records.set_index("scaffold_id", inplace=True)

                feature_start_column_id = start_column_name if start_column_name else "start"
                feature_end_column_id = end_column_name if end_column_name else "end"

            elif input_type == "gff":
                feature_df = CollectionGFF(in_file=input_file, parsing_mode="only_coordinates")

                feature_start_column_id = start_column_name if start_column_name else "start"
                feature_end_column_id = end_column_name if end_column_name else "end"

            elif input_type == "bed":
                feature_df = CollectionBED(in_file=input_file, parsing_mode="coordinates_only", format="bed")

                feature_start_column_id = "start"
                feature_end_column_id = "end"

            elif input_type == "bed_table":
                feature_df = CollectionBED(in_file=input_file, parsing_mode="complete", format="bed", header_in_file=False)
                feature_start_column_id = "start"
                feature_end_column_id = "end"

            elif input_type == "bedgraph":
                feature_df = CollectionBED(in_file=input_file, parsing_mode="all", format="bed")
                feature_df.records.columns = pd.Index(["start", "end", "value"])
                feature_df.records.index.name = "scaffold"
                feature_start_column_id = "start"
                feature_end_column_id = "end"
            elif input_type == "bed_track":
                feature_df = CollectionBED(
                    in_file=input_file,
                    parsing_mode="all",
                    format="bed_track",
                )
                print(feature_df.records)
                feature_df.records.columns = pd.Index(["start", "end", "value", "color"])
                feature_df.records.index.name = "scaffold"
                feature_start_column_id = "start"
                feature_end_column_id = "end"
                color_column_name = "color"
            elif input_type == "tab6":
                feature_df = CollectionBLAST(in_file=input_file, parsing_mode="complete")
                feature_df.records.reset_index(level="query_id", inplace=True)
                feature_start_column_id = start_column_name if start_column_name else "target_start"
                feature_end_column_id = end_column_name if end_column_name else "target_end"
            elif input_type == "tab6_colored":
                feature_df = CollectionBLAST(in_file=input_file, parsing_mode="complete", format="tab6_colored", header=header)
                feature_df.records.reset_index(level="query_id", inplace=True)
                feature_start_column_id = start_column_name if start_column_name else "target_start"
                feature_end_column_id = end_column_name if end_column_name else "target_end"
            else:
                raise ValueError("ERROR!!! Unrecognized input type ({}). ".format(input_type))
        except pd.errors.EmptyDataError:
            print(
                "Empty input file. Silent exit."
            )  # try-except added to handle case when input file is empty without raising exception. For use in snakemake
            exit(0)

        legend_df = pd.read_csv(legend, header=None, index_col=0, sep="\t") if legend else None

        # print(scaffold_white_list)
        # print(feature_df.records)
        scaffold_to_keep = StatsVCF.get_filtered_entry_list(
            feature_df.records.index.get_level_values(level=0).unique().to_list(), entry_white_list=scaffold_white_list
        )
        # print(scaffold_to_keep)
        # print(scaffold_to_keep)
        # remove redundant scaffolds
        # print(scaffold_white_list)
        # print(feature_df.records)
        # print(scaffold_to_keep)
        feature_df.records = feature_df.records[feature_df.records.index.isin(scaffold_to_keep)]
        # print("BBBBBBbb")
        # print(scaffold_white_list)
        # print(scaffold_ordered_list)
        if not scaffold_white_list.empty:
            scaffold_ordered_list = scaffold_ordered_list[scaffold_ordered_list.isin(pd.Series(scaffold_white_list).replace(chr_syn_dict))]
        # print("CCCC")
        # print(scaffold_ordered_list)
        # print(scaffold_to_keep)
        # print(pd.Series(scaffold_to_keep).replace(chr_syn_dict))
        chr_len_df = pd.read_csv(scaffold_length_file, sep="\t", header=None, names=("scaffold", "length"), index_col=0)
        chr_len_df.index = pd.Index(list(map(str, chr_len_df.index)))

        # print(chr_len_df)

        # print(feature_df.records)
        if scaffold_syn_file:
            chr_len_df.rename(index=chr_syn_dict, inplace=True)
            feature_df.records.rename(index=chr_syn_dict, inplace=True)
        if verbose:
            print(chr_syn_dict)
            print(feature_df.records)

        # print(chr_len_df)
        # print({"features": feature_df})
        # print(feature_df.records.columns)
        # print(feature_df.records)
        # print(chr_len_df)
        # print(scaffold_ordered_list)
        Visualization.draw_features(
            {"features": feature_df},
            chr_len_df,
            scaffold_ordered_list,
            output_prefix,
            legend=Visualization.feature_legend(legend_df, colormap=colormap),
            # legend_df=legend_df,
            centromere_df=centromere_df,
            highlight_df=highlight_file,
            figure_width=figure_width,
            figure_height_per_scaffold=figure_height_per_scaffold,
            figure_header_height=figure_header_height,
            dpi=300,
            # colormap=None, thresholds=None, colors=None, background=None,
            default_color=default_color,
            title=title,
            extensions=output_formats,
            feature_shape=feature_shape,
            feature_start_column_id=feature_start_column_id,
            feature_end_column_id=feature_end_column_id,
            feature_color_column_id=color_column_name,
            feature_length_column_id="length",
            subplots_adjust_left=subplots_adjust_left,
            subplots_adjust_bottom=subplots_adjust_bottom,
            subplots_adjust_right=subplots_adjust_right,
            subplots_adjust_top=subplots_adjust_top,
            show_track_label=not hide_track_label,
            show_trackgroup_label=True,
            subplot_scale=subplot_scale,
            track_group_scale=track_group_scale,
            stranded_tracks=stranded,
            rounded_tracks=rounded,
            stranded_end_tracks=stranded_end,
            fill_empty_tracks=fill_empty_tracks,
            empty_color=empty_color,
            xtick_fontsize=x_tick_fontsize,
            subplot_title_fontsize=title_fontsize,
            subplot_title_fontweight="bold",
            x_tick_type=x_tick_type,
            # xmax_multiplier=2,
            axes=ax,
        )

    def draw_psmc_plot(
        self,
        ax,
        diploid_data,
        round_data=None,
        n=100,
        colorlist=None,
        xlim=(50000, 23000000),
        xticks=[
            50000,
            70000,
            100000,
            150000,
            200000,
            300000,
            400000,
            500000,
            600000,
            800000,
            1000000,
            1500000,
            2000000,
            3000000,
            5000000,
            7000000,
            10000000,
            15000000,
            20000000,
        ],
        ylim=(0, 200),
        scale=10000,
        mu=2.2e-9,
        g=5,
        show_legend=True,
        legend_loc="upper right",
        legend_ncol=1,
    ):
        """
        Draw a PSMC plot using diploid and round data.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes object where the plot will be drawn.

        diploid_data : list of str
            A list of file paths containing the diploid PSMC data (each file with columns: 'x', 'y', 'z', 'w', 'h').
            The 'x' column represents time (in years ago), and 'y' represents the effective population size.

        round_data : list of str, optional
            A list of file paths containing the round PSMC data (default is None). Each file is expected to follow the
            same format as the diploid data. The round data will be visualized with a reduced alpha transparency.

        n : int, optional, default: 100
            Number of round data.

        colorlist : list of str or None, optional
            A list of colors for the plot. If None, colors are generated automatically. If a string is provided, it is
            interpreted as a seaborn color palette.

        xlim : tuple of (int, int), optional, default: (50000, 23000000)
            The x-axis limits (years ago). The x-axis is logarithmic.

        xticks : list of int, optional, default: [50000, 70000, 100000, ..., 20000000]
            The ticks to display on the x-axis (years ago), with corresponding labels adjusted according to the `scale`.

        ylim : tuple of (int, int), optional, default: (0, 200)
            The y-axis limits, representing the effective population size.

        scale : int, optional, default: 10000
            The scaling factor for the x-axis labels.

        mu : float, optional, default: 2.2e-9
            The mutation rate used for PSMC analysis (used in axis label formatting).

        g : int, optional, default: 5
            The generation time used for PSMC analysis (used in axis label formatting).

        show_legend : bool, optional, default: True
            Whether to display the legend on the plot.

        legend_loc : str, optional, default: "upper right"
            The location of the legend on the plot.

        legend_ncol : int, optional, default: 1
            The number of columns in the legend.
        """
        if colorlist is None:
            colorlist = distinctipy.get_colors(len(diploid_data))
        else:
            if type(colorlist) == str:
                colorlist = sns.color_palette(colorlist, len(diploid_data))

        for diploid, color in zip(diploid_data, colorlist):
            sample_name = diploid.split("/")[-1].split(".")[0]
            data = pd.read_csv(diploid, names=["x", "y", "z", "w", "h"], sep="\t")
            data = data[data["x"] >= 0]
            ax.step(data["x"], data["y"], where="post", color=color, linewidth=2, label=sample_name)

        if round_data is not None:
            for round, color in zip(round_data, colorlist):
                data = pd.read_csv(round, names=["x", "y", "z", "w", "h"], sep="\t")
                data = data[data["x"] >= 0]
                for i in range(0, len(data), len(data) // n):
                    ax.step(
                        data["x"].iloc[i : i + len(data) // n],
                        data["y"].iloc[i : i + len(data) // n],
                        where="post",
                        color=color,
                        linewidth=1,
                        alpha=0.1,
                    )

        ax.set_xscale("log")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ax.set_xticks(xticks)
        ax.set_xticklabels([int(x / scale) for x in xticks])

        scale_exp = int(np.log10(scale))  # Calculate exponent for the scale
        ax.set_xlabel(rf"Years Ago, $10^{scale_exp}$ ($\mu={mu:.1e}$, g={g})")
        ax.set_ylabel(rf"Effective population size, $10^{scale_exp}$")

        if title is not None:
            ax.set_title(r"$\mathit{" + title.replace(" ", r"\,") + "}$")

        if show_legend:
            ax.legend(loc=legend_loc, ncol=legend_ncol, frameon=True)
