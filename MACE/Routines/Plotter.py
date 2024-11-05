#!/usr/bin/env python
__author__ = "tomarovsky"

import re
from pathlib import Path

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
            file_path = Path(count)
            file_name = file_path.stem
            id = file_name.split(".")[0]
            df = pd.read_csv(count, sep="\t")

            if removed_chrX:
                # If removed_chrX is a list, remove the specified chromosomes; otherwise, remove a single string
                if isinstance(removed_chrX, list):
                    for chrX in removed_chrX:
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

        ax.spines[["right", "top"]].set_visible(False)
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

        ax.spines[["right", "top"]].set_visible(False)
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
        font_size=None,
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

        font_size : int or float, optional
            Font size for the plot text. If not provided, the default font size is used.

        Notes:
        ------
        - The function reads ROH data from each file, sorts it by tract length, and calculates the cumulative sum of
          tract lengths normalized by the genome length.
        - The plot uses a logarithmic x-axis to display a wide range of homozygous tract lengths.
        - A step plot is created for each sample, and a legend is added to differentiate between samples.
        """
        ax.spines[["right", "top"]].set_visible(False)
        ax.set_axisbelow(True)

        if font_size is not None:
            plt.rcParams.update({"font.size": font_size})

        if colors is None:
            colors = distinctipy.get_colors(len(data))
        else:
            if type(colors) == str:
                colors = sns.color_palette(colors, len(data))

        for count, f in enumerate(data):
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
        """
        Draws a horizontal bar plot to visualize BUSCO (Benchmarking Universal Single-Copy Orthologs) results
        for multiple species.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the plot.

        data : list of str
            A list of file paths to BUSCO short summary text files (short_summary_{species}.txt).

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

        Notes:
        ------
        - The function reads each BUSCO summary file, extracts the BUSCO scores for complete single-copy,
          complete duplicated, fragmented, and missing categories, and organizes them into a DataFrame.
        - A horizontal stacked bar plot is created for each species, with bar segments colored according
          to the BUSCO category.
        - Species names and BUSCO scores are annotated alongside the bars.
        - The plot's aesthetics are customized by hiding the top and right spines, setting axis labels,
          and positioning a legend.
        """
        data_list = []
        for file_name in data:
            with open(file_name, "r") as file:
                lines = file.readlines()
                target_line = lines[8]
                matches = re.findall(r"\d+\.\d+|\d+", target_line)
                numbers = [float(match) if "." in match else int(match) for match in matches]
                species_info = file_name.split("/")[-1][14:-4]
                species, species_id = species_info.rsplit("_", 1)
                line = [species, species_id, numbers[1], numbers[2], numbers[3], numbers[4], numbers[5]]
                data_list.append(line)

        # Create DataFrame
        df = pd.DataFrame(data_list, columns=["Species", "ID", "S", "D", "F", "M", "N"])

        # Customize plot
        ax.spines[["left", "right", "top"]].set_visible(False)
        ax.set_axisbelow(True)
        position = range(len(data))

        # Plot bars
        ax.barh(position, df["S"], height=0.8, label="Complete and single-copy BUSCOs (S)", color=colors[0])
        ax.barh(position, df["D"], height=0.8, left=df["S"], label="Complete and duplicated BUSCOs (D)", color=colors[1])
        ax.barh(position, df["F"], height=0.8, left=df["S"] + df["D"], label="Fragmented BUSCOs (F)", color=colors[2])
        ax.barh(position, df["M"], height=0.8, left=df["S"] + df["D"] + df["F"], label="Missing BUSCOs (M)", color=colors[3])

        # Add legend
        ax.legend(ncol=legend_ncol, loc=legend_loc, handlelength=0.8, frameon=False)

        # Set ticks and labels
        ax.set_yticks(range(len(df["Species"].tolist())))
        ax.set_yticklabels(["" for _ in range(len(df["Species"].tolist()))])
        # ax.set_xlim(0, 100)
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{i}%" for i in xticks])

        # Annotate species and BUSCO scores
        ax.text(-1, len(data), "Вид", va="center", ha="right", fontweight="semibold")
        ax.text(1, len(data), "ID", va="center", ha="left", fontweight="semibold")
        for i, row in df.iterrows():
            ax.text(xticks[0] - 1, i, f'{row["Species"]}', va="center", ha="right", fontweight="medium", style="italic")
            ax.text(-1, i, f'{row["Species"]}', va="center", ha="right", fontweight="medium", style="italic")
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
