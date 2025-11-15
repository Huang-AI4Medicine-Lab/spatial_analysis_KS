from __future__ import annotations

from pathlib import Path

from copy import copy
from types import MappingProxyType
from typing import Any, Mapping

import matplotlib as mpl
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns
import squidpy as sq
from anndata import AnnData
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import PathPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scanpy.pl._dotplot import DotPlot
from scipy.cluster import hierarchy as sch
from squidpy._constants._pkg_constants import Key

from squidpy.gr._utils import _assert_categorical_obs
from squidpy.pl._color_utils import Palette_t, _get_palette, _maybe_set_colors

try:
    from matplotlib import colormaps as cm
except ImportError:
    from matplotlib import cm

plt.rcParams['svg.fonttype'] = 'none'


def _reorder(values, order, axis=1):
    if axis == 0:
        values = values.iloc[order, :]
    elif axis == 1:
        values = values.iloc[:, order]
    else:
        raise ValueError("The axis parameter accepts only values 0 and 1.")
    return values


def _clip(values, min_threshold=None, max_threshold=None, new_min=None, new_max=None, new_middle=None):
    values_clipped = values.copy()
    if new_middle is not None:
        values_clipped[:] = new_middle
    if min_threshold is not None:
        values_clipped[values < min_threshold] = new_min if new_min is not None else min_threshold
    if max_threshold is not None:
        values_clipped[values > max_threshold] = new_max if new_max is not None else max_threshold
    return values_clipped


class MyDotPlot(DotPlot):
    """Modified version :class:`scanpy.pl.DotPlot`."""

    def _plot_size_legend(self, size_legend_ax: Axes):
        size_range = np.linspace(self.dot_min, self.size_threshold, 3)
        if self.dot_min == 0:
            size_range[0] += self.size_threshold / 10

        size = (size_range / (self.size_threshold - self.dot_min)) ** self.size_exponent
        size = size * (self.largest_dot - self.smallest_dot) + self.smallest_dot
        # plot size bar
        size_legend_ax.scatter(
            np.arange(len(size)) + 0.5,
            np.repeat(0, len(size)),
            s=size,
            c=["gray" if s < self.color_threshold else (0.705673158, 0.01555616, 0.150232812, 1.0) for s in size_range],
            edgecolor="black",
            linewidth=self.dot_edge_lw,
            zorder=100,
        )
        size_legend_ax.set_xticks(np.arange(len(size)) + 0.5)
        labels = [f"{np.round((x * self.max_value), decimals=2)}" for x in size_range]
        labels[-1] = f">{labels[-1]}"
        size_legend_ax.set_xticklabels(labels, fontsize="x-small")

        # remove y ticks and labels
        size_legend_ax.tick_params(axis="y", left=False, labelleft=False, labelright=False)

        # remove surrounding lines
        size_legend_ax.spines["right"].set_visible(False)
        size_legend_ax.spines["top"].set_visible(False)
        size_legend_ax.spines["left"].set_visible(False)
        size_legend_ax.spines["bottom"].set_visible(False)
        size_legend_ax.grid(False)

        ymax = size_legend_ax.get_ylim()[1]
        size_legend_ax.set_ylim(-1.05 - self.largest_dot * 0.003, 4)
        size_legend_ax.set_title(self.size_title, y=ymax + 0.45, size="small")

        xmin, xmax = size_legend_ax.get_xlim()
        size_legend_ax.set_xlim(xmin - 0.15, xmax + 0.5)

    # ToDo: need to find a way to get get_axes()['mainplot_ax'] without showing the plot
    def _rotate_xlabels(self, ax):
        ax = ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="center", minor=False)

def _dotplot(
    adata: AnnData,
    x_key: str,
    y_key: str,
    values: np.ndarray,
    abs_values: bool = False,
    size_threshold: tuple[float, float] | tuple[None, None] = (None, None),
    color_threshold: tuple[float, float] = (-1, 1),
    figsize: tuple[float, float] | None = None,
    cmap: str | Palette_t = "bwr",
    size_title: str = "log2 FC",
    dot_scale: float = 1,
    cluster_y: bool = True,
    **kwargs,
):
    values_color = _clip(
        values, min_threshold=color_threshold[0], max_threshold=color_threshold[1], new_min=-1, new_max=1
    )
    values_color[(values < 0) & (values > color_threshold[0])] = 0  # -0.3
    values_color[(values > 0) & (values < color_threshold[1])] = 0  # 0.3

    if cluster_y is True:
        order = sp.cluster.hierarchy.dendrogram(
            sp.cluster.hierarchy.linkage(values_color.T, method="complete"), no_plot=True
        )["leaves"]
        values = _reorder(values, order, axis=1)
        values_color = _reorder(values_color, order, axis=1)

    one_hot_encoded = pd.get_dummies(adata.obs[y_key])

    adata_obs = AnnData(one_hot_encoded, dtype=np.uint8, obs=adata.obs)

    values_size = _clip(values, size_threshold[0], size_threshold[1])
    values_size = pd.DataFrame(
        (mcolors.TwoSlopeNorm(vcenter=0, vmin=size_threshold[0], vmax=size_threshold[1])(values_size) - 0.5) * 2,
        columns=values_size.columns,
        index=values_size.index,
    )

    if abs_values:
        print("Warning: label for depletion/enrichment to be implemented.")
        values_size = np.abs(values_size)

    if figsize is None:
        figsize = (10, 10 * values.shape[1] / values.shape[0])

    dp = MyDotPlot(
        adata_obs,
        adata_obs.var_names,
        groupby=x_key,
        dot_color_df=values_color,
        dot_size_df=values_size,
        figsize=figsize,
        **kwargs,
    )
    dp.max_value = np.max(values_size.values)
    dp.color_threshold = color_threshold[1]
    dp.size_threshold = size_threshold[1]

    dp.swap_axes()
    dp = dp.style(
        cmap=cmap,
        largest_dot=dp.largest_dot * dot_scale,
        dot_edge_lw=DotPlot.DEFAULT_DOT_EDGELW,
    )
    dp = dp.legend(show_colorbar=False, size_title=size_title)
    return dp

def _proportion(adata, id_key, val_key, normalize=True):
    df = pd.pivot(adata.obs[[id_key, val_key]].value_counts().reset_index(), index=id_key, columns=val_key)
    df[df.isna()] = 0
    df.columns = df.columns.droplevel(0)
    if normalize:
        return df.div(df.sum(axis=1), axis=0)
    else:
        return df

def proportion(
    adata,
    group_key: str,
    label_key: str,
    groups: list | None = None,
    labels: list | None = None,
    rotation_xlabel: int = 45,
    ncols: int = 1,
    normalize: bool = True,
    palette: Palette_t = None,
    figsize: tuple[float, float] | None = None,
    dpi: int | None = None,
    legend: bool = True,
    right_align_labels: bool = False,
    remove_box: bool = True,
    new_xlabels: list | None = None,
    sort_by_label: int | None = None,
    ascending: bool = True,
    remove_xlabels: bool = False,
    save: str | Path | None = None,
    return_group_order: bool = False,
    xlabel_fontsize: int | None = None,  # Font size for x-axis labels
    ylabel_fontsize: int | None = None,  # Font size for y-axis labels
    title_fontsize: int | None = None,   # Font size for the title
    legend_fontsize: int | None = None,  # Font size for the legend
    xlabel_colors: dict | None = None,   # New parameter for coloring x-axis labels
    **kwargs,
) -> list | None:
    """
    Plot the proportion of labels in groups, with options to sort by specific labels, customize appearance, return group order, save the plot, and adjust font sizes.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix in which the observations are the groups and the variables are the labels.
    group_key : str
        Key in `adata.obs` that specifies the groups to plot on the x-axis (e.g., samples, clusters, conditions).
    label_key : str
        Key in `adata.obs` that specifies the labels to plot as stacked bars on the y-axis (e.g., cell types, clusters).
    groups : list, optional (default: None)
        List of specific groups to include in the plot. If None, all groups will be plotted.
    labels : list, optional (default: None)
        List of specific labels to include in the plot. If None, all labels will be plotted.
    rotation_xlabel : int, optional (default: 45)
        Rotation angle for the x-axis labels (in degrees). Useful for long group names.
    ncols : int, optional (default: 1)
        Number of columns for the legend. If you have many labels, increasing `ncols` can make the legend more compact.
    normalize : bool, optional (default: True)
        If True, the proportions (relative frequencies) of labels within each group are shown. If False, raw counts are shown.
    palette : Palette_t, optional (default: None)
        Categorical colormap for the labels. If None, colors from `adata.uns` corresponding to the `label_key` will be used.
    figsize : tuple[float, float], optional (default: None)
        Size of the figure (width, height) in inches.
    dpi : int, optional (default: None)
        Dots per inch (DPI) for the plot, controlling the resolution.
    legend : bool, optional (default: True)
        Whether to show the legend. If False, the legend is hidden.
    right_align_labels : bool, optional (default: False)
        If True, the x-axis labels will be right-aligned.
    remove_box : bool, optional (default: True)
        If True, removes the plot box (spines) from the figure for a cleaner look.
    new_xlabels : list, optional (default: None)
        Custom labels to use for the x-axis. If None, the default group names will be used.
    sort_by_label : int, optional (default: None)
        Column (label) to sort the groups by. Groups will be sorted by the counts or proportions of this label.
    ascending : bool, optional (default: True)
        Sorting order for `sort_by_label`. If True, groups will be sorted in ascending order.
    remove_xlabels : bool, optional (default: False)
        If True, the x-axis labels are removed from the plot.
    save : str or Path, optional (default: None)
        Path where the figure will be saved. If None, the plot is not saved.
    return_group_order : bool, optional (default: False)
        If True, the function returns the order of the groups in the final plot.
    xlabel_fontsize : int, optional (default: None)
        Font size for x-axis labels. If None, the default size is used.
    ylabel_fontsize : int, optional (default: None)
        Font size for y-axis labels. If None, the default size is used.
    title_fontsize : int, optional (default: None)
        Font size for the plot title. If None, the default size is used.
    legend_fontsize : int, optional (default: None)
        Font size for the legend. If None, the default size is used.
    xlabel_colors : dict, optional (default: None)
        A dictionary mapping group names (x-axis labels) to colors for the labels. If None, the default color is used.
    **kwargs
        Additional keyword arguments passed to `pandas.DataFrame.plot.bar` for customizing the bar plot.
    
    Returns
    -------
    list | None
        If `return_group_order=True`, returns the order of the groups in the final plot. Otherwise, returns None.
    """
    _assert_categorical_obs(adata, key=group_key)
    _assert_categorical_obs(adata, key=label_key)
    _maybe_set_colors(source=adata, target=adata, key=label_key, palette=palette)

    clusters = adata.obs[label_key].cat.categories
    palette = _get_palette(adata, cluster_key=label_key, categories=clusters)

    df = _proportion(adata=adata, id_key=group_key, val_key=label_key, normalize=normalize)

    # Sort columns based on the count of a specific label (column)
    if sort_by_label is not None:
        if sort_by_label in df.columns:
            df = df[df.columns].sort_values(by=sort_by_label, ascending=ascending)
        else:
            raise ValueError(f"Label '{sort_by_label}' not found in the DataFrame columns.")

    if groups is not None:
        df = df.loc[groups, :]

    if labels is not None:
        df = df.loc[:, labels]

    plt.figure(dpi=dpi)
    ax = df.plot.bar(stacked=True, figsize=figsize, color=palette, rot=rotation_xlabel, ax=plt.gca(), **kwargs)
    ax.grid(False)

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles[::-1], labels[::-1], loc="center left", ncol=ncols, bbox_to_anchor=(1.0, 0.5))
        if legend_fontsize is not None:
            for text in lgd.get_texts():
                text.set_fontsize(legend_fontsize)
    else:
        ax.legend().set_visible(False)

    if new_xlabels:
        ax.set_xticklabels(new_xlabels)

    if remove_xlabels:
        ax.set_xticklabels([])  # Remove x-axis labels

    if right_align_labels:
        ax.set_xticklabels(ax.get_xticklabels(), ha='right')

    # Apply font sizes if specified
    if xlabel_fontsize is not None:
        ax.set_xlabel(ax.get_xlabel(), fontsize=xlabel_fontsize)
        ax.tick_params(axis='x', labelsize=xlabel_fontsize)
    if ylabel_fontsize is not None:
        ax.set_ylabel(ax.get_ylabel(), fontsize=ylabel_fontsize)
        ax.tick_params(axis='y', labelsize=ylabel_fontsize)
    if title_fontsize is not None:
        ax.set_title(ax.get_title(), fontsize=title_fontsize)

    # Apply colors to x-axis labels based on `xlabel_colors`
    if xlabel_colors is not None:
        xticklabels = ax.get_xticklabels()
        for label in xticklabels:
            if label.get_text() in xlabel_colors:
                label.set_color(xlabel_colors[label.get_text()])

    if remove_box:
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

    if save:
        if legend:
            plt.savefig(save, bbox_extra_artists=(lgd, lgd), bbox_inches="tight")
        else:
            plt.savefig(save, bbox_inches="tight")

    if return_group_order:
        return df.index.tolist()  # Return the order of groups in the plot

    return None


def enrichment(
    adata: AnnData,
    group_key: str,
    label_key: str,
    size_threshold: float | None = None,
    color_threshold: float = 1,
    legend_title: str | None = None,
    dot_scale: float = 1,
    cluster_labels: bool = True,
    groups: list | None = None,
    labels: list | None = None,
    palette: Palette_t | matplotlib.colors.ListedColormap | None = None,
    figsize: tuple[float, float] | None = None,
    save: str | Path | None = None,
    **kwargs,
):
    """
    Plot a dotplot of the enrichment of `y_key` in `x_key`.

    This functions is based on a modified version of :func:`scanpy.pl.dotplot`.

    Parameters
    ----------
    %(adata)s
    group_key
        Key in :attr:`anndata.AnnData.obs` where groups are stored.
    label_key
        Key in :attr:`anndata.AnnData.obs` where labels are stored.
    size_threshold
        Threshold for the size of the dots. Enrichments with value above this threshold will have all the same size.
    color_threshold
        Threshold to mark enrichments as significant.
    legend_title
        Title for the size legend.
    dot_scale
        Scale of the dots.
    cluster_groups
        If `True`, display labels ordered according to hierarchical clustering.
    groups
        The groups for which to show the enrichment.
    labels
        The labels for which to show the enrichment.
    palette
        Colormap for the enrichment values.
    %(plotting)s
    kwargs
        Keyword arguments for :func:`matplotlib.pyplot.scatter`.
    """
    if f"{group_key}_{label_key}_enrichment" not in adata.uns:
        raise ValueError("Run cellcharter.gr.enrichment first.")

    if palette is None:
        palette = matplotlib.colors.LinearSegmentedColormap.from_list(
            "", [cm.get_cmap("coolwarm")(0), matplotlib.colors.to_rgb("darkgrey"), cm.get_cmap("coolwarm")(255)]
        )

    enrichment = adata.uns[f"{group_key}_{label_key}_enrichment"]["enrichment"]

    if labels is not None:
        enrichment = enrichment.loc[:, labels]

    if groups is not None:
        enrichment = enrichment.loc[groups]

    size_threshold = np.max(enrichment.values) if size_threshold is None else size_threshold

    dp = _dotplot(
        adata if labels is None else adata[adata.obs[label_key].isin(labels)],
        x_key=group_key,
        y_key=label_key,
        values=enrichment,
        abs_values=False,
        size_threshold=(-1, size_threshold),
        color_threshold=(0, color_threshold),
        figsize=figsize,
        cmap=palette,
        size_title=legend_title,
        dot_scale=dot_scale,
        cluster_y=cluster_labels,
        **kwargs,
    )
    if save:
        dp.savefig(save, bbox_inches="tight")
    else:
        dp.show()