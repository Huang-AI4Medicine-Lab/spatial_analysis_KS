from seaborn.matrix import _HeatMapper
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from seaborn.utils import despine, axis_ticklabels_overlap, relative_luminance, to_utf8

class _ScatterMapper(_HeatMapper):
    """
    Draw a scattermap plot, similar to heatmap plot, but use scatter dots instead of heatmap
    """

    def __init__(self, data,
                 marker, marker_size,
                 vmin, vmax, cmap, center, robust, cbar, cbar_kws,
                 xticklabels=True, yticklabels=True, mask=None,
                 size_legend=True, size_legend_kws=None,
                 size_title='-log10(FDR)',
                 xticklabel_rotation=0,
                 xticklabel_ha='right'):

        super(_ScatterMapper, self).__init__(
            data, vmin, vmax, cmap, center, robust, cbar=cbar, cbar_kws=cbar_kws,
            xticklabels=xticklabels, yticklabels=yticklabels, mask=mask,
            # Don't support annotation
            annot=False, fmt=None, annot_kws=None,
        )

        self.marker = marker
        self.size_title = size_title
    
        self.cbar_title = self.cbar_kws.pop('title', None)
        self.size_legend = size_legend
        self.size_legend_kws = size_legend_kws or {}
        
        if isinstance(marker_size, float) or isinstance(marker_size, int):
            self.marker_size = marker_size
            self.size_legend = False  # Disable size legend if uniform size
        elif isinstance(marker_size, pd.DataFrame):
            self.marker_size = marker_size.loc[self.data.index, self.data.columns].values
        else:
            self.marker_size = marker_size
        
        # Ensure vmin and vmax are properly set
        if vmin is None or vmax is None:
            if robust:
                vmin, vmax = np.percentile(self.data, [2, 98])
            else:
                vmin, vmax = self.data.min().min(), self.data.max().max()
        
        # Convert to float to avoid pandas Series comparison issues
        vmin = float(vmin)
        vmax = float(vmax)
        
        # Store the center value for colormap normalization
        self.center = center
        
        if center is not None:
            center = float(center)
            # For centered colormaps, ensure proper ordering of vmin, center, and vmax
            values = [vmin, center, vmax]
            values.sort()  # Sort in ascending order
            self.vmin, self.center, self.vmax = values
        else:
            self.vmin = vmin
            self.vmax = vmax
        
        # Set default colormap if not provided
        if cmap is None:
            if center is not None:
                self.cmap = "RdBu_r"
            else:
                self.cmap = "viridis"
        else:
            self.cmap = cmap

        self.xticklabel_rotation = xticklabel_rotation
        self.xticklabel_ha = xticklabel_ha

    def _plot_size_legend(self, size_legend_ax):
        """Plot the size legend in a style similar to Scanpy's dotplot."""
        # Get unique marker sizes
        unique_sizes = np.unique(self.marker_size)
        if len(unique_sizes) <= 1:
            return

        # Create a horizontal bar with multiple dots
        n_dots = 3  # Number of dots to show in the legend
        size_range = np.linspace(unique_sizes[0], unique_sizes[-1], n_dots)
        
        # Plot size bar
        size_legend_ax.scatter(
            np.arange(len(size_range)),
            np.repeat(0, len(size_range)),
            s=size_range,
            c='gray',
            edgecolor='black',
            linewidth=0.5,
            zorder=100
        )
        
        # Set x-ticks and labels
        size_legend_ax.set_xticks(np.arange(len(size_range)))
        labels = [f"{x:.0f}" for x in size_range]
        size_legend_ax.set_xticklabels(labels, fontsize='small', rotation=45)
        
        # Remove y-ticks and labels
        size_legend_ax.tick_params(axis='y', left=False, labelleft=False)
        
        # Remove surrounding lines
        size_legend_ax.spines['right'].set_visible(False)
        size_legend_ax.spines['top'].set_visible(False)
        size_legend_ax.spines['left'].set_visible(False)
        size_legend_ax.spines['bottom'].set_visible(False)
        
        # Explicitly remove grid
        size_legend_ax.grid(False)
        
        # Set title
        size_legend_ax.set_title(self.size_title, size='small', pad=10)
        
        # Adjust limits with padding
        xmin, xmax = size_legend_ax.get_xlim()
        size_legend_ax.set_xlim(xmin - 0.2, xmax + 0.2)
        
        # Set y-axis limits to be symmetric around 0
        y_range = 0.5  # Fixed small range for y-axis
        size_legend_ax.set_ylim(-y_range, y_range)

    def plot(self, ax, cax, kws, show_legends=True, legend_spacing=0.15, inner_border=True, show_grid=False):
        """Draw the scattermap on the provided Axes."""
        # Create a new figure with two subplots
        fig = ax.figure
        gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 1], hspace=0.05, wspace=0.05)
        
        # Create the main plot axis
        main_ax = fig.add_subplot(gs[:, 0])

        # Set white background for all axes and figure
        main_ax.set_facecolor('white')
        fig.patch.set_facecolor('white')

        # Remove all the Axes spines from main plot
        despine(ax=main_ax, left=True, bottom=True)

        # Remove the original axis since we're not using it
        ax.remove()

        # Draw the heatmap in the main axis
        data = self.plot_data
        range_y = np.arange(data.shape[0], dtype=int) + 0.5
        range_x = np.arange(data.shape[1], dtype=int) + 0.5
        x, y = np.meshgrid(range_x, range_y)

        # Use TwoSlopeNorm for centered colormaps, otherwise use Normalize
        if self.center is not None:
            norm = plt.cm.colors.TwoSlopeNorm(vmin=self.vmin, vcenter=self.center, vmax=self.vmax)
        else:
            norm = plt.cm.colors.Normalize(vmin=self.vmin, vmax=self.vmax)

        hmap = main_ax.scatter(x, y,
                          c=data,
                          marker=self.marker,
                          cmap=self.cmap,
                          norm=norm,
                          s=self.marker_size, **kws)
        
        self.hmap = hmap
        # Set the axis limits
        main_ax.set(xlim=(0, self.data.shape[1]), ylim=(0, self.data.shape[0]))
        
        # Add a light grey grid if requested
        if show_grid:
            main_ax.grid(True, linestyle='-', color='lightgrey', alpha=0.5, zorder=0)
        
        # Add a 1-point inner border to the main plot if requested
        if inner_border:
            from matplotlib.patches import Rectangle
            border = Rectangle((0, 0), self.data.shape[1], self.data.shape[0], 
                              fill=False, edgecolor='black', linewidth=1, zorder=100)
            main_ax.add_patch(border)
        
        # Add size legend if enabled and marker sizes vary
        legend_width = 0.15  # Default width
        if show_legends and self.size_legend and not isinstance(self.marker_size, (int, float)):
            # Check if we have more than one unique size
            unique_sizes = np.unique(self.marker_size)
            if len(unique_sizes) > 1:
                # Get parameters from size_legend_kws
                legend_width = self.size_legend_kws.get('width', 0.15)
                
                # Create a second axes for the size legend with 1:1.2 aspect ratio (width:height)
                size_ax = fig.add_axes([0.82, 0.1, legend_width, legend_width * 0.25])
                size_ax.set_facecolor('white')  # Set white background for size legend
                self._plot_size_legend(size_ax)
                # Remove y-axis from size legend
                size_ax.tick_params(axis='y', left=False, labelleft=False)

        # Possibly add a colorbar
        if show_legends and self.cbar:
            # Adjust colorbar position to be on top of the size legend
            cbar_ax = fig.add_axes([0.82, 0.35, legend_width, 0.02])
            cbar_ax.set_facecolor('white')  # Set white background for colorbar
            
            # Create a copy of cbar_kws without location parameter and labels
            cbar_kws_copy = self.cbar_kws.copy()
            cbar_kws_copy.pop('location', None)
            cbar_kws_copy.pop('xlabel', None)  # Remove xlabel from kwargs
            cbar_kws_copy.pop('ylabel', None)  # Remove ylabel from kwargs
            
            cb = fig.colorbar(hmap, cax=cbar_ax, ax=main_ax, orientation='horizontal', **cbar_kws_copy)
            cb.outline.set_linewidth(0)
            cb.ax.set_title(self.cbar_title, pad=10)
            
            # Show only min and max values on the colorbar
            cb.ax.tick_params(axis='x', bottom=True, labelbottom=True)
            cb.ax.tick_params(axis='y', left=False, labelleft=False)
            cb.set_ticks([self.vmin, self.vmax])
            cb.set_ticklabels([f'{self.vmin:.1f}', f'{self.vmax:.1f}'])
            
            # Remove all other ticks from colorbar
            cb.ax.tick_params(axis='x', which='minor', bottom=False)
            cb.ax.tick_params(axis='x', which='major', bottom=True)
            
            # If rasterized is passed to pcolormesh, also rasterize the
            # colorbar to avoid white lines on the PDF rendering
            if kws.get('rasterized', False):
                cb.solids.set_rasterized(True)

        # Add row and column labels to main plot
        if isinstance(self.xticks, str) and self.xticks == "auto":
            xticks, xticklabels = self._auto_ticks(main_ax, self.xticklabels, 0)
        else:
            xticks, xticklabels = self.xticks, self.xticklabels

        if isinstance(self.yticks, str) and self.yticks == "auto":
            yticks, yticklabels = self._auto_ticks(main_ax, self.yticklabels, 1)
        else:
            yticks, yticklabels = self.yticks, self.yticklabels

        main_ax.set(xticks=xticks, yticks=yticks)
        xtl = main_ax.set_xticklabels(xticklabels)
        ytl = main_ax.set_yticklabels(yticklabels, rotation="horizontal")

        # Possibly rotate them if they overlap
        fig.draw(fig.canvas.get_renderer())
        if axis_ticklabels_overlap(xtl):
            plt.setp(xtl, rotation=self.xticklabel_rotation, ha=self.xticklabel_ha)
        else:
            # Apply the rotation even if labels don't overlap
            plt.setp(xtl, rotation=self.xticklabel_rotation, ha=self.xticklabel_ha)
        if axis_ticklabels_overlap(ytl):
            plt.setp(ytl, rotation="horizontal")

        # Remove only axis labels from main plot
        main_ax.set_xlabel('')
        main_ax.set_ylabel('')

        # Annotate the cells with the formatted values
        if self.annot:
            # Get the data values and format them
            data = hmap.get_array()
            data = data.reshape(self.data.shape)
            
            # Create text annotations
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    # Skip masked values
                    if self.mask is not None and self.mask[i, j]:
                        continue
                        
                    # Format the value
                    val = data[i, j]
                    if isinstance(self.annot, bool):
                        text = f"{val:.2f}"
                    elif callable(self.annot):
                        text = self.annot(val)
                    else:
                        text = str(self.annot[i, j])
                    
                    # Add text annotation
                    text_color = 'white' if val < 0 else 'black'
                    main_ax.text(j, i, text,
                               horizontalalignment='center',
                               verticalalignment='center',
                               color=text_color)

        # Invert the y axis to show the plot in matrix form
        main_ax.invert_yaxis()
        
        return(hmap)

def scattermap(data,
               marker='o',
               marker_size=100,
               vmin=None, vmax=None, cmap=None, center=None, robust=False,
               linewidths=0, linecolor="white",
               cbar=True, cbar_kws=None, cbar_ax=None,
               square=False, xticklabels="auto", yticklabels="auto",
               mask=None, ax=None, fig_size=None,
               size_legend=True, size_legend_kws=None,
               show_legends=True, legend_spacing=0.15,
               col_order=None, save_pdf=None,
               size_title='-log10(FDR)',
               vcenter=None,
               inner_border=True,
               show_grid=False,
               xticklabel_rotation=0,
               xticklabel_ha='right',
               **kwargs):
    """Plot rectangular data as a color-encoded matrix.

    This function is similar to `sns.heatmap`, as it is an Axes-level function that will draw the
    heatmap into the currently-active Axes if none is provided to the ``ax`` argument.

    The main difference is that instead of drawing an actual heatmap with filled squares,
    this function will use the `plt.scatter` behind the scenes to draw a scatterplot-heatmap.

    The default is set to plot a grid of circles, however this can be changed via `marker`
    parameter.

    Parameters
    ----------
    data : rectangular dataset
        2D dataset that can be coerced into an ndarray. If a Pandas DataFrame
        is provided, the index/column information will be used to label the
        columns and rows.
    marker: string, optional
        Marker to use: any marker that `pyplot.scatter` supports. Defaults to circle.
    marker_size: int or rectangular dataset
        Either an integer to set the marker size of all data points to,
        or a 2D dataset (like in `data`) that sets individual point sizes.
        Defaults to 100.
    vmin, vmax : floats, optional
        Values to anchor the colormap, otherwise they are inferred from the
        data and other keyword arguments.
    cmap : matplotlib colormap name or object, or list of colors, optional
        The mapping from data values to color space. If not provided, the
        default will depend on whether ``center`` is set.
    center : float, optional
        The value at which to center the colormap when plotting divergant data.
        Using this parameter will change the default ``cmap`` if none is
        specified.
    vcenter : float, optional
        The value at which to center the colormap for diverging colormaps.
        This is particularly useful for colormaps like 'bwr' where you want
        zero to be the white point. If not provided, the center will be
        calculated as (vmin + vmax) / 2.
    robust : bool, optional
        If True and ``vmin`` or ``vmax`` are absent, the colormap range is
        computed with robust quantiles instead of the extreme values.
    linewidths : float, optional
        Width of the border lines that will surround the markers
    linecolor : color, optional
        Color of the border lines to the markers
    cbar : boolean, optional
        Whether to draw a colorbar.
    cbar_kws : dict of key, value mappings, optional
        Keyword arguments for `fig.colorbar`.
    cbar_ax : matplotlib Axes, optional
        Axes in which to draw the colorbar, otherwise take space from the
        main Axes.
    square : boolean, optional
        If True, set the Axes aspect to "equal" so each cell will be
        square-shaped.
    xticklabels, yticklabels : "auto", bool, list-like, or int, optional
        If True, plot the column names of the dataframe. If False, don't plot
        the column names. If list-like, plot these alternate labels as the
        xticklabels. If an integer, use the column names but plot only every
        n label. If "auto", try to densely plot non-overlapping labels.
    mask : boolean array or DataFrame, optional
        If passed, data will not be shown in cells where ``mask`` is True.
        Cells with missing values are automatically masked.
    ax : matplotlib Axes, optional
        Axes in which to draw the plot, otherwise use the currently-active
        Axes.
    fig_size : tuple, optional
        Tuple specifying the figure size (width, height).
    size_legend : bool, optional
        Whether to show a legend for marker sizes. Only applies when marker_size is not uniform.
    size_legend_kws : dict, optional
        Keyword arguments for the size legend.
    show_legends : bool, optional
        Whether to show both colorbar and size legend. If False, neither will be shown.
    legend_spacing : float, optional
        The spacing between the main plot and the legends (colorbar and size legend).
        Default is 0.15.
    col_order : list, optional
        List of column names to specify the order of columns in the plot.
        Must contain all column names from the data.
    save_pdf : str, optional
        If provided, save the plot as an editable PDF with the specified filename.
        The PDF will be compatible with Adobe Illustrator and maintain all fonts.
    size_title : str, optional
        Title for the size legend.
    inner_border : bool, optional
        Whether to add a 1-point inner border around the plot. Default is True.
    show_grid : bool, optional
        Whether to show a light grey grid inside the plot. Default is False.
    xticklabel_rotation : int, optional
        Rotation angle in degrees for x-axis labels. Default is 0.
    kwargs : other keyword arguments
        All other keyword arguments are passed to ``ax.pcolormesh``.

    Returns
    -------
    ax : matplotlib Axes
        Axes object with the heatmap.

    See also
    --------
    clustermap : Plot a matrix using hierachical clustering to arrange the
                 rows and columns.
    """
    # Set figure size if provided
    if fig_size is not None:
        plt.figure(figsize=fig_size)

    # Reorder columns if specified
    if col_order is not None:
        if isinstance(data, pd.DataFrame):
            data = data[col_order]
        else:
            data = pd.DataFrame(data).iloc[:, [list(data.columns).index(col) for col in col_order]]

    # Use vcenter if provided, otherwise use center
    if vcenter is not None:
        center = vcenter

    # Initialize the plotter object
    plotter = _ScatterMapper(data,
                             marker, marker_size,
                             vmin, vmax, cmap, center, robust,
                             cbar, cbar_kws, xticklabels,
                             yticklabels, mask,
                             size_legend, size_legend_kws,
                             size_title,
                             xticklabel_rotation,
                             xticklabel_ha)

    # Add the pcolormesh kwargs here
    kwargs["linewidths"] = linewidths
    kwargs["edgecolor"] = linecolor

    # Draw the plot and return the Axes
    if ax is None:
        ax = plt.gca()
    ax.set_aspect("auto")  # Change aspect ratio to be more flexible
    plotter.plot(ax, cbar_ax, kwargs, show_legends, legend_spacing, inner_border, show_grid)
    
    # Save as editable PDF if requested
    if save_pdf is not None:
        # Ensure all text elements are properly rendered
        for ax in plt.gcf().get_axes():
            for text in ax.texts:
                text.set_path_effects([plt.matplotlib.patheffects.Normal()])
        
        # Save with proper font handling
        plt.savefig(save_pdf, 
                   format='pdf',
                   dpi=300,
                   bbox_inches='tight',
                   pad_inches=0.1,
                   transparent=True,
                   metadata={'Creator': 'Matplotlib', 'Producer': 'Matplotlib'},
                   facecolor='white',
                   edgecolor='none')
    
    return ax, plotter