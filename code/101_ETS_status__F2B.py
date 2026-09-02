import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

def bot_bar_plot(
    data,
    palette=None,
    lrot=0,
    figsize=(5, 5),
    title='',
    ax=None,
    order=None,
    stars=False,
    percent=False,
    pvalue=False,
    p_digits=5,
    legend=True,
    xl=True,
    offset=-0.1,
    linewidth=0,
    align='center',
    bar_width=0.9,
    edgecolor=None,
    hide_grid=True,
    draw_horizontal=False,
    plot_all_borders=True,
):
    """
    Plot a stacked bar plot based on contingency table
    :param data: pd.DataFrame, contingency table for plotting. Each element of index corresponds to a bar.
    :param palette: dict, palette for plotting. Keys are unique values from groups, entries are color hexes
    :param lrot: float, rotation angle of bar labels in degrees
    :param figsize: (float, float), figure size in inches
    :param title: str, plot title
    :param ax: matplotlib axis, axis to plot on
    :param order: list, what order to plot the stacks of each bar in. Contains column labels of "data"
    :param stars: bool, whether to use the star notation for p value instead of numerical value
    :param percent: bool, whether to normalize each bar to 1
    :param pvalue: bool, whether to add the p value (chi2 contingency test) to the plot title.
    :param p_digits: int, number of digits to round the p value to
    :param legend: bool, whether to plot the legend
    :param xl: bool, whether to plot bar labels (on x axis for horizontal plot, on y axis for vertical plot)
    :param hide_grid: bool, whether to hide grid on plot
    :param draw_horizontal: bool, whether to draw horizontal bot bar plot
    :param plot_all_borders: bool, whether to plot top and right border
    :return: matplotlib axis
    """
    from matplotlib.ticker import FuncFormatter

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if pvalue:
        from scipy.stats import chi2_contingency

        chi2_test_data = chi2_contingency(data)
        p = chi2_test_data[1]
        if title is not False:
            title += '\n' + get_pvalue_string(p, p_digits, stars=stars)

    if percent:
        c_data = data.apply(lambda x: x * 1.0 / x.sum(), axis=1)
        if title:
            title = '% ' + title
        ax.set_ylim(0, 1)
    else:
        c_data = data

    c_data.columns = [str(x) for x in c_data.columns]

    if order is None:
        order = c_data.columns
    else:
        order = [str(x) for x in order]

    if palette is None:
        c_palette = lin_colors(pd.Series(order))

        if len(order) == 1:
            c_palette = {order[0]: blue_color}
    else:
        c_palette = {str(k): v for k, v in palette.items()}

    if edgecolor is not None:
        edgecolor = [edgecolor] * len(c_data)

    kind_type = 'bar'
    if draw_horizontal:
        kind_type = 'barh'

    c_data[order].plot(
        kind=kind_type,
        stacked=True,
        position=offset,
        width=bar_width,
        color=pd.Series(order).map(c_palette).values,
        ax=ax,
        linewidth=linewidth,
        align=align,
        edgecolor=edgecolor,
    )

    ax = bot_bar_plot_prettify_axis(ax, c_data, legend, draw_horizontal, xl, lrot, title, hide_grid, plot_all_borders)

    if percent:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y)))

    return ax

def bot_bar_plot_prettify_axis(ax, c_data, legend, draw_horizontal, xl, lrot, title, hide_grid, plot_all_borders):
    """
    Change some properties of bot_bar_plot ax
    :return: prettified axis
    """

    if legend:
        ax.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.1)
    else:
        ax.legend_.remove()

    if not draw_horizontal:
        ax.set_xticks(np.arange(len(c_data.index)) + 0.5)
        if xl:
            ax.set_xticklabels(c_data.index, rotation=lrot)
        else:
            ax.set_xticklabels([])
    else:
        ax.set_yticks(np.arange(len(c_data.index)) + 0.5)
        if xl:
            ax.set_yticklabels(c_data.index, rotation=lrot)
        else:
            ax.set_yticklabels([])

    if title is not False:
        ax.set_title(title)

    if hide_grid:
        ax.grid(False)

    sns.despine(ax=ax)

    if plot_all_borders:
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)

    return ax


df = pd.read_csv('~/ETS_events.tsv', 
                 sep='\t', 
                 header=0, 
                 index_col=0, 
                 comment=None)

pal = {'RNA': '#3d5a80', 'CNA': '#e0fbfc', 'SV': '#98c1d9',
       'RNA and SV': '#e7b4a5', 'CNA and RNA': '#ee6c4d',
       'CNA and SV': '#6b0f1a', 'all three': '#2a211b'}


fig, ax = plt.subplots(figsize=(2.5, 3))

bot_bar_plot(df, palette=pal, percent=True, plot_all_borders=False, bar_width=0.8, ax=ax)

plt.xticks(fontsize=7, fontname='Arial')
plt.yticks(fontsize=7, fontname='Arial')
plt.xlabel('ETS-fusion gene', fontsize=7, fontname='Arial')
plt.ylabel('ETS-positive samples (%)', fontsize=7, fontname='Arial')

plot_filename = "~/ETS_events_perc.png"
plt.savefig(plot_filename, format="png", pad_inches=0.3, bbox_inches='tight', dpi=1200, facecolor='white')