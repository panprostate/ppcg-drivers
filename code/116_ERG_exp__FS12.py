import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import mannwhitneyu


## Helper functions used by boxplot_with_pvalue

def to_common_samples(objects):
    common_samples = objects[0].index

    for obj in objects[1:]:
        common_samples = common_samples.intersection(obj.index)

    return [
        obj.loc[common_samples]
        for obj in objects
    ]


def define_ax_figsize(ax):
    figure = ax.get_figure()
    bbox = ax.get_window_extent().transformed(
        figure.dpi_scale_trans.inverted()
    )

    return bbox.width, bbox.height


def get_pvalue_string(pvalue, p_digits=3, stars=True):
    if stars:
        if pvalue < 0.0001:
            return "****"
        if pvalue < 0.001:
            return "***"
        if pvalue < 0.01:
            return "**"
        if pvalue < 0.05:
            return "*"

        return "ns"

    return f"p={pvalue:.{p_digits}f}"


## Boxplot function with pairwise Mann–Whitney U tests

def boxplot_with_pvalue(
    data,
    grouping,
    title="",
    ax=None,
    figsize=None,
    swarm=True,
    p_digits=3,
    stars=True,
    violin=False,
    palette=None,
    order=None,
    y_min=None,
    y_max=None,
    s=7,
    p_fontsize=16,
    xlabel=None,
    **kwargs,
):
    cdata, cgrouping = to_common_samples(
        [
            data.dropna(),
            grouping.dropna(),
        ]
    )

    valid_samples = cdata.notna() & cgrouping.notna()

    cdata = cdata.loc[valid_samples]
    cgrouping = cgrouping.loc[valid_samples]

    if cgrouping.nunique() < 2:
        raise ValueError(
            f"Fewer than two groups were provided: "
            f"{cgrouping.nunique()}"
        )

    if order is None:
        order = cgrouping.unique().tolist()
    else:
        order = [
            group_name
            for group_name in order
            if group_name in cgrouping.unique()
        ]

    if ax is None:
        if figsize is None:
            figsize = (1.2 * len(order), 4)

        _, ax = plt.subplots(figsize=figsize)

    if violin:
        sns.violinplot(
            y=cdata,
            x=cgrouping,
            ax=ax,
            palette=palette,
            order=order,
            **kwargs,
        )

        swarm = False

    else:
        sns.boxplot(
            y=cdata,
            x=cgrouping,
            ax=ax,
            palette=palette,
            order=order,
            fliersize=0,
            **kwargs,
        )

    if swarm:
        sns.swarmplot(
            y=cdata,
            x=cgrouping,
            ax=ax,
            color=".25",
            order=order,
            s=s,
        )

    pvalues = []

    for group_1, group_2 in zip(order[:-1], order[1:]):
        samples_group_1 = cgrouping.index[
            cgrouping == group_1
        ]

        samples_group_2 = cgrouping.index[
            cgrouping == group_2
        ]

        if len(samples_group_1) > 0 and len(samples_group_2) > 0:
            try:
                pvalue = mannwhitneyu(
                    cdata.loc[samples_group_1],
                    cdata.loc[samples_group_2],
                    alternative="two-sided",
                ).pvalue

            except ValueError:
                pvalue = 1.0

        else:
            pvalue = 1.0

        pvalues.append(pvalue)

    if y_max is None:
        y_max = cdata.max()

    if y_min is None:
        y_min = cdata.min()

    effective_size = y_max - y_min

    if effective_size == 0:
        effective_size = 1

    plot_y_limits = (
        y_min - effective_size * 0.15,
        y_max + effective_size * 0.20,
    )

    if p_digits > 0 and len(pvalues) > 0:
        pvalue_line_y = y_max + effective_size * 0.05

        if figsize is None:
            figsize = define_ax_figsize(ax)

        pvalue_text_y = (
            pvalue_line_y
            + 0.25 * effective_size / figsize[1]
        )

        for position, pvalue in enumerate(pvalues):
            pvalue_string = get_pvalue_string(
                pvalue=pvalue,
                p_digits=p_digits,
                stars=stars,
            )

            local_text_y = pvalue_text_y

            if pvalue_string == "-":
                local_text_y += (
                    0.1 * effective_size / figsize[1]
                )

            bar_fraction = str(
                0.25
                / 2.0
                / (figsize[0] / float(len(order)))
            )

            ax.annotate(
                "",
                xy=(position + 0.1, pvalue_line_y),
                xycoords="data",
                xytext=(position + 0.9, pvalue_line_y),
                textcoords="data",
                arrowprops={
                    "arrowstyle": "-",
                    "ec": "#000000",
                    "connectionstyle": (
                        f"bar,fraction={bar_fraction}"
                    ),
                },
            )

            ax.text(
                position + 0.5,
                local_text_y,
                pvalue_string,
                fontsize=p_fontsize,
                horizontalalignment="center",
                verticalalignment="center",
            )

    ax.set_title(title)
    ax.set_ylim(plot_y_limits)

    if xlabel is not None:
        ax.set_xlabel(xlabel)

    return ax


## Settings

genes = ["ERG", "ETV1", "ETV4", "ETV5", "FLI1"]

group_names = [
    "WT",
    "sv",
    "del",
    "del/rna",
    "rna",
    "del/sv",
    "rna/sv",
    "sv/del/rna",
]

group_order = [
    "WT",
    "del",
    "rna",
    "del/sv",
    "sv",
    "del/rna",
    "rna/sv",
    "sv/del/rna",
]

palette = {
    "WT": "#bdbdbd",
    "rna": "#3d5a80",
    "del": "#e0fbfc",
    "sv": "#98c1d9",
    "rna/sv": "#e7b4a5",
    "del/rna": "#ee6c4d",
    "del/sv": "#6b0f1a",
    "sv/del/rna": "#2a211b",
}

pretty_labels = {
    "WT": "WT",
    "rna": "RNA",
    "del": "CNA",
    "sv": "SV",
    "rna/sv": "RNA and SV",
    "del/rna": "CNA and RNA",
    "del/sv": "CNA and SV",
    "sv/del/rna": "All three",
}


## Prepare annotation and expression data

ann_plot = ann.copy()
exp_plot = exp.copy()

ann_plot.index = ann_plot.index.astype(str)
exp_plot.columns = exp_plot.columns.astype(str)
exp_plot.index = exp_plot.index.astype(str)

ann_plot = ann_plot.loc[
    ~ann_plot.index.duplicated(keep="first")
].copy()

if exp_plot.columns.duplicated().any():
    exp_plot = exp_plot.T.groupby(level=0).mean().T

if exp_plot.index.duplicated().any():
    exp_plot = exp_plot.groupby(level=0).mean()


## Convert alteration columns to Boolean masks

def make_binary_mask(series):
    numeric_series = pd.to_numeric(
        series,
        errors="coerce",
    )

    return numeric_series.fillna(0).eq(1)


## Assign mutually exclusive alteration groups

def assign_gene_group(annotation, gene):
    rna_col = f"{gene}_rna"
    del_col = f"{gene}_del"
    sv_col = f"{gene}_sv"

    if rna_col in annotation.columns:
        rna_mask = make_binary_mask(
            annotation[rna_col]
        )
    else:
        rna_mask = pd.Series(
            False,
            index=annotation.index,
        )

    if del_col in annotation.columns:
        del_mask = make_binary_mask(
            annotation[del_col]
        )
    else:
        del_mask = pd.Series(
            False,
            index=annotation.index,
        )

    if sv_col in annotation.columns:
        sv_mask = make_binary_mask(
            annotation[sv_col]
        )
    else:
        sv_mask = pd.Series(
            False,
            index=annotation.index,
        )

    gene_group = pd.Series(
        "WT",
        index=annotation.index,
        dtype="object",
        name=f"{gene}_group",
    )

    gene_group.loc[
        rna_mask & ~del_mask & ~sv_mask
    ] = "rna"

    gene_group.loc[
        del_mask & ~rna_mask & ~sv_mask
    ] = "del"

    gene_group.loc[
        sv_mask & ~rna_mask & ~del_mask
    ] = "sv"

    gene_group.loc[
        del_mask & rna_mask & ~sv_mask
    ] = "del/rna"

    gene_group.loc[
        del_mask & ~rna_mask & sv_mask
    ] = "del/sv"

    gene_group.loc[
        ~del_mask & rna_mask & sv_mask
    ] = "rna/sv"

    gene_group.loc[
        del_mask & rna_mask & sv_mask
    ] = "sv/del/rna"

    return pd.Series(
        pd.Categorical(
            gene_group,
            categories=group_order,
            ordered=True,
        ),
        index=annotation.index,
        name=f"{gene}_group",
    )


## Add alteration groups to annotation

for gene in genes:
    ann_plot[f"{gene}_group"] = assign_gene_group(
        annotation=ann_plot,
        gene=gene,
    )


## Store sample IDs and group counts

group_samples = {}

output_df = pd.DataFrame(
    0,
    index=genes,
    columns=group_names,
    dtype=int,
)

for gene in genes:
    group_col = f"{gene}_group"

    group_samples[gene] = {}

    for group_name in group_names:
        sample_ids = ann_plot.index[
            ann_plot[group_col]
            .astype("string")
            .eq(group_name)
        ].tolist()

        group_samples[gene][group_name] = sample_ids

        output_df.loc[
            gene,
            group_name,
        ] = len(sample_ids)


## Check group assignments

for gene in genes:
    number_grouped = output_df.loc[gene].sum()

    if number_grouped != len(ann_plot):
        raise ValueError(
            f"{gene}: {number_grouped} samples were assigned, "
            f"but ann contains {len(ann_plot)} samples."
        )

display(output_df)


## Prepare ERG expression and groups

gene_to_plot = "ERG"

if gene_to_plot not in exp_plot.index:
    raise KeyError(
        f"{gene_to_plot} is not present in exp."
    )

erg_expr = pd.to_numeric(
    exp_plot.loc[gene_to_plot],
    errors="coerce",
)

erg_expr.name = "ERG_expression"

erg_group = ann_plot[
    f"{gene_to_plot}_group"
].copy()

erg_group.name = "ERG_group"


## Match annotation and expression samples

common_samples = erg_expr.index.intersection(
    erg_group.index
)

erg_expr = erg_expr.loc[common_samples]
erg_group = erg_group.loc[common_samples]

valid_samples = (
    erg_expr.notna()
    & erg_group.notna()
)

erg_expr = erg_expr.loc[valid_samples]
erg_group = erg_group.loc[valid_samples]


## Remove empty groups

plot_order = [
    group_name
    for group_name in group_order
    if (
        erg_group.astype("string")
        == group_name
    ).any()
]

plot_palette = {
    group_name: palette[group_name]
    for group_name in plot_order
}


## Count samples in the ERG plot

erg_plot_counts = (
    erg_group
    .astype("string")
    .value_counts()
    .reindex(plot_order)
    .fillna(0)
    .astype(int)
)

display(
    erg_plot_counts.to_frame(
        "Number of samples"
    )
)


## Plot ERG expression

fig, ax = plt.subplots(
    figsize=(5, 4)
)

boxplot_with_pvalue(
    data=erg_expr,
    grouping=erg_group,
    ax=ax,
    order=plot_order,
    palette=plot_palette,
    swarm=True,
    violin=False,
    p_digits=3,
    s=2,
    stars=True,
    xlabel="ERG alteration group",
)

ax.set_ylabel("ERG expression")

ax.set_xticks(
    range(len(plot_order))
)

ax.set_xticklabels(
    [
        (
            f"{pretty_labels[group_name]}"
            f"\n(n={erg_plot_counts[group_name]})"
        )
        for group_name in plot_order
    ],
    rotation=45,
    ha="right",
)

plt.tight_layout()
plt.show()