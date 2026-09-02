import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


## Extract unique donor IDs

def get_donors(values):
    return {
        donor.strip()
        for value in values.dropna().astype(str)
        for donor in value.split(",")
        if donor.strip()
    }


## Plot targetable events

def plot_targetable_events(
    df,
    total_patients=961,
    levels=None,
    figsize=(6, 2.5),
):
    df = df.copy()

    df["level of evidence"] = pd.to_numeric(
        df["level of evidence"],
        errors="coerce",
    )

    if levels is not None:
        if not isinstance(levels, (list, tuple, set)):
            levels = [levels]

        df = df[
            df["level of evidence"].isin(levels)
        ].copy()

    if df.empty:
        raise ValueError("No rows remain after filtering.")

    level_colors = {
        1: "#56a754",
        2: "#5276c8",
        4: "#3f3f3f",
    }

    rows = []

    all_n = len(get_donors(df["donor id"]))

    rows.append(
        {
            "Gene": "All features",
            "Level": 0,
            "Cancer": "All",
            "N": all_n,
            "Percent": all_n / total_patients * 100,
        }
    )

    for gene in df.index.unique():
        gene_df = df.loc[[gene]]
        n = len(get_donors(gene_df["donor id"]))

        levels_gene = (
            gene_df["level of evidence"]
            .dropna()
            .astype(int)
        )

        cancers_gene = (
            gene_df["cancer"]
            .dropna()
            .astype(str)
        )

        rows.append(
            {
                "Gene": gene,
                "Level": (
                    levels_gene.min()
                    if len(levels_gene)
                    else 0
                ),
                "Cancer": (
                    "Prostate cancer"
                    if "Prostate cancer" in cancers_gene.values
                    else cancers_gene.iloc[0]
                    if len(cancers_gene)
                    else "Unknown"
                ),
                "N": n,
                "Percent": n / total_patients * 100,
            }
        )

    plot_df = pd.DataFrame(rows)

    plot_df = pd.concat(
        [
            plot_df.iloc[[0]],
            plot_df.iloc[1:].sort_values(
                "N",
                ascending=False,
            ),
        ],
        ignore_index=True,
    )

    fig, ax = plt.subplots(figsize=figsize)

    x_positions = [
        i * 0.8
        for i in range(len(plot_df))
    ]

    for x, (_, row) in zip(
        x_positions,
        plot_df.iterrows(),
    ):
        hatch = (
            "///"
            if row["Cancer"] == "Prostate cancer"
            else ""
        )

        if row["Gene"] == "All features":
            ax.bar(
                x,
                row["Percent"],
                width=0.7,
                color="#c85252",
            )

            ax.bar(
                x,
                100 - row["Percent"],
                width=0.7,
                bottom=row["Percent"],
                color="lightgrey",
            )

        else:
            ax.bar(
                x,
                row["Percent"],
                width=0.7,
                color=level_colors.get(
                    row["Level"],
                    "lightgrey",
                ),
                hatch=hatch,
            )

        ax.text(
            x,
            row["Percent"] + 1,
            str(row["N"]),
            ha="center",
            va="bottom",
            fontsize=6,
        )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(
        plot_df["Gene"],
        rotation=90,
        fontsize=7,
    )

    ax.set_ylabel(
        "Patients (%)",
        fontsize=7,
    )

    ax.set_ylim(0, 100)
    ax.tick_params(axis="y", labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [
        mpatches.Patch(
            color="#c85252",
            label="With targetable events",
        ),
        mpatches.Patch(
            color="lightgrey",
            label="Without targetable events",
        ),
    ]

    for level in sorted(
        df["level of evidence"]
        .dropna()
        .astype(int)
        .unique()
    ):
        handles.append(
            mpatches.Patch(
                color=level_colors.get(
                    level,
                    "lightgrey",
                ),
                label=f"Level of evidence {level}",
            )
        )

    handles.append(
        mpatches.Patch(
            facecolor="none",
            edgecolor="black",
            hatch="///",
            label="Prostate cancer",
        )
    )

    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1, 1),
        prop={"size": 7},
    )

    plt.tight_layout()

    return fig, ax, plot_df