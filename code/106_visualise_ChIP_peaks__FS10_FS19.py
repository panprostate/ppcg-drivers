import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle, FancyArrow, Arc
import matplotlib.pyplot as plt
import pyBigWig
from matplotlib import rcParams

rcParams['font.family'] = 'Arial'
rcParams['font.size'] = 7


def load_gff_pd(file, limits=None, ordered=False):
    names = [
        "chrom",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "phase",
    ]
    types = [str, str, str, int, int, str, str, str]

    ordered = ordered and limits is None
    if limits is None:
        unlimited = True
        start = end = None
    else:
        unlimited = False
        start, end = limits

    dfcols = []
    with open(file) as f:
        for line in f:
            stripline = line.strip().split("\t")
            if stripline[0][0] == "#":  
                continue
            if unlimited or (int(stripline[4]) > start and int(stripline[3]) < end):
                core_dict = {n: t(s) for n, s, t in zip(names, stripline[:-1], types)}
                attr_dict = dict(
                    [i.strip('"') for i in item.strip().split("=", 1)]
                    for item in stripline[-1].split(";")
                    if item.strip()
                )
                core_dict.update(attr_dict) 
                dfcols.append(core_dict)
            if ordered and int(stripline[3]) > end:
                break

    df = pd.DataFrame(dfcols)
    print("DataFrame created with shape:", df.shape)  
    return df 


def plot_genes_with_loops_bigwig(gene_data, loop_data, bigwig_files, gene1_name=None, gene2_name=None):
    
    """
    Plot gene structures and chromatin loops along with H3K4me3, H3K27ac, CTCF, and BRD4 ChIP-seq signals.

    :param gene_data: pd.DataFrame, containing gene annotations with columns ['gene', 'chrom', 'start', 'end', 'type',    'strand'].
    :param loop_data: pd.DataFrame, containing chromatin loop data with columns ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'].
    :param bigwig_files: list of str, paths to four bigWig files for H3K4me3, H3K27ac, CTCF, and BRD4 signals.
    :param gene1_name: str, name of the first gene to be plotted.
    :param gene2_name: str, name of the second gene to be plotted (optional).
    :return: None, displays the plot.
    """
    
    if bigwig_files is None or len(bigwig_files) != 4:
        print("Please provide a list of four bigWig file paths for H3K4me3, H3K27ac, CTCF, and BRD4.")
        return

    bw_H3K4me3 = pyBigWig.open(bigwig_files[0])
    bw_H3K27ac = pyBigWig.open(bigwig_files[1])
    bw_CTCF = pyBigWig.open(bigwig_files[2])
    bw_BRD4 = pyBigWig.open(bigwig_files[3])

    genes = []
    if gene1_name:
        gene1_filtered = gene_data[gene_data['gene'] == gene1_name]
        if not gene1_filtered.empty:
            genes.append((gene1_name, gene1_filtered))
        else:
            print(f"No data found for gene: {gene1_name}")

    if gene2_name:
        gene2_filtered = gene_data[gene_data['gene'] == gene2_name]
        if not gene2_filtered.empty:
            genes.append((gene2_name, gene2_filtered))
        else:
            print(f"No data found for gene: {gene2_name}")

    if not genes:
        print("No valid gene data provided.")
        return

    gene_chrom = genes[0][1]['chrom'].iloc[0]
    plot_start = min(g[1]['start'].min() for g in genes) - 10000
    plot_end = max(g[1]['end'].max() for g in genes) + 30000

    loops_filtered = loop_data[
        (loop_data['chrom1'] == gene_chrom) &
        (loop_data['chrom2'] == gene_chrom) &
        (loop_data['start1'] >= plot_start) & (loop_data['end1'] <= plot_end) &
        (loop_data['start2'] >= plot_start) & (loop_data['end2'] <= plot_end)
    ]

    fig, axs = plt.subplots(6, 1, figsize=(6, 5), sharex=True, 
                            gridspec_kw={'height_ratios': [0.3, 0.5, 0.3, 0.3, 0.3, 0.3]})
    fig.subplots_adjust(hspace=0.005)

    for _, loop in loops_filtered.iterrows():
        start1, end1 = loop['start1'] / 1e6, loop['end1'] / 1e6
        start2, end2 = loop['start2'] / 1e6, loop['end2'] / 1e6
        mid1, mid2 = (start1 + end1) / 2, (start2 + end2) / 2
        width = abs(mid2 - mid1)

        axs[0].hlines(y=0.1, xmin=start1, xmax=end1, color="black", linewidth=1)
        axs[0].hlines(y=0.1, xmin=start2, xmax=end2, color="black", linewidth=1)
        axs[0].add_patch(Arc(((mid1 + mid2) / 2, 0.1), width, 1.5, theta1=0, theta2=180, edgecolor='black', lw=1, alpha=0.7))
    
    axs[0].set_ylabel('chromatin\nloops', rotation=0)
    axs[0].set_ylim(0, 1.5)
    axs[0].set_yticks([])

    y_positions = [0.7, 0.3]
    for idx, (gene_name, gene_filtered) in enumerate(genes):
        y_position = y_positions[idx]
        axs[1].hlines(y=y_position, xmin=gene_filtered['start'].min() / 1e6, 
                      xmax=gene_filtered['end'].max() / 1e6, color='black', linewidth=1)
        for _, exon in gene_filtered[gene_filtered['type'] == 'exon'].iterrows():
            axs[1].add_patch(Rectangle(
                (exon['start'] / 1e6, y_position - 0.05),  
                (exon['end'] - exon['start']) / 1e6,
                0.1,
                color='black'
            ))
        axs[1].text((gene_filtered['start'].min() + gene_filtered['end'].max()) / 2 / 1e6, y_position - 0.15, 
                    gene_name, ha='center', va='top', fontsize=7, style='italic')
    
    axs[1].set_ylim(0, 1)
    axs[1].set_yticks([])
    axs[1].set_xlim(plot_start / 1e6, plot_end / 1e6)

    bigwig_bin_centers = np.linspace(plot_start, plot_end, plot_end - plot_start) / 1e6

    def plot_bigwig_data(bw, ax, color, label):
        values = bw.values(gene_chrom, plot_start, plot_end, numpy=True)
        values = np.nan_to_num(values)
        ax.fill_between(bigwig_bin_centers, values, color=color, alpha=1)
        ax.set_ylabel(label, rotation=0)
        ax.set_yticks([])
        ax.spines['right'].set_visible(True)

    plot_bigwig_data(bw_H3K4me3, axs[2], 'darkgreen', 'H3K4me3')
    plot_bigwig_data(bw_H3K27ac, axs[3], 'darkorange', 'H3K27ac')
    plot_bigwig_data(bw_CTCF, axs[4], 'purple', 'CTCF')
    plot_bigwig_data(bw_BRD4, axs[5], 'brown', 'BRD4')
    axs[5].set_xlabel('Position (Mb)')

    plt.tight_layout()
    
#     plt.savefig("~/drivers_{}_{}_new.png".format(gene1_name or "", gene2_name or ""), 
#                     format="png", pad_inches=.3, bbox_inches='tight', dpi=600)
    
    plt.show()

    bw_H3K4me3.close()
    bw_H3K27ac.close()
    bw_CTCF.close()
    bw_BRD4.close()

#read loops data
loops = pd.read_csv('~/loops.tsv', 
                    sep='\t', 
                    comment=None,
                    index_col=False,
                    header=0
                    )
loops = loops.iloc[:,:7]  
    
#read GFF files    
gff_df = load_gff_pd('~/gencode.v19.annotation.gff3')

g=['FOXA1', 'MIPOL1', 'NEAT1', 'MALAT1', 'SLC30A4', 'SLC45A3', 'ELK4']
gff_df = gff_df[gff_df.gene_name.notna()]
gff_df = gff_df[gff_df.gene_name.isin(g)]
gff_df = gff_df[['chrom', 'start', 'end', 'gene_name', 
                   'feature','strand', 'gene_id']].rename(columns={'feature':'type'})

tr_ids = ['ENSG00000129514.4', #FOXA1
          'ENSG00000151338.14', #MIPOL1
         'ENSG00000245532.4',  #NEAT1,
          'ENSG00000251562.3', #MALAT1
         'ENSG00000104154.5', #SLC30A4
          'ENSG00000158715.5', #SLC45A3
          'ENSG00000158711.9' #ELK4
         ]
gene_ref = gff_df[gff_df['gene_id'].isin(tr_ids)]
gene_ref = gene_ref.rename(columns={'gene_name':'gene'})
    

# paths to bigwig files
bigwig_files = [
    "/home/dkiriy/data/refs/H3K4me3_GSE96019_ENCFF710JYO_signal_p-value_hg19.bigWig",
    "/home/dkiriy/data/refs/H3K27ac_GSE105424_ENCFF994HZA_signal_p-value_hg19.bigWig",
    "/home/dkiriy/data/refs/CTCF_GSE123210_ENCFF785IOE_signal_p-value_hg19.bigWig",
    "/home/dkiriy/data/refs/BRD4_GSE209886_LNCaP.bigWig"
]

# plotting
plot_genes_with_loops_bigwig(gene_data=gene_ref, loop_data=loops,
                                 gene1_name='NEAT1', 
                                 gene2_name='MALAT1', 
                                 bigwig_files=bigwig_files)

plot_genes_with_loops_bigwig(gene_data=gene_ref, loop_data=loops,
                                 gene1_name='MIPOL1', 
                                 gene2_name='FOXA1', 
                                 bigwig_files=bigwig_files)