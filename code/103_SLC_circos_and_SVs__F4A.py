import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
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




def plot_genes_with_chromhmm_bigwig(gene_data, chromhmm_data, svs_data, svs1_data, bigwig_files, 
                                    hijack_start, hijack_end,
                                    gene1_name=None, gene2_name=None, sv_bin_size=10000):
    
    """
    Plot gene structures, chromatin states, structural variations, and enhancer hijacking region
    along with H3K4me3, H3K27ac, and BRD4 ChIP-seq signals using bigWig files.

    :param gene_data: pd.DataFrame, containing gene annotations with columns ['gene', 'chrom', 'start', 'end', 'type'].
    :param chromhmm_data: pd.DataFrame, containing ChromHMM chromatin state annotations with columns ['chrom', 'start', 'end', 'state'].
    :param svs_data: pd.DataFrame, containing structural variants with columns ['chrom1', 'start1', 'chrom2', 'start2'].
    :param svs1_data: pd.DataFrame, subset of structural variants with the same format as svs_data.
    :param bigwig_files: list of str, paths to three bigWig files for H3K4me3, H3K27ac, and BRD4 signals.
    :param hijack_start: int, start position of enhancer hijacking region from CESAM results.
    :param hijack_end: int, end position of enhancer hijacking region from CESAM results.
    :param gene1_name: str, name of the first gene to be plotted.
    :param gene2_name: str, name of the second gene to be plotted (optional).
    :param sv_bin_size: int, bin size in base pairs for aggregating structural variant counts (default: 10,000 bp).
    :return: None, displays the plot.
    """
 
    if bigwig_files is None or len(bigwig_files) != 3:
        print("Please provide a list of three bigWig file paths for H3K4me3, H3K27ac, and BRD4.")
        return

    bw_H3K4me3 = pyBigWig.open(bigwig_files[0])
    bw_H3K27ac = pyBigWig.open(bigwig_files[1])
    bw_BRD4 = pyBigWig.open(bigwig_files[2])

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
    plot_start = min(g[1]['start'].min() for g in genes) - 80000
    plot_end = max(g[1]['end'].max() for g in genes) + 80000

    chromhmm_filtered = chromhmm_data[
        (chromhmm_data['chrom'] == gene_chrom) &
        (chromhmm_data['start'] >= plot_start) &
        (chromhmm_data['end'] <= plot_end)
    ]

    svs1_filtered_1 = svs1_data[
        (svs1_data['chrom1'] == gene_chrom) &
        (svs1_data['start1'] >= plot_start) &
        (svs1_data['start1'] <= plot_end)
    ]

    svs1_filtered_2 = svs1_data[
        (svs1_data['chrom2'] == gene_chrom) &
        (svs1_data['start2'] >= plot_start) &
        (svs1_data['start2'] <= plot_end)
    ]

    all_svs_positions = np.concatenate([svs1_filtered_1['start1'].values, svs1_filtered_2['start2'].values])

    bin_edges = np.arange(plot_start, plot_end, sv_bin_size)
    svs1_hist, _ = np.histogram(all_svs_positions, bins=bin_edges)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    svs_filtered = svs_data[(svs_data['chrom1'] == gene_chrom) | (svs_data['chrom2'] == gene_chrom)]
    svs_start = hijack_start #45789579 #svs_filtered['start1'].min() #hijacking region
    svs_end = hijack_end #45812683 #svs_filtered['end1'].max() hijacking region

    fig, axs = plt.subplots(7, 1, figsize=(3, 3), sharex=True,
                            gridspec_kw={'height_ratios': [0.1, 0.05, 0.05, 0.15, 0.05, 0.05, 0.05]})
    
    fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05, hspace=0)

    y_positions = [0.7, 0.3]
    for idx, (gene_name, gene_filtered) in enumerate(genes):
        y_position = y_positions[idx]
        axs[0].hlines(y=y_position, xmin=gene_filtered['start'].min() / 1e6,
                      xmax=gene_filtered['end'].max() / 1e6, color='black', linewidth=1)
        for _, exon in gene_filtered[gene_filtered['type'] == 'exon'].iterrows():
            axs[0].add_patch(Rectangle(
                (exon['start'] / 1e6, y_position - 0.03),
                (exon['end'] - exon['start']) / 1e6,
                0.06,
                color='black'
            ))
        
        axs[0].text((gene_filtered['start'].min() + gene_filtered['end'].max()) / 2 / 1e6, 
                    y_position - 0.08, 
                    gene_name, ha='center', va='top', fontsize=7, style='italic')

    axs[0].set_ylim(0, 1)
    axs[0].set_yticks([])
    axs[0].set_xlim(plot_start / 1e6, plot_end / 1e6)

    #SV counts
    axs[1].bar(bin_centers / 1e6, svs1_hist, width=sv_bin_size / 1e6, color='black', alpha=1, edgecolor='black')
    axs[1].set_ylabel('SV counts', fontsize=7, rotation=0, labelpad=15)
    axs[1].set_ylim(0, max(svs1_hist) * 1.1)
    axs[1].tick_params(axis='y', labelsize=7)

    #ChromHMM states
    if not chromhmm_filtered.empty:
        for _, row in chromhmm_filtered.iterrows():
            axs[2].add_patch(Rectangle(
                (row['start'] / 1e6, 0.2),
                (row['end'] - row['start']) / 1e6,
                0.6,
                color='#ffc100', alpha=0.8
            ))
        axs[2].set_ylabel('ChromHMM', fontsize=7, rotation=0, labelpad=15)
        axs[2].set_yticks([])
    else:
        print("⚠ No ChromHMM data found in the selected region!")

    #Enhancer hijacking region
    axs[3].hlines(0.5, svs_start / 1e6, svs_end / 1e6, color='red', linewidth=5)
    axs[3].set_ylabel('Enhancer hijacking', fontsize=7, rotation=0, labelpad=15)
    axs[3].set_yticks([])

    def plot_bigwig_data_peaks(bw, ax, color, label, chrom, start, end, bin_size=500):
        values = bw.values(chrom, start, end, numpy=True)
        values = np.nan_to_num(values)  # Convert NaNs to zero
        if len(values) > 0:
            binned_values = []
            binned_positions = []
            for i in range(0, len(values), bin_size):
                bin_start = start + i
                bin_end = min(bin_start + bin_size, end)
                binned_positions.append((bin_start + bin_end) / 2)
                binned_values.append(np.nanmean(values[i:i + bin_size]))

            ax.fill_between(np.array(binned_positions) / 1e6, binned_values, color=color, alpha=1)
            ax.set_ylabel(label, fontsize=7, rotation=0, labelpad=15)
            ax.set_yticks([])
            ax.set_xlim(start / 1e6, end / 1e6)
        else:
            print(f"⚠ Warning: No data in BigWig track {label}!")

    plot_bigwig_data_peaks(bw_H3K4me3, axs[4], "darkgreen", "H3K4me3", gene_chrom, plot_start, plot_end, bin_size=200)
    plot_bigwig_data_peaks(bw_H3K27ac, axs[5], "darkorange", "H3K27ac", gene_chrom, plot_start, plot_end, bin_size=200)
    plot_bigwig_data_peaks(bw_BRD4, axs[6], "brown", "BRD4", gene_chrom, plot_start, plot_end, bin_size=200)
    axs[6].set_xlabel('Position (Mb)', fontsize=7)

    axs[0].spines['bottom'].set_visible(False)
    axs[-1].spines['top'].set_visible(False)
    for ax in axs[1:-1]:
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        
    #plt.savefig(f"/home/dkiriy/data/plots/SLC30A4.png", format="png", bbox_inches='tight', dpi=1200)
    plt.show()

    bw_H3K4me3.close()
    bw_H3K27ac.close()
    bw_BRD4.close()

    
    
# read Pomerantz data    
po = pd.read_csv('~/Pomerantz.bed',
                    sep='\t', 
                    comment=None,
                    index_col=False,
                    header=None)[[0,1,2,3]].rename(columns={0:'chrom',1:'start',2:'end',3:'state'})

po['state'].replace({
   'Active_prostate_lineage-specific_promoter':'Act_PSP',
    'Bivalent_poised_promoter':'Bival_P',
     'Active_non-prostate_lineage_promoter':'Act_nPSP',
      'Active_prostate_lineage-specific_enhancer':'Act_PSE',
       'Active_non-prostate_lineage_enhancer':'Act_nPSE',
        'Primed_non-prostate_lineage-specific_enhancer':'Primed_nPSE',
        'Primed_prostate_lineage_enhancer':'Primed_PSE',
         'Bivalent_poised_enhancer':'Bival_E',
          'Primed_prostate_lineage_enhancer':'Primed_PSE',
           'Heterochromatin':'Het',
            'Repressed_chromatin':'Rep_chrom'
}, inplace=True)

# choose only Enhancer regions
po = po[po.state.isin(['Act_PSE', 
                       'Act_nPSE', 
                      ])]    


# read samples annotation
main_track = pd.read_csv('~/annotation.tsv',
                         sep='\t', 
                         comment=None,
                         index_col=False,
                         header=0,
                         )
main_track = main_track[main_track['selected_one_primary_per_donor']==True]


#read SV data
svs = pd.read_csv('~/SVs.tsv', 
                   sep='\t', 
                   comment=None,
                   index_col=False,
                   header=0
                  )
svs['chrom1'] = 'chr'+svs['chrom1'].astype(str)
svs['chrom2'] = 'chr'+svs['chrom2'].astype(str)
svs = svs[svs.PPCG_Sample_ID.isin(main_track['WGS_AssayID'])]


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

#paths to bigwig files
bigwig_files = [
    "~/H3K4me3.bigWig",
    "~/H3K27ac.bigWig",
    "~/BRD4.bigWig"
]

plot_genes_with_chromhmm_bigwig(
    gene_data=gene_ref, 
    chromhmm_data=po, 
    svs_data=svs,
    svs1_data=svs,
    bigwig_files=bigwig_files,
    hijack_start = 45789379, 
    hijack_end = 45812683,
    gene1_name='SLC30A4',
    sv_bin_size=200
)

plot_genes_with_chromhmm_bigwig(
    gene_data=gene_ref, 
    chromhmm_data=po, 
    svs_data=svs,
    svs1_data=svs,
    bigwig_files=bigwig_files,
    hijack_start = 205617386,
    hijack_end = 205646830,
    gene1_name='SLC45A3',
    gene2_name='ELK4',
    sv_bin_size=200 
)