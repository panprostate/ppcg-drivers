#!/usr/bin/env python
from __future__ import division, print_function
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import re
import glob

base_dir = "."
sv_dir = "./svs_per_pid/"

df_meta = pd.read_csv(base_dir + "/metadata/WGS_master_tracking_sheet.tsv", sep="\t")
df_meta_primary_select = df_meta.query("selected_one_sample_per_donor == True and SampleType == 'TUM'")
wgs_assayid = df_meta_primary_select["WGS_AssayID"].unique()

oncogenes = ["ERG", "ETV4", "ETV1", "FLI1", "SKIL", "MYC", "ETV5", "AR", "ETS1", "ELK4"]

df_collect = pd.DataFrame()
min_sv_distance = 25000
max_gene_dist = 200000

for i in glob.glob(sv_dir + "*_DNA.ppcg_consensus_annotated.somatic.sv.genes.txt"):
    wgs_assay = re.split(r".ppcg", i.split("/")[-1])[0]
    ppcg_donor_id = re.split(r"[a-z]_DNA", i.split("/")[-1])[0]

    if wgs_assay not in wgs_assayid:
        continue
    df_i_raw = pd.read_csv(i, sep="\t")
    df_i = df_i_raw[((df_i_raw["chrom1"] != df_i_raw["chrom2"]) \
                    | ((df_i_raw["start2"] - df_i_raw["start1"]) > min_sv_distance)) \
                    & (df_i_raw["dist_gene1"].abs() < max_gene_dist) \
                    & (df_i_raw["dist_gene2"].abs() < max_gene_dist)]
    if df_i.empty:
        continue
    df_oncogene_filter = df_i[(df_i["name_gene1"].isin(oncogenes)) | (df_i["name_gene2"].isin(oncogenes))]
    if df_oncogene_filter.empty:
        continue
    df_oncogene_filter = df_oncogene_filter.assign(merge=df_oncogene_filter.\
                 apply(lambda x: x["name_gene1"] + ":" + x["name_gene2"] if x["name_gene1"] != x["name_gene2"] else np.nan, axis=1))

    all_genes = set([i for i in set(df_oncogene_filter["merge"].dropna().values.ravel()) if i != "."])
    partner_gene = set([j for i in all_genes for j in i.split(":") if j not in oncogenes])
    df_tmp = pd.DataFrame({"donor_id": ppcg_donor_id, "partner_genes": ','.join(map(str, partner_gene))}, index=[ppcg_donor_id])
    df_collect = pd.concat([df_collect, df_tmp])

## get number of SV oncogene partner genes
partner_genes = set([j for i in df_collect.partner_genes.unique() for j in i.split(",")])
partner_gene_count = dict()
for partner in partner_genes:
    partner_gene_count[partner] = len(df_collect[df_collect.partner_genes.str.contains(partner)].donor_id.unique())

slc_partner_gene_count = {k: v for k, v in partner_gene_count.items() if "SLC" in k}

df_slc_order_occur = pd.DataFrame(slc_partner_gene_count, index=["n"]).T

####################################################################################################
### GTEX ###
gtex_dir = "/maps/projects/bricsoftware/data/genome_reference/expression/hg38/gtex"
df_tissue = pd.read_csv(gtex_dir + "/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz", sep="\t", skiprows=2)

def threshold_rna(df, variance_threshold=0.2, expression_threshold=1):
    df_hivar = df[df.iloc[:, 2:].var(axis=1) > variance_threshold]
    df_hivar_exp = df_hivar[df_hivar.iloc[:, 2:].mean(axis=1) > expression_threshold]
    return df_hivar_exp

df_tissue_hivar_exp = threshold_rna(df_tissue, variance_threshold=0.2, expression_threshold=1)

df_tissue_slc = df_tissue_hivar_exp.query("Description.str.contains('^SLC')", engine='python').copy()
df_tissue_slc["median_expression"] = df_tissue_slc.iloc[:, 2:].median(axis=1)

## fold change only for gtex data - compare to average of other tissues
df_tissue_slc["prostate_fc"] = df_tissue_slc["Prostate"].apply(lambda x: max(x, 1)) / df_tissue_slc["median_expression"].apply(lambda x: max(x, 1))
df_tissue_slc["prostate_fc_log2"] = np.log2(df_tissue_slc["prostate_fc"])

df_tissue_slc_sorted = df_tissue_slc.sort_values(by="prostate_fc", ascending=False)

df_tissue_slc_sorted = pd.merge(df_tissue_slc_sorted, df_slc_order_occur, left_on="Description", right_index=True, how="left")
df_tissue_slc_sorted["n"] = df_tissue_slc_sorted["n"].fillna(0)

## lollipop plot
plt.rcParams.update({'figure.autolayout': True, 'font.size': 8, 'font.family': 'sans-serif'})

genes2color = df_slc_order_occur.index.tolist()

df_2_plot = df_tissue_slc_sorted.copy()
plotname = "prostate_slc_expression_break_lollipop.pdf"

figsize = (7, 4)
fig, ax = plt.subplots(figsize=figsize, dpi=300)  # Increase the dpi for higher resolution
ax.tick_params(axis='y', labelsize=7)
ax.vlines(x=df_2_plot["Description"], ymin=0, ymax=df_2_plot["prostate_fc_log2"], color='grey', linewidth=0.5, alpha=0.5)
ax.scatter(x=df_2_plot["Description"], y=df_2_plot["prostate_fc_log2"], color='grey', s=0.1, alpha=0.5)
ax.scatter(x=df_2_plot["Description"], y=df_2_plot["prostate_fc_log2"], color='firebrick', s=df_2_plot["n"] * 4)
ax.set_xlim(left=-3)  # Adjust the x-axis limit to bring points closer to the y-axis
for i, row in df_2_plot.iterrows():
    if row["Description"] in genes2color:
        ax.text(row["Description"], row["prostate_fc_log2"] + 0.6,
                row["Description"],
                rotation=45,
                fontsize=6,
                ha='left', va='center')
ax.set_xticks([])
# Add legend with annotated size (single legend, laid out horizontally)
for size in [1, 10, 20, 30]:
    ax.scatter([], [], c='firebrick', s=size * 4, label=str(size))
legend = ax.legend(scatterpoints=1, frameon=False, ncol=4, columnspacing=1, handletextpad=0.3,
                    title='', bbox_to_anchor=(0.55, 1), loc='upper left')
legend.set_title('Number of SLC enhancer\nrearranged tumours')
sns.despine()
# Add labels and title
ax.set_xlabel('SLC genes (sorted by expression fold-change in prostate)')
ax.set_ylabel('Prostate fold-change gene expression (log2)')
ax.set_title('SLC gene expression in prostate vs other tissues')
plt.savefig("/figures/" + plotname, dpi=600)
plt.close()
