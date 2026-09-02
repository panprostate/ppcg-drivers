source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library("gtools")
library("RColorBrewer")

#
# BI driver genes
#
load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.005
MIN_N_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)
PVALUE_CUTOFF = 0.05

load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes = c(cgc2024, prostate_cancer_genes)


#
# biallelic patient lists
#
BI_results = read.delim(
		paste0("DATA_USED__", this_timestamp, "/biallelic/Daria_2024-10-25/PPCG_Drivers/Tables/biallelic_inactivation_pvals_from_curveball_251024.tsv"),
		stringsAsFactors = FALSE)
BI_results$is_cancer_gene = BI_results$gene %in% cancer_genes

BI_patients = lapply(BI_results$pats, function(x) strsplit(gsub("\\'|\\[|\\]", "", x), ", ")[[1]])
names(BI_patients) = BI_results$gene

# exclude mets if present
load(file = pff("tracking_sheet.rsav"))
unique_donor_ids = unique(tracking_sheet$PPCG_Donor_ID)
BI_patients = lapply(BI_patients, intersect, unique_donor_ids)
BI_results$n_hits = lapply(BI_patients, length)

# select non-met patients
which_BI_genes_select = which(BI_results$pvalue < PVALUE_CUTOFF & BI_results$n_hits >= MIN_N_PATIENTS)
BI_results = BI_results[which_BI_genes_select,]
BI_patients = BI_patients[which_BI_genes_select]

BI_results$patients = sapply(BI_patients, paste, collapse = ",")
BI_results = BI_results[, c("gene", "pvalue", "n_hits", "is_cancer_gene", "patients")]		
colnames(BI_results) = 	c("gene", "pval", "n_patients", "is_cancer_gene", "patients")	
save(BI_results, file = pff("BI_results.rsav"))


#
# plot genes grouped by event type
#
MIN_DRIVER_FREQ = 0.01
MIN_N_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)

load(pff("blacklisted_tumor_ids.rsav"))
load(pff("distinct_tumor_ids.rsav"))

genes_to_plot = BI_results[BI_results$is_cancer_gene == TRUE, "gene"]
BI_patient_sets = read.delim(
		paste0("DATA_USED__", this_timestamp, "/", "biallelic/Daria_2024-10-25/PPCG_Drivers/Tables/biallelic_inactivation_251024.tsv"),
		stringsAsFactors = FALSE)
BI_patient_sets = BI_patient_sets[BI_patient_sets$gene %in% genes_to_plot, ]
BI_patient_sets = BI_patient_sets[BI_patient_sets$donor_id %in% unique_donor_ids,]


# keep only 1% mutated genes
gene_mut_count = by(BI_patient_sets$donor_id, BI_patient_sets$gene, function(x) length(unique(x)))
genes_to_select = names(which(gene_mut_count > MIN_N_PATIENTS))
BI_patient_sets = BI_patient_sets[BI_patient_sets$gene %in% genes_to_select,]

# top genes first, most frequent events first
genes_by_muts = names(sort(by(BI_patient_sets$donor_id, BI_patient_sets$gene, function(x) length(unique(x)))))
BI_patient_sets$gene = factor(BI_patient_sets$gene, levels = genes_by_muts)
events_by_muts = names(sort(by(BI_patient_sets$sample_id, BI_patient_sets$hit_type, function(x) length(unique(x)))))
BI_patient_sets$hit_type = factor(BI_patient_sets$hit_type, levels = events_by_muts)

# create colors using Rcolorbrewer for BI event types
bi_colors_set1 = brewer.pal(length(events_by_muts), "Set3")
names(bi_colors_set1) = rev(events_by_muts)

# label pvalues as stars
gene_to_pval = structure(names = BI_results$gene, BI_results$pval)
BI_patient_sets$pval_stars = gene_to_pval[BI_patient_sets$gene]
BI_patient_sets$pval_stars[duplicated(BI_patient_sets$gene)] = NA

# label cancer genes
BI_patient_sets$is_cancer_gene = NA
BI_patient_sets$is_cancer_gene[BI_patient_sets$gene %in% cgc2024] = "cancer gene"
BI_patient_sets$is_cancer_gene[BI_patient_sets$gene %in% prostate_cancer_genes] = "prostate cancer gene"
BI_patient_sets$is_cancer_gene[duplicated(BI_patient_sets$gene)] = NA
cancer_gene_symbol = c("prostate cancer gene" = 1, "cancer gene" = 4)

save(BI_patient_sets, file = pff("BI_patient_sets.rsav"))

plt = ggplot(BI_patient_sets, aes(gene, fill = hit_type)) + 
	geom_bar() +
	coord_flip() + 
	plot_theme() + 
	scale_fill_manual(values = bi_colors_set1) +
	geom_point(aes(gene, y = -5, shape = is_cancer_gene)) + 
	scale_shape_manual(values = cancer_gene_symbol) + 
	geom_text(aes(label = stars.pval(pval_stars)), y = 20) + 
	theme(legend.position = "right") +
	ggtitle("genes with biallelic inactivation", 
			paste0("N=", length(unique(BI_patient_sets$gene)), "; P<", PVALUE_CUTOFF, "; n_pat>=", ceiling(MIN_N_PATIENTS)))
	

fname = pff("figures/BI_gene_list.pdf")
ggsave(plt, file = fname)
file_open_call2(fname)



#
# another plot showing all BI genes
#

MIN_DRIVER_FREQ = 0.005
MIN_N_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)

load(pff("BI_results.rsav"))
load(pff("blacklisted_tumor_ids.rsav"))
load(pff("distinct_tumor_ids.rsav"))
load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))


genes_to_plot = BI_results[BI_results$n_patients >= MIN_N_PATIENTS, "gene"]
BI_patient_sets = read.delim(
		paste0("DATA_USED__", this_timestamp, "/", "biallelic/Daria_2024-10-25/PPCG_Drivers/Tables/biallelic_inactivation_251024.tsv"),
		stringsAsFactors = FALSE)
BI_patient_sets = BI_patient_sets[BI_patient_sets$gene %in% genes_to_plot, ]

# exclude mets if present
load(file = pff("tracking_sheet.rsav"))
unique_donor_ids = unique(tracking_sheet$PPCG_Donor_ID)
BI_patient_sets = BI_patient_sets[BI_patient_sets$donor_id %in% unique_donor_ids,]

# keep only 1% mutated genes
gene_mut_count = by(BI_patient_sets$donor_id, BI_patient_sets$gene, function(x) length(unique(x)))
genes_to_select = names(which(gene_mut_count > MIN_N_PATIENTS))
BI_patient_sets = BI_patient_sets[BI_patient_sets$gene %in% genes_to_select,]

genes_by_muts = names(sort(by(BI_patient_sets$donor_id, BI_patient_sets$gene, function(x) length(unique(x)))))
BI_patient_sets$gene = factor(BI_patient_sets$gene, levels = genes_by_muts)
events_by_muts = names(sort(by(BI_patient_sets$sample_id, BI_patient_sets$hit_type, function(x) length(unique(x)))))
BI_patient_sets$hit_type = factor(BI_patient_sets$hit_type, levels = events_by_muts)

gene_to_pval = structure(names = BI_results$gene, BI_results$pval)
BI_patient_sets$pval_stars = gene_to_pval[BI_patient_sets$gene]
BI_patient_sets$pval_stars[duplicated(BI_patient_sets$gene)] = NA

BI_patient_sets$is_cancer_gene = NA
BI_patient_sets$is_cancer_gene[BI_patient_sets$gene %in% cgc2024] = "cancer gene"
BI_patient_sets$is_cancer_gene[BI_patient_sets$gene %in% prostate_cancer_genes] = "prostate cancer gene"
BI_patient_sets$is_cancer_gene[duplicated(BI_patient_sets$gene)] = NA

cancer_gene_symbol = c("prostate cancer gene" = 1, "cancer gene" = 4)

# three facets for genes
genes_facet1 = genes_by_muts[1:100]
genes_facet2 = genes_by_muts[101:200]
genes_facet3 = genes_by_muts[201:length(genes_by_muts)]

BI_patient_sets$gene_facet = NA
BI_patient_sets$gene_facet[BI_patient_sets$gene %in% genes_facet1] = "facet3"
BI_patient_sets$gene_facet[BI_patient_sets$gene %in% genes_facet2] = "facet2"
BI_patient_sets$gene_facet[BI_patient_sets$gene %in% genes_facet3] = "facet1"


plt = ggplot(BI_patient_sets, aes(gene, fill = hit_type)) + 
	geom_bar() +
	coord_flip() + 
	plot_theme() + 
	facet_wrap(~gene_facet, nrow = 1, scale = "free") + 
	scale_fill_manual(values = bi_colors_set1) +
	geom_point(aes(gene, y = 0, shape = is_cancer_gene)) + 
	scale_shape_manual(values = cancer_gene_symbol) + 
	geom_text(aes(label = stars.pval(pval_stars)), y = 4) + 
	theme(legend.position = "right", axis.text.y = element_text(size = 5)) +
	ggtitle("genes with biallelic inactivation", 
			paste0("N=", length(unique(BI_patient_sets$gene)), "; P<", PVALUE_CUTOFF, "; n_pat>=", ceiling(MIN_N_PATIENTS)))
	

fname = pff(c("figures/BI_gene_list_full__",  length(genes_to_plot), ".pdf"))
ggsave(plt, file = fname, width = 8.5, height = 8.5)
file_open_call2(fname)