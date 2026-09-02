source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(gplots)
library(ActivePathways)
library(survival)
library(survminer)
library(patchwork)
library(gtools)

load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
known_cancer_genes = c(cgc2024, prostate_cancer_genes)


load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
MIN_MUT_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)
FDR_CUTOFF = 0.05
FDR_SCORE_CAP = 16

load(pff("patient_cluster_membership.rsav"))
load(pff("patient_sets_for_drivers.rsav"))
patient_sets_for_drivers = patient_sets_for_drivers[sapply(patient_sets_for_drivers, length) >= MIN_MUT_PATIENTS]

all_RGAs = names(patient_sets_for_drivers)
all_clusters = unique(patient_cluster_membership$group_id)
todo_list = as.matrix(expand.grid(all_RGAs, all_clusters))

# for each RGA test if it is enriched in a given cluster compared to other clusters combined
test_rga_cluster_enrichment = function(i, todo_list, all_patients, patient_cluster_membership, patient_sets_for_drivers) {

	this_rga = todo_list[i, 1]
	this_cluster = todo_list[i, 2]
	
	patients_in_cluster = patient_cluster_membership[patient_cluster_membership$group_id == this_cluster, "patient"]
	patients_with_rga = patient_sets_for_drivers[[this_rga]]
	
	fish_test = fisher.test(
			factor(all_patients %in% patients_in_cluster, levels = c(TRUE, FALSE)), 
			factor(all_patients %in% patients_with_rga, levels = c(TRUE, FALSE)), 
			alt = "greater")
	pval = fish_test$p.value
	or = fish_test$estimate
	
	n_in_cluster = length(patients_in_cluster)
	n_with_rga = length(patients_with_rga)
	shared_patients = intersect(patients_in_cluster, patients_with_rga)
	n_overlap = length(shared_patients)
	
	fract_of_cluster_patients = n_overlap / n_in_cluster
	fract_of_rga_patients = n_overlap / n_with_rga
	
	data.frame(this_rga, this_cluster, pval, or, n_in_cluster, n_with_rga, n_overlap, fract_of_cluster_patients, fract_of_rga_patients, 
			stringsAsFactors = FALSE)
}


rgas_enriched_in_clusters = do.call(rbind, lapply(1:nrow(todo_list), 
		test_rga_cluster_enrichment, todo_list, all_patients, patient_cluster_membership, patient_sets_for_drivers))
rgas_enriched_in_clusters$fdr = p.adjust(rgas_enriched_in_clusters$pval, method = "fdr")


load(pff("gr_gene_coords.rsav"))
# add cytoband labels to genes
cytobands = read.delim("elements/cytoBand.txt", stringsAsFactors = FALSE, header = FALSE)
cytobands$id = paste0(cytobands$V1, cytobands$V4)
genes_here = unique(gsub("(.+)__(.+)", "\\1", rgas_enriched_in_clusters$this_rga))
gr_genes_here = gr_gene_coords[gr_gene_coords$symbol %in% genes_here]
gene2cytoband = get_gene_cytobands(gr_genes_here, cytobands)
gene2cytoband = structure(names = gene2cytoband$gene, gene2cytoband$cytoband)
rgas_enriched_in_clusters$cytoband = gene2cytoband[gsub("(.+)__(.+)", "\\1", rgas_enriched_in_clusters$this_rga)]
rgas_enriched_in_clusters$cytoband[grep("__gain|__loss|__hAMP", rgas_enriched_in_clusters$this_rga, invert = TRUE)] = ""
rgas_enriched_in_clusters$this_rga2 = paste0(gsub("__", " ", rgas_enriched_in_clusters$this_rga), " (", rgas_enriched_in_clusters$cytoband, ")")
rgas_enriched_in_clusters$this_rga2 = gsub(" \\(\\)", "", rgas_enriched_in_clusters$this_rga2)

save(rgas_enriched_in_clusters, file = pff("rgas_enriched_in_clusters.rsav"))


# visualise RGA/cluster associations as piecharts
# select RGAs that are significant in at least one cluster, plot their signals in all clusters
rgas_with_any_enrichments = names(which(by(rgas_enriched_in_clusters$fdr, rgas_enriched_in_clusters$this_rga2, min) < FDR_CUTOFF))
rgas_enriched_in_clusters = rgas_enriched_in_clusters[rgas_enriched_in_clusters$this_rga2 %in% rgas_with_any_enrichments,]
rgas_enriched_in_clusters = 
		rgas_enriched_in_clusters[order(rgas_enriched_in_clusters$this_cluster, -log10(rgas_enriched_in_clusters$fdr), decreasing = T),]
rgas_ordered_by_cluster_fraction = as.character(unique(rgas_enriched_in_clusters$this_rga2))
rgas_enriched_in_clusters$this_rga2 = factor(rgas_enriched_in_clusters$this_rga2, levels = rgas_ordered_by_cluster_fraction)
rgas_enriched_in_clusters$this_cluster2 = gsub("cluster_", "c", rgas_enriched_in_clusters$this_cluster)
rgas_enriched_in_clusters$this_cluster2 = factor(rgas_enriched_in_clusters$this_cluster2, 
		levels = rev(sort(unique(rgas_enriched_in_clusters$this_cluster2))))

# capping, colors
rgas_enriched_in_clusters$enrichment_score = -log10(rgas_enriched_in_clusters$fdr) * sign(log2(rgas_enriched_in_clusters$or))
rgas_enriched_in_clusters$enrichment_score = pmax(rgas_enriched_in_clusters$enrichment_score, -FDR_SCORE_CAP)
rgas_enriched_in_clusters$enrichment_score = pmin(rgas_enriched_in_clusters$enrichment_score, +FDR_SCORE_CAP)
rgas_enriched_in_clusters$is_significant = c("NO", "YES")[c(rgas_enriched_in_clusters$fdr < FDR_CUTOFF) + 1]


plt_title = paste0("clusters: ", SELECT_CLUSTER_METHOD, "__", SELECT_N_CLUSTERS, "; FDR=", FDR_CUTOFF)
plt = ggplot(rgas_enriched_in_clusters, 
		aes(factor(1), fract_of_cluster_patients, fill = enrichment_score, label = stars.pval(fdr), color = is_significant)) +
		geom_bar(stat = "identity") +
		scale_fill_gradient2(low = "darkcyan", high = "darkred", mid = "white") + 
		scale_color_manual(values = c("YES" = "black", "NO" = "lightgrey")) + 
		scale_y_continuous(NULL, breaks = NULL, labels = NULL, lim = c(0, 1)) +
		scale_x_discrete(NULL, breaks = NULL, labels = NULL) +
		coord_polar(theta = "y") + 
		plot_theme() + 
		theme(strip.background = element_blank(), 
				strip.text.x = element_text(angle = 90), panel.spacing = unit(0.01, "lines"), 
				legend.position = "bottom") +
		facet_grid(this_cluster2 ~ this_rga2, switch = "both") + 
		ggtitle(NULL, plt_title)

fname = pff(c("figures/cluster_genetics_fractions__", SELECT_CLUSTER_METHOD, "__", SELECT_N_CLUSTERS, ".pdf"))
ggsave(plt, file = fname, width = 8.5, height = 5)
file_open_call2(fname)

cluster_genetics_select_for_piecharts = rgas_enriched_in_clusters
save(cluster_genetics_select_for_piecharts, file = pff("cluster_genetics_select_for_piecharts.rsav"))
