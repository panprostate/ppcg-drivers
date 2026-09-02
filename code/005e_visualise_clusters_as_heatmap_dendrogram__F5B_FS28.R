source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(gplots)
library(proxy)
library(RColorBrewer)
library(patchwork)
library(ggsignif)
library(cluster)


load(file = pff("ccp_euc3.rsav"))
load(file = pff("driver_matrix_for_clustering.rsav"))

dist_funcs = list(euclidean.dist = function(x) dist(x, method = "euclidean"))
clust_funcs = list(wardD2 = function(x) hclust(x, method = "ward.D2"))
clust_func_name = gsub("(.+)__(.+)", "\\1", SELECT_CLUSTER_METHOD)
dist_func_name = gsub("(.+)__(.+)", "\\2", SELECT_CLUSTER_METHOD)
distfun = dist_funcs[[dist_func_name]]
clustfun = clust_funcs[[clust_func_name]]

# based on the consensus clustering
cons_tree = ccp_euc3[[SELECT_N_CLUSTERS]]$consensusTree
this_dendrogram = as.dendrogram(cons_tree)

# colors on the left - tumor/cluster annotations
load(pff("patient_cluster_membership.rsav"))
cluster_no = structure(as.numeric(gsub("cluster_", "", patient_cluster_membership[,2])), names = patient_cluster_membership[,1])
cluster_no = cluster_no[rownames(driver_matrix_for_clustering)]
cluster_colors_side = CLUSTER_COLORS[cluster_no]

# colors at the top - RGA alteration types
alt_type_color_legend = c(
	"hAMP" = "darkred",
	"BI" = "darkblue",
	"SNV_CDS" = "gold",
	"SNV_NC" = "darkorange",
	"SV" = "darkgreen",
	"gain" = "brown2",
	"loss" = "steelblue"
)
alt_types = gsub("(.+)__(.+)", "\\2", colnames(driver_matrix_for_clustering))
alteration_colors_side = alt_type_color_legend[alt_types]



fname_col = pff(c("CCP_heatmap_for_clusters__", SELECT_CLUSTER_METHOD, "_k", SELECT_N_CLUSTERS, ".pdf"))
pdf(fname_col, width = 100, height = 100)

heatmap_colored = heatmap.2(
		driver_matrix_for_clustering,
		distfun = distfun,
		hclustfun = clustfun, 
		Rowv = this_dendrogram,
		col = c("white", "darkcyan"), 		
		trace = "none",
		key = FALSE,
		main = paste0(SELECT_CLUSTER_METHOD, "(CCP); k=", SELECT_N_CLUSTERS),
		RowSideColors = cluster_colors_side, 
		ColSideColors = alteration_colors_side, 
		margins = c(16,8)
)
		
legend("left", title = "Clusters", legend = paste0("cl", 1:SELECT_N_CLUSTERS), 
       fill = CLUSTER_COLORS[1:SELECT_N_CLUSTERS], cex=0.8, box.lty=0)
legend("top", title = "Alterations", legend = names(alt_type_color_legend), 
       fill = alt_type_color_legend, cex=0.8, box.lty=0)

dev.off()
file_open_call2(fname_col)
save(heatmap_colored, file = pff(c("heatmap_colored__", SELECT_CLUSTER_METHOD, "_k", SELECT_N_CLUSTERS, "_CCP.rsav")))


# add sidebars using clinical info
load(pff("patient_clinical_data.rsav")) 
patients_in_heatmap = colnames(heatmap_colored$carpet)
patient_clinical_data = patient_clinical_data[patient_clinical_data$patient %in% patients_in_heatmap,]
patient_with_cluster_no = data.frame(patient = names(cluster_no), cluster_no = cluster_no, stringsAsFactors = FALSE)
patient_clinical_data = merge(patient_clinical_data, patient_with_cluster_no, by = "patient")
patient_clinical_data$patient = factor(patient_clinical_data$patient, levels = patients_in_heatmap)

plt_cluster_no = ggplot(patient_clinical_data, aes(factor(1), patient, fill = factor(cluster_no))) + 
		geom_tile() + 
		plot_theme() + 
		scale_fill_manual("k", values = CLUSTER_COLORS[sort(unique(cluster_no))]) + 
		theme(axis.text.y = element_text(size = 2, color = "black"))

plt_gleason = ggplot(patient_clinical_data, aes(factor(1), patient, fill = factor(Gleason_group))) + 
		geom_tile() + 
		plot_theme() + 
		scale_fill_brewer("grade", palette = "OrRd", na.value = "grey") + 
		theme(axis.text.y = element_text(size = 2, color = "black"))
		
plt_stage = ggplot(patient_clinical_data, aes(factor(1), patient, fill = factor(T_stage))) + 
		geom_tile() + 
		plot_theme() + 
		scale_fill_brewer("grade", palette = "PuRd", na.value = "grey") + 
		theme(axis.text.y = element_text(size = 2, color = "black"))

plt_PSA = ggplot(patient_clinical_data, aes(factor(1), patient, fill = PSA_log2)) + 
		geom_tile() + 
		plot_theme() + 
		scale_fill_distiller("PSA", na.value = "grey", palette = "Reds", direction = 1) + 
		theme(axis.text.y = element_text(size = 2, color = "black"))
		
plt_metsbio = ggplot(patient_clinical_data, aes(factor(1), patient, fill = new_metastatic_biology_indicator)) + 
		geom_tile() + 
		plot_theme() + 
		scale_fill_manual("metsbio", values = c("mets_biol" = "darkred", "no_mets_biol" = "white"), na.value = "grey") + 
		theme(axis.text.y = element_text(size = 2, color = "black"))

plt_eo = ggplot(patient_clinical_data, aes(factor(1), patient, fill = early_onset)) + 
		geom_tile() + 
		plot_theme() + 
		scale_fill_manual("EO", values = c("TRUE" = "cornflowerblue", "FALSE" = "white"), na.value = "grey") + 
		theme(axis.text.y = element_text(size = 2, color = "black"))

plt_pga = ggplot(patient_clinical_data, aes(factor(1), patient, fill = PGA)) + 
		geom_tile() + 
		plot_theme() + 
		scale_fill_distiller("PGA", na.value = "grey", palette = "Blues", direction = 1) + 
		theme(axis.text.y = element_text(size = 2, color = "black"))

plt_tmb = ggplot(patient_clinical_data, aes(factor(1), patient, fill = GTMB_log1p)) + 
		geom_tile() + 
		plot_theme() + 
		scale_fill_distiller("TMB", na.value = "grey", palette = "Blues", direction = 1) + 
		theme(axis.text.y = element_text(size = 2, color = "black"))

plt_combined = plt_cluster_no | plt_gleason | plt_stage | plt_metsbio | plt_PSA | plt_pga | plt_eo | plt_tmb

fname = pff(c("figures/heatmap_sidebars__", SELECT_CLUSTER_METHOD, "_k", SELECT_N_CLUSTERS, ".pdf"))
ggsave(plt_combined, file = fname, height = 15, width = 16)
file_open_call2(fname)

