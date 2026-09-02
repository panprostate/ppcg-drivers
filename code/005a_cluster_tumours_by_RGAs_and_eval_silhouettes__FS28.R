source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(gplots)
library(proxy)
library(RColorBrewer)
library(cluster)
library(ggsignif)
library(patchwork)
library(vegan)
library(ConsensusClusterPlus)


load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
MIN_N_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)

# driver lists and PGA from clinical data
load(pff("patient_clinical_data.rsav"))
load(pff("patient_sets_for_drivers.rsav"))
patient_sets_for_drivers = patient_sets_for_drivers[sapply(patient_sets_for_drivers, length) >= MIN_N_PATIENTS]
patients_covered = unique(unlist(patient_sets_for_drivers))

drivers_stacked = stack(patient_sets_for_drivers)
colnames(drivers_stacked) = c("patient", "driver")
drivers_stacked$mut_type = gsub("(.+)__(.+)", "\\2", drivers_stacked$driver)

# patient/rga binary matrix
# if multiple drivers per gene/patient found, collapse to one
driver_matrix = dcast(patient ~ driver, data = drivers_stacked, fun.aggregate = length, value.var = "mut_type")
driver_matrix[,-1][driver_matrix[,-1] > 1] = 1
driver_matrix_for_clustering = as.matrix(driver_matrix[,-1])
driver_matrix_for_clustering = t(apply(driver_matrix_for_clustering, 1, as.numeric))
rownames(driver_matrix_for_clustering) = driver_matrix[,1]
colnames(driver_matrix_for_clustering) = colnames(driver_matrix)[-1]
save(driver_matrix_for_clustering, file = pff("driver_matrix_for_clustering.rsav"))


# perform consensus clustering to evaluate clustering in bootstrap, avoiding ties
plots_fname_tag = pff("figures/CCP_wardD2_euclidean_10000reps_pItem08_pFeature08_maxk10")
driver_mat_for_cpp = t(driver_matrix_for_clustering)
ccp_euc3 = ConsensusClusterPlus(
  d = driver_mat_for_cpp,  # samples in columns, features in rows
  maxK = 10,
  reps = 10000,
  pItem = 0.8,
  pFeature = 0.8,
  clusterAlg = "hc",
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2",
  distance = "euclidean",
  seed = 123,
  title = plots_fname_tag,
  plot = "pdf",
  writeTable = FALSE,
  verbose = TRUE
)
ccp_euc3_icl = calcICL(ccp_euc3)
save(ccp_euc3, file = pff("ccp_euc3.rsav"))
save(ccp_euc3_icl, file = pff("ccp_euc3_icl.rsav"))

# calculate cluster silhouettes
dist_funcs = list(euclidean.dist = function(x) dist(x, method = "euclidean"))
clust_funcs = list(wardD2 = function(x) hclust(x, method = "ward.D2"))
clust_func_name = gsub("(.+)__(.+)", "\\1", SELECT_CLUSTER_METHOD)
dist_func_name = gsub("(.+)__(.+)", "\\2", SELECT_CLUSTER_METHOD)
distfun = dist_funcs[[dist_func_name]]
clustfun = clust_funcs[[clust_func_name]]


load(pff("ccp_euc3.rsav"))
load(pff("driver_matrix_for_clustering.rsav"))
load(pff("ccp_euc3_icl.rsav"))

all_silhouettes = NULL
n_clusters_span = length(ccp_euc3)
for (k in 2:n_clusters_span) {
	cat(k, " ")
	dist_mat_here = distfun(driver_matrix_for_clustering)
	cons_clusters_here = ccp_euc3[[k]]$consensusClass
	silhouette_scores = silhouette(cons_clusters_here, dist_mat_here)
	silhouette_here = data.frame(cl = SELECT_CLUSTER_METHOD, k = paste0("k", k), silhouette_scores, stringsAsFactors = FALSE)
	all_silhouettes = rbind(all_silhouettes, silhouette_here)
}			
			

# visualise and stat.compare cluster sihlouette widths
compare_grps = list(
	c("k3", "k4"),
	c("k4", "k5"),
	c("k5", "k6"),
	c("k6", "k7"),
	c("k7", "k8"),
	c("k8", "k9"),
	c("k9", "k10")
)

all_silhouettes$k = factor(all_silhouettes$k, levels = unique(all_silhouettes$k))
all_silhouettes$mean_val = ave(all_silhouettes$sil_width, interaction(all_silhouettes$cl, all_silhouettes$k), FUN = mean)
all_silhouettes$mean_val_label = signif(all_silhouettes$mean_val, 2)
all_silhouettes$mean_val_label[duplicated(as.character(interaction(all_silhouettes$cl, all_silhouettes$k)))] = NA

plt2 = ggplot(all_silhouettes, aes(k, sil_width, fill = mean_val, label = mean_val_label)) + 
		geom_boxplot(notch = T) + 
		plot_theme() + 
		facet_wrap(~cl, nrow = 1) + 
		scale_fill_gradient2(low = "steelblue", high = "darkred") + 
		geom_text(size = 2, y= min(all_silhouettes$sil_width, na.rm = TRUE)) + 
		geom_signif(aes(y = sil_width), test = "wilcox.test", test.args = list(alternative = "greater"),
				comparisons = compare_grps, map_signif_level = TRUE, textsize = 2)

fname = pff("cluster_silhouettes_scores_CPP.pdf")
ggsave(plt2, file = fname, width = 5, height = 5)
file_open_call2(fname)


# visualise consensus scores of clustering
consensus_scores = ccp_euc3_icl$itemConsensus
consensus_scores$k = factor(paste0("k", consensus_scores$k), levels = paste0("k", sort(unique(consensus_scores$k))))
consensus_scores$cl = SELECT_CLUSTER_METHOD
consensus_scores$mean_val = ave(consensus_scores$itemConsensus, interaction(consensus_scores$cl, consensus_scores$k), FUN = mean)
consensus_scores$mean_val_label = signif(consensus_scores$mean_val, 2)
consensus_scores$mean_val_label[duplicated(as.character(interaction(consensus_scores$cl, consensus_scores$k)))] = NA

plt2 = ggplot(consensus_scores, aes(k, itemConsensus, fill = mean_val, label = mean_val_label)) + 
		geom_boxplot(notch = T) + 
		plot_theme() + 
		facet_wrap(~cl, nrow = 1) + 
		scale_fill_gradient2(low = "steelblue", high = "darkred") + 
		geom_text(size = 2, y= min(consensus_scores$itemConsensus, na.rm = TRUE)) + 
		geom_signif(aes(y = itemConsensus), test = "wilcox.test", test.args = list(alternative = "greater"),
				comparisons = compare_grps, map_signif_level = TRUE, textsize = 2)

fname = pff("cluster_itemConsensus_scores_CPP.pdf")
ggsave(plt2, file = fname, width = 5, height = 5)
file_open_call2(fname)


# review the data and decide the optimal number of clusters
N_CLUSTERS = 4
members_selected_clustering = ccp_euc3[[N_CLUSTERS]]$consensusClass

patient_cluster_membership = data.frame(names(members_selected_clustering), members_selected_clustering, stringsAsFactors = FALSE)
patient_cluster_membership[,2] = paste0("cluster_", patient_cluster_membership[,2])
colnames(patient_cluster_membership) = c("patient", "group_id")
rownames(patient_cluster_membership) = NULL
save(patient_cluster_membership, file = pff("patient_cluster_membership.rsav"))

patient_cluster_certainty = ccp_euc3_icl$itemConsensus[ccp_euc3_icl$itemConsensus$k == N_CLUSTERS,]
patient_cluster_certainty = c(by(patient_cluster_certainty$itemConsensus, patient_cluster_certainty$item, max))
save(patient_cluster_certainty, file = pff("patient_cluster_certainty.rsav"))
