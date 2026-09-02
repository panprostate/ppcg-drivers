source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(gtools)
library(patchwork)

load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
MIN_N_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)

# merge drivers across alteration types; keep 1%+ frequent
load(file = pff("patient_sets_for_drivers.rsav"))
stacked_rgas = data.frame(as.matrix(stack(patient_sets_for_drivers)), stringsAsFactors = FALSE)
stacked_rgas$mut_type = gsub("(.+)__(.+)", "\\2", stacked_rgas$ind)
stacked_rgas$gene = gsub("(.+)__(.+)", "\\1", stacked_rgas$ind)

collapsed_rgas_list = by(stacked_rgas$values, stacked_rgas$gene, unique)
collapsed_rgas_list = collapsed_rgas_list[sapply(collapsed_rgas_list, length) > MIN_N_PATIENTS]
stacked_rgas2 = data.frame(as.matrix(stack(collapsed_rgas_list)), stringsAsFactors = FALSE) 

# genes x patients
rga_matrix = dcast(stacked_rgas2, ind ~ values, value.var = "ind", fun.aggregate = length)
rownames(rga_matrix) = rga_matrix[,1]
rga_matrix_memosort = memoSort( rga_matrix[,-1])

patient_order = colnames(rga_matrix_memosort)
gene_order = rev(rownames(rga_matrix_memosort))

colnames(stacked_rgas2) = c("patient", "gene")
stacked_rgas2$is_mutated = 1

stacked_rgas2$patient = factor(stacked_rgas2$patient, levels = patient_order)
stacked_rgas2$gene = factor(stacked_rgas2$gene, levels = gene_order)

n_pat = length(unique(stacked_rgas2$patient))
ggtitl = paste0("all RGAs (n=", length(gene_order),"; min ", 100*MIN_DRIVER_FREQ, "%); ", length(patient_order), " tumors" )
plt_oncoprint = ggplot(stacked_rgas2, aes(patient, gene, fill = is_mutated)) + 
		geom_tile() + 
		scale_x_discrete("patients", breaks = NULL, labels = NULL) +
		plot_theme() +
		theme(axis.text.y = element_text(size = 6, color = "black")) + 
		ggtitle(NULL, ggtitl)
		
fname = pff("figures/RGA_oncoprint_for_ALL.pdf")
ggsave(plt_oncoprint, file = fname)
file_open_call2(fname)






#
# oncoprint for ETS- / SPOP- tumors ; ask which RGAs are enriched in this subset
#
rga_matrix_ESneg = rga_matrix

ETS_neg_samples = colnames(rga_matrix_ESneg)[which(rga_matrix_ESneg["ETS",] == 0)]
SPOP_neg_samples = colnames(rga_matrix_ESneg)[which(rga_matrix_ESneg["SPOP",] == 0)]
ETS_SPOP_neg_samples = intersect(ETS_neg_samples, SPOP_neg_samples)
MIN_N_PATIENTS_ESneg = MIN_DRIVER_FREQ * length(ETS_SPOP_neg_samples)


# include only ETS-neg SPOP-neg samples
# keep only the genes with min 1% freq in this subset of samples
rga_matrix_ESneg = rga_matrix_ESneg[,colnames(rga_matrix_ESneg) %in% c("ind", ETS_SPOP_neg_samples)]
select_RGAs_ESneg = names(which(apply(rga_matrix_ESneg[,-1], 1, sum) > MIN_N_PATIENTS_ESneg))
rga_matrix_ESneg = rga_matrix_ESneg[rga_matrix_ESneg$ind %in% select_RGAs_ESneg, ]
rga_matrix_ESneg_memosort = memoSort( rga_matrix_ESneg[,-1])

patient_order = colnames(rga_matrix_ESneg_memosort)
gene_order = rev(rownames(rga_matrix_ESneg_memosort))

stacked_rgas2_ESneg = stacked_rgas2
colnames(stacked_rgas2_ESneg) = c("patient", "gene", "is_mutated")
stacked_rgas2_ESneg$gene = as.character(stacked_rgas2_ESneg$gene)

stacked_rgas2_ESneg = stacked_rgas2_ESneg[stacked_rgas2_ESneg$patient %in% patient_order,]
stacked_rgas2_ESneg = stacked_rgas2_ESneg[stacked_rgas2_ESneg$gene %in% gene_order,]
stacked_rgas2_ESneg$patient = factor(stacked_rgas2_ESneg$patient, levels = patient_order)
stacked_rgas2_ESneg$gene = factor(stacked_rgas2_ESneg$gene, levels = gene_order)
n_pat = length(unique(stacked_rgas2_ESneg$patient))

# stat analysis: which drivers are especially enriched in this SPOP-neg / ETS-neg set
ETSPOP_minus_patients = unique(stacked_rgas2_ESneg$patient)


test_ETSPOP_RGAs = function(gene, ETSPOP_minus_patients, all_patients, collapsed_rgas_list) {
	
	patients_with_rga = collapsed_rgas_list[[gene]]
	ft = fisher.test(all_patients %in% ETSPOP_minus_patients, all_patients %in% patients_with_rga)
	pval = ft$p.value
	or = ft$estimate
	n_intersect = length(intersect(ETSPOP_minus_patients, patients_with_rga))
	n_rga = length(patients_with_rga)
	
	pct_etspop_pos = length(setdiff(patients_with_rga, ETSPOP_minus_patients)) / length(setdiff(all_patients, ETSPOP_minus_patients))
	pct_etspop_neg = n_intersect / length(ETSPOP_minus_patients)
	obs_exp_FC = pct_etspop_neg / pct_etspop_pos

	data.frame(gene, n_rga, n_intersect, pct_etspop_pos, pct_etspop_neg, obs_exp_FC, or, pval, stringsAsFactors = FALSE)
}

genes_to_check = unique(as.character(stacked_rgas2_ESneg$gene))
ETSPOP_minus_stats = do.call(rbind, lapply(genes_to_check, test_ETSPOP_RGAs, ETSPOP_minus_patients, all_patients, collapsed_rgas_list))
ETSPOP_minus_stats$fdr = p.adjust(ETSPOP_minus_stats$pval, method = "fdr")
save(ETSPOP_minus_stats, file = pff("ETSPOP_minus_stats.rsav"))


ETSPOP_minus_stats$gene = factor(ETSPOP_minus_stats$gene, levels = levels(stacked_rgas2_ESneg$gene))

LOG2_ODDS_CAP = 2
ETSPOP_minus_stats$or_log2_cap = log2(ETSPOP_minus_stats$or)
ETSPOP_minus_stats$or_log2_cap[ETSPOP_minus_stats$or_log2_cap > LOG2_ODDS_CAP] = LOG2_ODDS_CAP
ETSPOP_minus_stats$or_log2_cap[ETSPOP_minus_stats$or_log2_cap < -LOG2_ODDS_CAP] = -LOG2_ODDS_CAP

ggtitl = paste0("ETS-/SPOP- RGAs (n=", length(gene_order),"; min ", 100*MIN_DRIVER_FREQ, "%); ", length(patient_order), " tumors" )

plt_oncoprint = ggplot(stacked_rgas2_ESneg, aes(patient, gene, fill = is_mutated)) + 
		geom_tile() + 
		scale_x_discrete("patients", breaks = NULL, labels = NULL) +
		plot_theme() +
		theme(axis.text.y = element_text(size = 6, color = "black"), legend.position = "bottom") + 
		ggtitle(NULL, ggtitl)
		
plt_stats = ggplot(ETSPOP_minus_stats, aes(factor(1), gene, fill = or_log2_cap, label = stars.pval(fdr))) + 
		geom_tile() + 
		geom_text() + 
		plot_theme() + 
		scale_fill_gradient2(mid = "white", high = "darkred", low = "darkblue") + 
		theme(axis.text.y = element_text(size = 6, color = "black"), legend.position = "bottom")


plt_combined = (plt_stats + plt_oncoprint) + plot_layout(widths = c(1, 8))

fname = pff("figures/RGA_oncoprint_for_ETSneg_SPOPneg.pdf")
ggsave(plt_combined, file = fname, height = 9, width = 9)
file_open_call2(fname)



#
# oncoprint for ETS factors and SLC30A4
#
MIN_FUSIONS = 2

# take first or second partner from "G1:G2"
get_fusion_partners = function(i, dat_ETS, which_to_take) {

	fusions = unique(gsub("(.+):(.+)", paste0("\\", which_to_take), setdiff(unlist(strsplit(unlist(dat_ETS[i, -1:-9]), s = ";|,")), "FALSE")))
	fusions
}

# tgt_gene = "SLC30A4"; fusions_gene = "TMPRSS2"
test_mut_exclusivity = function(fusions_gene, tgt_gene, fusions_matrix, all_patients) {
	tgt_patients = names(which(sapply(fusions_matrix[tgt_gene, -1, drop = TRUE], function(x) any(x>0))))
	fusions_patients = names(which(sapply(fusions_matrix[fusions_gene, -1, drop = TRUE], function(x) any(x>0))))
	ft = fisher.test(all_patients %in% tgt_patients, all_patients %in% fusions_patients, alt = "l")
	
	data.frame(tgt_gene, fusions_gene, pval = ft$p.value, or = ft$estimate, stringsAsFactors = FALSE)
}


# ets fusions are encoded as PARTNER:ETS_GENE; splitting these will be element 1 and element 2

for (i in 1:2) {
	
	gene_type = c("ETS_partners", "ETS_genes")[i]
	
	# combine ETS fusion partners into a stack
	# get either fusion partner or ETS gene 
	load(file = pff("dat_ETS.rsav"))
	fusions_list = lapply(1:nrow(dat_ETS), get_fusion_partners, dat_ETS, i)
	names(fusions_list) = dat_ETS$PPCG_Donor_ID
	
	# gene,patient pairs
	# remove sparse pairs
	fusions_stack = as.matrix(stack(fusions_list))
	nonsparse_genes = names(which(table(fusions_stack[,1]) >= MIN_FUSIONS))
	fusions_stack = fusions_stack[fusions_stack[,1] %in% nonsparse_genes,]
	
	# add all SLC SV patients to the stack
	SLC_RGA_patients = patient_sets_for_drivers[["SLC30A4__SV"]]
	fusions_stack_plus_SLC = data.frame(unique(rbind(fusions_stack, cbind("SLC30A4", SLC_RGA_patients))), stringsAsFactors = FALSE)  
	
	# fusions x patients
	fusions_matrix = dcast(fusions_stack_plus_SLC, values ~ ind, value.var = "ind", fun.aggregate = length)
	rownames(fusions_matrix) = fusions_matrix[,1]
	fusions_matrix_memosort = memoSort( fusions_matrix[,-1])
	patient_order = colnames(fusions_matrix_memosort)
	gene_order = rev(rownames(fusions_matrix_memosort))
	
	colnames(fusions_stack_plus_SLC) = c("gene", "patient")
	fusions_stack_plus_SLC$has_fusion = 1
	fusions_stack_plus_SLC$patient = factor(fusions_stack_plus_SLC$patient, levels = patient_order)
	fusions_stack_plus_SLC$gene = factor(fusions_stack_plus_SLC$gene, levels = gene_order)
	
	SLC_fusion_stats = do.call(rbind, lapply(as.character(unique(fusions_stack_plus_SLC$gene)), 
			test_mut_exclusivity, "SLC30A4", fusions_matrix, all_patients))
	SLC_fusion_stats$fdr = p.adjust(SLC_fusion_stats$pval, method = "fdr")
	SLC_fusion_stats$fusions_gene = factor(SLC_fusion_stats$fusions_gene, levels = gene_order)	
	
	n_pat = length(unique(fusions_stack_plus_SLC$patient))
	ggtitl = paste0(gene_type, " (n=", length(gene_order),");  ", length(patient_order), " tumours; min_fusions=", MIN_FUSIONS)
	
	plt_oncoprint = ggplot(fusions_stack_plus_SLC, aes(patient, gene, fill = has_fusion)) + 
			geom_tile() + 
			scale_x_discrete("patients", breaks = NULL, labels = NULL) +
			plot_theme() +
			theme(axis.text.y = element_text(size = 6, color = "black")) + 
			ggtitle(NULL, ggtitl)
	
	plt_stats = ggplot(SLC_fusion_stats, aes(factor(1), fusions_gene, label = stars.pval(fdr))) + 
			geom_text() + 
			plot_theme() +
			theme(axis.text.y = element_text(size = 1, color = "black"))
	
	plt_combined = plt_oncoprint + plt_stats + plot_layout(widths = c(3, 1))
			
	fname = pff(c("figures/SLC30A4_ERGfusion_oncoprint_", gene_type, ".pdf"))
	ggsave(plt_combined, file = fname, width = 12, height = 3)
	file_open_call2(fname)
}	
