source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(ggrepel)
library(scales)
library(patchwork)
library(RColorBrewer)

MIN_COHORT_FREQUENCY = 0.01
load(file = pff("all_patients.rsav"))
MIN_N_PATIENTS  = MIN_COHORT_FREQUENCY * length(all_patients)

get_gene_BI_or_AMP = function(gene, full_results, mut) {
	results_here = full_results[full_results$gene == gene,, drop = FALSE]
	if (nrow(results_here) == 0) {
		return(NULL)
	}
	patient = unique(unlist(strsplit(results_here$patients, split = ",")))
	data.frame(gene, patient, mut, event_type = mut, stringsAsFactors = FALSE)
}

get_gene_CDS_SNVs = function(gene, dnds_results_annotations_RGA) {
	
	trunc_impacts = c("Nonsense", "indel_delfrshift", "Essential_Splice", "indel_insfrshift", "Stop_loss")
	subst_impacts = c("Missense", "indel_delinframe")
	
	this_cds_variants = dnds_results_annotations_RGA[dnds_results_annotations_RGA$gene == gene,]
	if (nrow(this_cds_variants) == 0) {
		return(NULL)	
	}
	
	patient_trunc = unique(this_cds_variants[this_cds_variants$impact %in% trunc_impacts, "sampleID"])
	patient_subst = unique(this_cds_variants[this_cds_variants$impact %in% subst_impacts, "sampleID"])

	events = list(
			"truncating" = patient_trunc, 
			"substitution" = patient_subst
	)
	events = as.matrix(stack(events)) 
	colnames(events) = c("patient", "mut")
	data.frame(gene = gene, events, event_type = "SNV_indel", stringsAsFactors = FALSE)
}

get_gene_NC_SNVs = function(gene, results_NC_SNV_drivers, var2el_NC_SNV) {

	# element has to be considered a driver for this analysis
	this_elements = strsplit(results_NC_SNV_drivers[results_NC_SNV_drivers$annots_MAIN == gene, "element_ids"], split = ",")
	this_elements = unique(unlist(this_elements))
	if (length(this_elements) == 0) {
		return(NULL)
	}
	
	this_muts = do.call(rbind, var2el_NC_SNV[this_elements])
	this_patients = unique(this_muts$patient)
	
	events = list("noncoding" = this_patients)
	events = as.matrix(stack(events)) 
	colnames(events) = c("patient", "mut")
	data.frame(gene = gene, events, event_type = "SNV_indel", stringsAsFactors = FALSE)
}

get_gene_SVs = function(gene, results_SV_drivers, var2el_SV, ETS_SV_status) {
	
	if (gene == "ETS") {
		dfr = data.frame(gene, patient = ETS_SV_status$patient, mut = paste0("SV", tolower(ETS_SV_status$SV_type)), event_type = "SV", 
				stringsAsFactors = FALSE)
		return(dfr)
	}
	
	# element has to be considered a driver for this analysis
	this_elements = strsplit(results_SV_drivers[results_SV_drivers$annots_MAIN == gene, "element_ids"], split = ",")
	this_elements = unique(unlist(this_elements))
	if (length(this_elements) == 0) {
		return(NULL)
	}

	this_muts = do.call(rbind, var2el_SV[this_elements])
	rownames(this_muts) = NULL
	this_muts[this_muts$top_mut_signt %in% c("h2hINV", "t2tINV"), "top_mut_signt"] = "INV"
	
	patient_events = c(by(this_muts$top_mut_signt, this_muts$patient, function(x) paste(sort(unique(x)), collapse=",")))
	patient_events[grep(",", patient_events)] = "cplx"
	patient_events = sapply(patient_events, function(x) paste0("SV", tolower(x)))

	events = split(names(patient_events), patient_events)	
	events = as.matrix(stack(events)) 
	colnames(events) = c("patient", "mut")
	
	data.frame(gene = gene, events, event_type = "SV", stringsAsFactors = FALSE)
}

get_combined_events = function(gene, BI_results, AMP_results, dnds_results_annotations_RGA,
		results_NC_SNV_drivers, var2el_NC_SNV, 
		results_SV_drivers, var2el_SV, 
		ETS_SV_status) {
	
	cat(gene, " ")
	
	BI_events = get_gene_BI_or_AMP(gene, BI_results, "BI")
	AMP_events = get_gene_BI_or_AMP(gene, AMP_results, "hAmp")
	NC_SNV_events = get_gene_NC_SNVs(gene, results_NC_SNV_drivers, var2el_NC_SNV)
	CDS_SNV_events = get_gene_CDS_SNVs(gene, dnds_results_annotations_RGA)
	SV_events = get_gene_SVs(gene, results_SV_drivers, var2el_SV, ETS_SV_status)
	
	all_events = rbind(BI_events, AMP_events, NC_SNV_events, CDS_SNV_events, SV_events)
	all_events_collapse = c(by(all_events$mut, all_events$patient, function(x) paste(sort(x), collapse = "__")))
	all_events = data.frame(gene, patient = names(all_events_collapse), mut = all_events_collapse, stringsAsFactors = FALSE)
	unique(all_events)
}


load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
prostate_cancer_genes = c(prostate_cancer_genes, "ETS")
cancer_genes = unique(c(cgc2024, prostate_cancer_genes))

# separate lists for BI
load(pff("BI_results.rsav"))
BI_results = BI_results[BI_results$is_cancer_gene,]

# separate lists for AMP
load(pff("gene_CNA_module_dfr.rsav"))
AMP_results = gene_CNA_module_dfr[gene_CNA_module_dfr$annot == "high_gain",]

# this provides protein-coding (CDS) SNV driver list and CDS SNVs associated with each driver element
load(pff("results_signf_merged_annot.rsav"))
results_CDS_SNV_drivers = results_signf_merged_annot[results_signf_merged_annot$mut_type == "SNV_CDS",]
rm(results_signf_merged_annot)

# dnds results are the basis of CDS mutations, but we need to focus on genes that were selected by DNDS
load(pff("dnds_results_annotations_all.rsav"))
dnds_results_annotations_RGA = dnds_results_annotations_all[dnds_results_annotations_all$gene %in% results_CDS_SNV_drivers$annots_MAIN,] 
rm(dnds_results_annotations_all)

# this provides non-coding (NC) SNV driver list and NC SNVs associated with each driver element
load(pff("results_signf_merged_annot.rsav"))
results_NC_SNV_drivers = results_signf_merged_annot[results_signf_merged_annot$mut_type == "SNV_NC",]

load(pff("variants_to_elements.rsav"))
var2el_NC_SNV = variants_to_elements[["SNV_indel"]]
rm(variants_to_elements, results_signf_merged_annot)
gc()



# this provides SV driver list and SVs associated with each driver element
load(pff("results_signf_merged_annot.rsav"))
results_SV_drivers = results_signf_merged_annot[results_signf_merged_annot$mut_type %in% c("SV"),]

load(pff("variants_to_elements.rsav"))
var2el_SV = variants_to_elements[["SV"]]
rm(variants_to_elements, results_signf_merged_annot)
gc()

# this provides ets-altered genes as a separate set
load(file = pff("ETS_SV_status.rsav"))


# collect all RGA genes from the alteration types
selected_genes = unique(c(
		BI_results$gene, 
		AMP_results$gene, 
		results_CDS_SNV_drivers$annots_MAIN, 
		results_NC_SNV_drivers$annots_MAIN, 
		results_SV_drivers$annots_MAIN))

gene_combined_events = do.call(rbind, lapply(selected_genes, get_combined_events, 
		BI_results, AMP_results, dnds_results_annotations_RGA,
		results_NC_SNV_drivers, var2el_NC_SNV, 
		results_SV_drivers, var2el_SV, 
		ETS_SV_status))


gene_freqs = sapply(by(gene_combined_events$patient, gene_combined_events$gene, unique), length)
genes_to_select = names(gene_freqs[gene_freqs >= MIN_N_PATIENTS])

gene_combined_events = gene_combined_events[gene_combined_events$gene %in% genes_to_select,]

events_by_table = rev(sort(table(gene_combined_events$mut)))
events_by_table = data.frame(mut = names(events_by_table), freq = c(events_by_table), ANNOT = names(events_by_table), stringsAsFactors = FALSE)
rownames(events_by_table) = NULL

fname = pff(c("figures/006b_events_by_table.csv"))
write.csv(events_by_table, file = fname)
file_open_call2(fname)
# curate manually 
fname_new = pff(c("figures/006b_events_by_table_EDIT.csv"))
events_by_table_EDIT = read.csv(file = fname_new, stringsAsFactors = FALSE)
events_by_table_EDIT = events_by_table_EDIT[,-1]


# add the manually curated event types to the matrix
gene_combined_events_EDIT = merge(gene_combined_events, events_by_table_EDIT, by = "mut", all = TRUE)
gene_combined_events_EDIT = gene_combined_events_EDIT[gene_combined_events_EDIT$gene %in% genes_to_select, ]

# order genes by total cohort coverage
genes_by_freq = names(sort(table(gene_combined_events_EDIT$gene)))
gene_combined_events_EDIT$gene = factor(gene_combined_events_EDIT$gene, levels = genes_by_freq)

mut_colors = c(
		"BI"      = "darkblue",
		"SVdel"     = "royalblue1" ,
		"SVcplx"    = "darkolivegreen1" ,
		"SVtra"     = "darkolivegreen" ,
		"noncoding" = "darkorange" ,	
		"substitution"     = "gold" ,
		"SVinv"     = "darkkhaki" ,
		"hAmp"      = "darkred",
		"truncating"  = "azure4" ,
		"SVdup"     = "brown2"
)


n_genes = length(unique(gene_combined_events_EDIT$gene))
plt_title = paste0("min_freq = ", MIN_COHORT_FREQUENCY * 100, "%; n = ", n_genes)

# label known cancer genes
gene_combined_events_EDIT$is_cancer_gene = NA
gene_combined_events_EDIT[as.character(gene_combined_events_EDIT$gene) %in% cgc2024, "is_cancer_gene"] = "CGC"
gene_combined_events_EDIT[as.character(gene_combined_events_EDIT$gene) %in% prostate_cancer_genes, "is_cancer_gene"] = "PrCa"
gene_combined_events_EDIT[duplicated(as.character(gene_combined_events_EDIT$gene)), "is_cancer_gene"] = NA
cancer_gene_symbol = c("PrCa" = 1, "CGC" = 4)


plt_events = ggplot(gene_combined_events_EDIT, aes(y = gene, fill = ANNOT)) + 
		geom_bar(position = "fill") + 
		scale_fill_manual(values = mut_colors) +
		geom_point(data = gene_combined_events_EDIT, aes(shape = is_cancer_gene), x = -0.02) +
		scale_shape_manual("gene", values = cancer_gene_symbol) + 
		plot_theme() + 
		theme(legend.position = "right", 
				legend.key.size = unit(0.5, 'cm'), 
				legend.text = element_text(size = 6), 
				axis.text.y = element_text(size = 6)) + 
		ggtitle(NULL, plt_title)
		
# cohort freq
dfr_cohort_freq = signif(sort(table(gene_combined_events_EDIT$gene)) / length(all_patients), 2) * 100
dfr_cohort_freq = data.frame(gene = names(dfr_cohort_freq), freq = c(dfr_cohort_freq), stringsAsFactors = FALSE)
dfr_cohort_freq = dfr_cohort_freq[dfr_cohort_freq$gene %in% genes_to_select, ]

plt_cohort_freq = ggplot(gene_combined_events_EDIT, aes(y = gene)) + 
		geom_bar() + 
		plot_theme() +
		geom_text(data = dfr_cohort_freq, aes(x = 200, y = gene, label = freq)) + 
		scale_y_discrete(NULL, labels = NULL) 

save(gene_combined_events_EDIT, file = pff("gene_combined_events_EDIT.rsav"))


# collect additional clinical variables for each driver
all_drivers = as.character(unique(gene_combined_events_EDIT$gene))
all_drivers = setdiff(all_drivers, NA)
load(file = pff("patient_clinical_data.rsav"))


# this_driver = all_drivers[[i]]; clin_var = "GleasonGroup_4plus"
get_clinvars = function(this_driver, gene_combined_events_EDIT, clin_var, patient_clinical_data) {
	
	this_patients = gene_combined_events_EDIT[gene_combined_events_EDIT$gene == this_driver, "patient"]
	clin_vals = patient_clinical_data[patient_clinical_data$patient %in% this_patients, clin_var]
	dfr = data.frame(this_driver, clin_var, clin_vals, stringsAsFactors = FALSE)
	dfr = dfr[!is.na(dfr$clin_vals),]
	dfr
}

# for each driver get distribution of Gleason GG
drivers_gleasons = do.call(rbind, lapply(all_drivers, get_clinvars, gene_combined_events_EDIT, "Gleason_group", patient_clinical_data))
drivers_gleasons$this_driver = factor(drivers_gleasons$this_driver, levels = genes_by_freq)

plt_gleason_freq = ggplot(drivers_gleasons, aes(this_driver, fill = factor(clin_vals))) + 
		geom_bar(position = "fill") + 
		coord_flip() + 
		scale_fill_brewer("grade", palette = "OrRd", na.value = "grey") + 
		plot_theme() + 
		scale_x_discrete(NULL, labels = NULL, breaks = NULL) + 
		scale_y_continuous (NULL, lim = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%")) +
		theme(legend.position = "right", axis.text.y = element_text(size = 4), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size = 6),
				axis.text.x = element_text(angle = 0, hjust = 0)) + 
		ggtitle(NULL, "Grade")


# for each driver get distribution of T stage
drivers_stage = do.call(rbind, lapply(all_drivers, get_clinvars, gene_combined_events_EDIT, "T_stage", patient_clinical_data))
drivers_stage$this_driver = factor(drivers_stage$this_driver, levels = genes_by_freq)

plt_stage_freq = ggplot(drivers_stage, aes(this_driver, fill = factor(clin_vals))) + 
		geom_bar(position = "fill") + 
		coord_flip() + 
		scale_fill_brewer("stage", palette = "PuRd", na.value = "grey") + 
		plot_theme() + 
		scale_x_discrete(NULL, labels = NULL, breaks = NULL) + 
		scale_y_continuous (NULL, lim = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%")) +
		theme(legend.position = "right", axis.text.y = element_text(size = 4), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size = 6),
				axis.text.x = element_text(angle = 0, hjust = 0)) + 
		ggtitle(NULL, "Stage")


# for each driver get mets biology
drivers_metsbio = do.call(rbind, lapply(all_drivers, get_clinvars, gene_combined_events_EDIT, "new_metastatic_biology_indicator", patient_clinical_data))
drivers_metsbio$this_driver = factor(drivers_metsbio$this_driver, levels = genes_by_freq)
drivers_metsbio$clin_vals = factor(drivers_metsbio$clin_vals, levels = c("no_mets_biol", "mets_biol"))

plt_metsbio_freq = ggplot(drivers_metsbio, aes(this_driver, fill = factor(clin_vals))) + 
		geom_bar(position = "fill") + 
		coord_flip() + 
		scale_fill_manual(NULL, values = c("mets_biol" = "darkred", "no_mets_biol" = "lightgrey")) + 
		plot_theme() + 
		scale_x_discrete(NULL, labels = NULL, breaks = NULL) + 
		scale_y_continuous (NULL, lim = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%")) +
		theme(legend.position = "right", axis.text.y = element_text(size = 4), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size = 6), 
				axis.text.x = element_text(angle = 0, hjust = 0)) + 
		ggtitle(NULL, "MetsBio")


# for each driver get mets biology
drivers_EO = do.call(rbind, lapply(all_drivers, get_clinvars, gene_combined_events_EDIT, "early_onset", patient_clinical_data))
drivers_EO$this_driver = factor(drivers_EO$this_driver, levels = genes_by_freq)
drivers_EO$clin_vals = factor(drivers_EO$clin_vals, levels = c("FALSE", "TRUE"))


plt_EO_freq = ggplot(drivers_EO, aes(this_driver, fill = factor(clin_vals))) + 
		geom_bar(position = "fill") + 
		coord_flip() + 
		scale_fill_manual(NULL, values = c("TRUE" = "cornflowerblue", "FALSE" = "lightgrey")) + 
		plot_theme() + 
		scale_x_discrete(NULL, labels = NULL, breaks = NULL) + 
		scale_y_continuous (NULL, lim = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%")) +
		theme(legend.position = "right", axis.text.y = element_text(size = 4), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size = 6),
				axis.text.x = element_text(angle = 0, hjust = 0)) + 
		ggtitle(NULL, "EO")
		
		
						
combo_plot = (plt_events | plt_cohort_freq | plt_gleason_freq | plt_stage_freq | plt_metsbio_freq | plt_EO_freq) +
		plot_layout(widths = unit(c(0.6, 0.1, 0.1, 0.1, 0.1, 0.1), 'null')) +
		plot_annotation(title = plt_title)

fname = pff("selected_driver_mutations_barplot.pdf")
ggsave(combo_plot, file = fname, width = 13, height = 7)
file_open_call2(fname)
