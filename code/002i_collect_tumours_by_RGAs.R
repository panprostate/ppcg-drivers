source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")

load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.005
MIN_N_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)
PVALUE_CUTOFF = 0.05

get_gene_pval_samples = function (gene, mut_types, results_signf_merged_annot) {
	
	mut_types_collapse = paste(mut_types, collapse = ",")	
	res_this = results_signf_merged_annot[
				results_signf_merged_annot$annots_MAIN == gene & 
				results_signf_merged_annot$mut_type %in% mut_types,]
	if (nrow(res_this) == 0) {
		dfr = data.frame(gene, mut_type = mut_types_collapse, fdr = 1, n_patients = 0, concat_patients = "", stringsAsFactors = FALSE)
		return(dfr)
	}
				
	fdr = min(res_this$fdr)
	mut_patients =  unique(unlist(strsplit(res_this[,"patient_ids"], s=",")))
	n_patients = length(mut_patients)
	concat_patients = paste(mut_patients, collapse = ",")
	data.frame(gene, mut_type = mut_types_collapse, fdr, n_patients, concat_patients, stringsAsFactors = FALSE)
}



# samples by mutation type for each gene
strplit_samples = function(this_res, mut_type) {
	samples = strsplit(this_res[this_res$mut_type == mut_type, "concat_patients"], split = ",")[[1]]
	if (length(samples) == 0) {
		return(NULL)	
	}
	cbind(samples, mut_type)
}

get_sample_mut_types = function(gene, merged_pvals_counts) {

	this_res = merged_pvals_counts[merged_pvals_counts$gene == gene,]
		
	samples_CDS = strplit_samples(this_res, "SNV_CDS")
	samples_NC = strplit_samples(this_res, "SNV_NC")
	samples_SV = strplit_samples(this_res, "SV")
	
	samples_combined = rbind(samples_CDS, samples_NC, samples_SV)
	samples_combined = data.frame(gene, samples_combined, stringsAsFactors = FALSE)
	samples_combined$mut_type2 = gsub("SNV_CDS|SNV_NC", "SNV", samples_combined$mut_type)
	samples_combined
}



load(pff("results_signf_merged_annot.rsav"))
all_elements = unique(results_signf_merged_annot$annots_MAIN)

results_SV = do.call(rbind, 
		lapply(all_elements, get_gene_pval_samples, "SV", results_signf_merged_annot))
results_NC = do.call(rbind, 
		lapply(all_elements, get_gene_pval_samples, "SNV_NC", results_signf_merged_annot))
results_CDS = do.call(rbind, 
		lapply(all_elements, get_gene_pval_samples, "SNV_CDS", results_signf_merged_annot))

merged_pvals_counts = rbind(results_SV, results_NC, results_CDS)
merged_pvals_counts = merged_pvals_counts[order(merged_pvals_counts$fdr),]
save(merged_pvals_counts, file = pff("merged_pvals_counts.rsav"))


# SNV, SV lists
merged_sample_lists = do.call(rbind, lapply(all_elements, get_sample_mut_types, merged_pvals_counts))
driver_patient_sets = split(merged_sample_lists$samples, paste(merged_sample_lists$gene, merged_sample_lists$mut_type, sep = "__"))


# for CNA events, take the annotated genes from modules/peaks and take high-confidence segments and tumors
FDR_CUTOFF = 0.05

get_cna_patients = function(annot, CNA_modules_annotated, gene_CNA_enriched, annot_tag, fdr_cutoff) {
	genes_index_modules = CNA_modules_annotated[CNA_modules_annotated$annot == annot, "TGT_GENE"]
	# sometimes two genes per module, split them
	genes_index_modules = unlist(strsplit(genes_index_modules, split = ","))
	gene_CNA_enriched_here = gene_CNA_enriched[
			gene_CNA_enriched$annot == annot & 
			gene_CNA_enriched$gene %in% genes_index_modules & 
			gene_CNA_enriched$fdr < fdr_cutoff,]
	
	CNA_patients_sets = strsplit(gene_CNA_enriched_here$patients, s = ",")
	names(CNA_patients_sets) = paste0(gene_CNA_enriched_here$gene, annot_tag)
	list(CNA_patients_sets, gene_CNA_enriched_here)
}


load(pff("gene_CNA_enriched.rsav"))
load(pff("CNA_modules_w_tgt_genes.rsav"))
load(pff("GISTIC_peaks_w_tgt_genes.rsav"))

CNA_modules_w_tgt_genes$annot = gsub("(.+)__(.+)__(.+)", "\\3", CNA_modules_w_tgt_genes$id)
GISTIC_peaks_w_tgt_genes$TGT_GENE = GISTIC_peaks_w_tgt_genes$TGT_gene
GISTIC_modules_combined = rbind(
		unique(GISTIC_peaks_w_tgt_genes[, c("TGT_GENE", "annot")]), 
		unique(CNA_modules_w_tgt_genes[, c("TGT_GENE", "annot")]))
GISTIC_modules_combined = GISTIC_modules_combined[GISTIC_modules_combined$TGT_GENE != "",]

focalAMP_res = get_cna_patients("high_gain", GISTIC_modules_combined, gene_CNA_enriched, "__hAMP", FDR_CUTOFF)
focalAMP_patients_sets = focalAMP_res[[1]]
focalAMP_gene_dfr = focalAMP_res[[2]]

loss_res = get_cna_patients("all_loss", GISTIC_modules_combined, gene_CNA_enriched, "__loss", FDR_CUTOFF)
loss_patients_sets = loss_res[[1]]
loss_gene_dfr = loss_res[[2]]

gain_res = get_cna_patients("all_gain", GISTIC_modules_combined, gene_CNA_enriched, "__gain", FDR_CUTOFF)
gain_patients_sets = gain_res[[1]]
gain_gene_dfr = gain_res[[2]]

gene_CNA_module_dfr = rbind(focalAMP_gene_dfr, loss_gene_dfr, gain_gene_dfr)
save(gene_CNA_module_dfr, file = pff("gene_CNA_module_dfr.rsav"))


# biallelic lists - keep only driver genes here
load(file = pff("BI_results.rsav"))
BI_results = BI_results[BI_results$is_cancer_gene,]
biallelic_patients_sets = strsplit(BI_results$patients, split = ",")
names(biallelic_patients_sets) = paste0(BI_results$gene, "__BI")

patient_sets_for_drivers = c(driver_patient_sets, focalAMP_patients_sets, biallelic_patients_sets, loss_patients_sets, gain_patients_sets)
patient_sets_for_drivers = patient_sets_for_drivers[sapply(patient_sets_for_drivers, length) >= MIN_N_PATIENTS]
save(patient_sets_for_drivers, file = pff("patient_sets_for_drivers.rsav"))


get_secondary_targets = function(i, results_signf_merged_annot, genes_from_column) {
	cat(i, " ")

	tgt_genes = strsplit(results_signf_merged_annot[i, genes_from_column], split = ",")[[1]]
	patients_involved = strsplit(results_signf_merged_annot[i, "patient_ids"], split = ",")[[1]]
	annot_gene_main = results_signf_merged_annot[i, "annots_MAIN"]
	mut_type = results_signf_merged_annot[i, "mut_type"]
	tgt_genes = setdiff(tgt_genes, c("", NA, annot_gene_main))
	if (length(tgt_genes) == 0){
		return(NULL)
	}
	
	res = lapply(tgt_genes, function(g) patients_involved)
	names(res) = paste(tgt_genes, mut_type, annot_gene_main, sep = "__")
	res = data.frame(as.matrix(stack(res)), stringsAsFactors = FALSE)
	colnames(res) = c("patient", "id")
	res
}


load(file = pff("results_signf_merged_annot.rsav"))
results_signf_merged_annot = results_signf_merged_annot[results_signf_merged_annot$annots_loops != "", ]
results_signf_merged_annot = results_signf_merged_annot[results_signf_merged_annot$mut_type != "SNV_CDS", ]
patient_sets_for_loops = do.call(rbind, lapply(1:nrow(results_signf_merged_annot), get_secondary_targets, results_signf_merged_annot, "annots_loops"))
patient_sets_for_loops = split(patient_sets_for_loops$patient, patient_sets_for_loops$id)
save(patient_sets_for_loops, file = pff("patient_sets_for_loops.rsav"))


load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes_combined = unique(c(cgc2024, prostate_cancer_genes))


load(pff("all_patients.rsav"))
load(file = pff("patient_sets_for_drivers.rsav"))
patient_sets_for_drivers = patient_sets_for_drivers[rev(order(sapply(patient_sets_for_drivers, length)))]

n_spls =  sapply(patient_sets_for_drivers, length)
pct_spls = 100*signif(n_spls / length(all_patients), 2)
gene = gsub("(.+)__(.+)", "\\1", names(n_spls))
mut_type = gsub("(.+)__(.+)", "\\2", names(n_spls))

dfr = data.frame(gene, mut_type, event = names(n_spls), pct_spls, n_spls, stringsAsFactors = FALSE)
rownames(dfr) = NULL

fname = pff("figures/patient_sets_for_drivers__summary_table.csv")
write.csv(dfr, fname)
file_open_call2(fname)
