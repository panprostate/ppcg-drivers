source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")


##
# get BI lists from mets for the P M comparison
##
load(file = pff("allmets_tumor_ids.rsav"))

# select the genes we already select for the BI list
load(file = pff("BI_results.rsav"))
genes_with_BI = BI_results[BI_results$is_cancer_gene, "gene"]
BI_patient_sets = read.delim(
		paste0("DATA_USED__", this_timestamp, "/", "biallelic/Daria_2026-04-30/Bialellic_all_samples_Aprl2026.tsv"),
		stringsAsFactors = FALSE)

BI_patient_sets_mets = BI_patient_sets[BI_patient_sets$gene %in% genes_with_BI & BI_patient_sets$sample_id %in% allmets_tumor_ids, ]
BI_patient_sets_mets = split(BI_patient_sets_mets$sample_id, BI_patient_sets_mets$gene)
BI_patient_sets_mets = BI_patient_sets_mets[genes_with_BI]
names(BI_patient_sets_mets) = genes_with_BI
save(BI_patient_sets_mets, file = pff("BI_patient_sets_mets.rsav"))
