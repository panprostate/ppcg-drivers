source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(qs)
library(SummarizedExperiment)

load(pff("all_patients.rsav"))
load(pff("tracking_sheet.rsav"))
WGS_preferred_samples = tracking_sheet$WGS_AssayID[tracking_sheet$selected_one_sample_per_donor]

mrna_fname = paste0("DATA_USED__", this_timestamp, "/RNA_2023-04-06/Counts_Release_1.2/PPCG_RNAseq_Release1.1.qs")
qdat = qread(mrna_fname)
expr_mat = assays(qdat)[['RUV_III_PRPS']]
sample_meta = colData(qdat)

# select preferred samples from WGS and those with RNAseq
sample_meta = colData(qdat)
which_RNAseq = sample_meta$Platform.RNAseq..Affymetrix..Beadchip...... %in% c("RNA_seq", "RNAseq")
which_preferred_PPCG_sample = sample_meta$Matching_WGS_Sample.ID.ADD.NA. %in% WGS_preferred_samples
which_select = which_RNAseq & which_preferred_PPCG_sample

sample_meta = sample_meta[which_select,]
expr_mat = expr_mat[,which_select]

# remove PPCG donor duplicates: keep the first random sample
which_not_duplic = !duplicated(sample_meta$PPCG_Donor_ID)

sample_meta = sample_meta[which_not_duplic,]
expr_mat = expr_mat[,which_not_duplic]

# donor ID assigned to columns
colnames(expr_mat) = sample_meta$PPCG_Donor_ID
gene_ids = strsplit(rownames(qdat), s = "\\|\\|")
which_has_symbol = sapply(gene_ids, length) == 2
expr_mat = expr_mat[which_has_symbol,]

rownames(expr_mat) = gsub("(.+)\\|\\|(.+)", "\\2", rownames(expr_mat))
expr_mat = expr_mat[!duplicated(rownames(expr_mat)),]
expr_sample_meta = sample_meta
save(expr_mat, file = pff("expr_mat.rsav"))
save(expr_sample_meta, file = pff("expr_sample_meta.rsav"))


##
# get lncRNA matrix separately
##
mrna_fname = paste0("DATA_USED__", this_timestamp, "/PPCG_RNAseq_rel1.2_tier2_2025-01-16/PPCG_RNAseqData_Tier2.qs")
qdat = qread(mrna_fname)
expr_mat = assays(qdat)[['normalized']]
sample_meta = colData(qdat)

# select preferred samples from WGS and those with RNAseq
which_RNAseq = sample_meta$Platform.RNAseq..Affymetrix..Beadchip...... %in% c("RNA_seq", "RNAseq")
which_preferred_PPCG_sample = sample_meta$Matching_WGS_Sample.ID.ADD.NA. %in% WGS_preferred_samples
which_select = which_RNAseq & which_preferred_PPCG_sample
sample_meta = sample_meta[which_select,]
expr_mat = expr_mat[,which_select]
# remove PPCG donor duplicates: keep the first random sample
which_not_duplic = !duplicated(sample_meta$PPCG_Donor_ID)
sample_meta = sample_meta[which_not_duplic,]
expr_mat = expr_mat[,which_not_duplic]

# donor ID assigned to columns
colnames(expr_mat) = sample_meta$PPCG_Donor_ID
gene_ids = strsplit(rownames(qdat), s = "\\|\\|")
which_has_symbol = sapply(gene_ids, length) == 2
expr_mat = expr_mat[which_has_symbol,]

rownames(expr_mat) = gsub("(.+)\\|\\|(.+)", "\\2", rownames(expr_mat))
expr_mat = expr_mat[!duplicated(rownames(expr_mat)),]
save(expr_mat, file = pff("expr_mat_with_LNCs.rsav"))
