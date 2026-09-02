source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library("dndscv")
library("gtools")

FDR_SELECT = 0.05

load(file = pff("prepared_variants__SNV_indel.rsav"))

mutations = prepared_variants
mutations = mutations[, c("patient", "chr", "pos1", "ref", "alt")]
colnames(mutations) = c("sampleID", "chr", "pos", "ref", "alt")
mutations$pos = as.numeric(mutations$pos)
mutations$chr = gsub("chr", "", mutations$chr)

dnds_results_all = dndscv(mutations)
save(dnds_results_all, file = pff("dnds_results_all.rsav"))

dnds_results = dnds_results_all$sel_cv[dnds_results_all$sel_cv$qglobal_cv < FDR_SELECT, ]
save(dnds_results, file = pff("dnds_results.rsav"))

# impact of mutations - selected genes
# label variant annotations, remove silent SNVs
dnds_results_annotations = dnds_results_all$annotmuts[dnds_results_all$annotmuts$gene %in% dnds_results$gene,]
dnds_results_annotations = dnds_variant_filter_label(dnds_results_annotations)
save(dnds_results_annotations, file = pff("dnds_results_annotations.rsav"))

# impact of mutations - all genes
# label variant annotations, remove silent SNVs
dnds_results_annotations_all = dnds_results_all$annotmuts
dnds_results_annotations_all = dnds_variant_filter_label(dnds_results_annotations_all)
save(dnds_results_annotations_all, file = pff("dnds_results_annotations_all.rsav"))

#
# run dn/ds on variants from mets to make annotation consistent with primary tumors
#
load(file = pff("variants_mets__SNV_indel.rsav"))
mutations = variants_mets__SNV_indel[, c("tumor_id", "chr", "pos1", "ref", "alt")]
colnames(mutations) = c("sampleID", "chr", "pos", "ref", "alt")
mutations$pos = as.numeric(mutations$pos)
mutations$chr = gsub("chr", "", mutations$chr)

dnds_results_mets_all = dndscv(mutations)
save(dnds_results_mets_all, file = pff("dnds_results_mets_all.rsav"))

dnds_results_mets_annotations_all = dnds_results_mets_all$annotmuts
dnds_results_mets_annotations_all = dnds_variant_filter_label(dnds_results_mets_annotations_all)
save(dnds_results_mets_annotations_all, file = pff("dnds_results_mets_annotations_all.rsav"))

#
# plot mutations for primary samples
#
load(file = pff("dnds_results_annotations.rsav")) 
load(file = pff("dnds_results.rsav"))

dnds_results_annotations$patient_gene = paste(dnds_results_annotations$sampleID, dnds_results_annotations$gene, sep = "__")

# collapse into "multiple" if several types of mutations found per tumor
muts_patient_gene = sapply(by(dnds_results_annotations$impact, dnds_results_annotations$patient_gene, unique), paste, collapse = "__")
muts_patient_gene = data.frame(patient_gene = names(muts_patient_gene), annot = muts_patient_gene, stringsAsFactors = FALSE)
muts_patient_gene$patient = gsub("(.+)__(.+)", "\\1", muts_patient_gene$patient_gene)
muts_patient_gene$gene = gsub("(.+)__(.+)", "\\2", muts_patient_gene$patient_gene)
muts_patient_gene$annot [grep("__", muts_patient_gene$annot)] = "multiple"

gene_fdrs = structure(names = dnds_results$gene, dnds_results$qglobal_cv)

# order genes and annotations by frequency
gene_ordering = names(sort(table(muts_patient_gene$gene)))
muts_patient_gene$gene = factor(muts_patient_gene$gene, levels = gene_ordering)
annot_ordering = rev(names(sort(table(muts_patient_gene$annot))))
muts_patient_gene$annot = factor(muts_patient_gene$annot, levels = annot_ordering)

muts_patient_gene$fdr = gene_fdrs[as.character(muts_patient_gene$gene)]
muts_patient_gene$fdr_stars = stars.pval(muts_patient_gene$fdr)

plt = ggplot(muts_patient_gene, aes(gene, fill = annot, label = fdr_stars)) +
		geom_bar() +
		plot_theme() + 
		geom_text(y = -1) +
		scale_fill_brewer(palette = "Set1") + 
		coord_flip() + 
		ggtitle(paste0("Protein-coding drivers, DNDS_CV"), paste0("FDR < ", FDR_SELECT, "; n = ", length(unique(muts_patient_gene$gene))))
		
fname = pff("figures/DnDs_proteincoding_driver_genes.pdf")
ggsave(plt, file = fname, width = 6.3)
file_open_call2(fname)



