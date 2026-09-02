source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(ggrepel)
library(patchwork)

FDR_CUTOFF = 0.05

load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes = unique(c(cgc2024, prostate_cancer_genes))


test_binom_enrichment = function(gene, gene_annots, n_patients, p_expected) {
	
	if (runif(1) < 0.001) cat("o")
	patients = unique(gene_annots[gene_annots$gene == gene, "patient"])
	n_patients_here = length(patients)
	
	mean_rel_cn = mean(gene_annots[gene_annots$gene == gene, "rel_cn"])
	
	pval = pbinom(n_patients_here - 1, n_patients, p_expected, lower.tail = FALSE)
	frac_CNA = n_patients_here / n_patients
	fc = frac_CNA / p_expected
	
	data.frame(gene, pval, frac_CNA, fc, n_patients = n_patients_here, mean_rel_cn, patients = paste(patients, collapse = ","), 
			stringsAsFactors = FALSE)
}


compute_gene_CNA_enrichment = function(CNA_tag, genes_with_full_CNAs, n_genes, n_patients, cancer_genes) {
	
	genes_with_CNAs_here = genes_with_full_CNAs[genes_with_full_CNAs$tag == CNA_tag,]
	n_events = nrow(unique(genes_with_CNAs_here[, c("gene", "patient")]))
	p_CNAs_expected = n_events / (n_genes * n_patients)
		
	genes_with_CNAs_enriched = do.call(rbind, mclapply(gr_genes$symbol, test_binom_enrichment, genes_with_CNAs_here, n_patients, p_CNAs_expected, 
			mc.cores = 8))
	genes_with_CNAs_enriched$annot = CNA_tag
	genes_with_CNAs_enriched$fdr = p.adjust(genes_with_CNAs_enriched$pval, method = "fdr")
	genes_with_CNAs_enriched = genes_with_CNAs_enriched[order(genes_with_CNAs_enriched$fdr),]
	
	genes_with_CNAs_enriched$is_cancer_gene = genes_with_CNAs_enriched$gene %in% cancer_genes
	genes_with_CNAs_enriched$p_event = p_CNAs_expected
	genes_with_CNAs_enriched
}




load(file = pff("genes_with_full_CNAs.rsav"))

# background gene list
load(file = pff("gr_prepared_elements.rsav"))
gr_genes = gr_prepared_elements[grep("^CDSgene", gr_prepared_elements$id)]
n_genes = length(gr_genes)

# background patient list
load(file = pff("all_patients.rsav"))
n_patients = length(unique(all_patients))

genes_enriched__high_gain = compute_gene_CNA_enrichment("high_gain", genes_with_full_CNAs, n_genes, n_patients, cancer_genes)
genes_enriched__full_loss = compute_gene_CNA_enrichment("full_loss", genes_with_full_CNAs, n_genes, n_patients, cancer_genes)
genes_enriched__all_gain = compute_gene_CNA_enrichment("all_gain", genes_with_full_CNAs, n_genes, n_patients, cancer_genes)
genes_enriched__all_loss = compute_gene_CNA_enrichment("all_loss", genes_with_full_CNAs, n_genes, n_patients, cancer_genes)

gene_CNA_enriched = rbind(genes_enriched__high_gain, genes_enriched__full_loss, genes_enriched__all_gain, genes_enriched__all_loss)

# select FDR-significant genes only
gene_CNA_enriched = gene_CNA_enriched[gene_CNA_enriched$fdr < FDR_CUTOFF,]
save(gene_CNA_enriched, file = pff("gene_CNA_enriched.rsav"))

fname = pff("figures/gene_CNA_enriched.csv")
write.csv(gene_CNA_enriched, file = fname)
file_open_call2(fname)

