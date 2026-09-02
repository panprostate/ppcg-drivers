source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(igraph)


FDR_CUTOFF_FOR_CNA_ENRICH = 0.05
FDR_CUTOFF_FOR_CNA_COOCCUR = 0.05
N_MCAPPLY_CORES = 8

load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes = unique(c(cgc2024, prostate_cancer_genes))


# testing co-occurrence of same CNA segments affecting a pair of genes
test_coCNA = function(i, gene_pairs, gene_focal_events, all_CNA_ids) {
	
	if (i %% 1000 == 1) cat("O")

	gene1 =  gene_pairs[i, 1]
	gene2 =  gene_pairs[i, 2]

	gene1_cnas = gene_focal_events[[gene1]]
	gene2_cnas = gene_focal_events[[gene2]]
	
	pval = fisher.test(all_CNA_ids %in% gene1_cnas, all_CNA_ids %in% gene2_cnas, alt = "greater")$p.value
	
	fraction_2in1 = sum(gene2_cnas %in% gene1_cnas) / length(gene2_cnas)
	fraction_1in2 = sum(gene1_cnas %in% gene2_cnas) / length(gene1_cnas)
	
	n_gene1 = length(gene1_cnas)
	n_gene2 = length(gene2_cnas)
	n_gene_both = length(intersect(gene1_cnas, gene2_cnas))

	data.frame(gene1, gene2, pval, fraction_2in1, fraction_1in2, n_gene1, n_gene2, n_gene_both, stringsAsFactors = FALSE)
}

# formats the co-altered gene modules into a common table 
collapse_module_focalCNAs = function(module, focalAMP_enriched_signf, annot, cgc2024, prostate_cancer_genes) {
	
	focalCNAs_here = focalAMP_enriched_signf[focalAMP_enriched_signf$gene %in% module,]
	focalCNAs_here = focalCNAs_here[order(focalCNAs_here$pval),]

	genes_here = focalCNAs_here$gene
	cancer_genes_here = intersect(genes_here, cgc2024)
	prostate_cancer_genes_here = intersect(genes_here, prostate_cancer_genes)

	concat_genes_here = paste(genes_here, collapse = ",")
	cgc_genes = paste(cancer_genes_here, collapse = ",")
	prca_genes = paste(prostate_cancer_genes_here, collapse = ",")
	TGT_GENE = prca_genes 
	
	# representative is either the first (topP) cancer gene, or if not, just the first gene
	repping_gene = genes_here[1]
	gene_as_id = concat_genes_here
	if (length(cancer_genes_here) > 0) {
		repping_gene = cancer_genes_here[1]
		gene_as_id = cgc_genes
	}
	# take pvals, FC etc for that repping gene
	repping_record = focalCNAs_here[focalCNAs_here$gene == repping_gene,][1,, drop = FALSE]
	pval = repping_record$pval[1]
	fdr = repping_record$fdr[1]
	fc = repping_record$fc[1]
	frac_CNA = repping_record$frac_CNA[1]
	mean_rel_cn = repping_record$mean_rel_cn[1]
	patients = unique(unlist(strsplit(repping_record$patient, split = ",")))
	n_patients = length(patients)
	patients = paste(patients, collapse = ",")

	data.frame(gene = gene_as_id, pval, frac_CNA, fc, n_patients, mean_rel_cn, 
			patients, annot, fdr, all_genes = concat_genes_here,
			cgc_genes, prca_genes, TGT_GENE,
			stringsAsFactors = FALSE)
}

get_coCNA_modules = function(annot, cna_annots_to_select, gene_CNA_enriched, gene_cna_annots, gr_genes) {

	gene_enriched_signf = gene_CNA_enriched[gene_CNA_enriched$fdr < FDR_CUTOFF_FOR_CNA_ENRICH & gene_CNA_enriched$annot == annot,]
	gene_cna_annots_here = gene_cna_annots[gene_cna_annots$annot %in% cna_annots_to_select,]
	
	all_genes = unique(gene_enriched_signf$gene)
	gene_pairs = t(combn(all_genes, 2))
	cat(annot, " ", nrow(gene_pairs), " before\n")
	
	# gene pairs only between the same chr
	gene_to_chr = structure(as.character(seqnames(gr_genes)), names = gr_genes$symbol)
	gene_chrs = cbind(gene_to_chr[gene_pairs[,1]], gene_to_chr[gene_pairs[,2]])
	gene_pairs = gene_pairs[gene_chrs[,1] == gene_chrs[,2],]
	
	all_CNA_ids = unique(gene_cna_annots_here$CNA_id)
	
	gene_focal_events = split(gene_cna_annots_here$CNA_id, gene_cna_annots_here$gene)
	gene_focal_events = sapply(gene_focal_events, unique)
	
	cat(annot, " ", nrow(gene_pairs), " after:\n")
	all_gene_coCNAs = do.call(rbind, mclapply(1:nrow(gene_pairs), test_coCNA, gene_pairs, gene_focal_events, all_CNA_ids, 
			mc.cores = N_MCAPPLY_CORES, mc.preschedule = TRUE))
			
	all_gene_coCNAs$fdr = p.adjust(all_gene_coCNAs$pval, method = "fdr")
	all_gene_coCNAs_signf = all_gene_coCNAs[all_gene_coCNAs$fdr < FDR_CUTOFF_FOR_CNA_COOCCUR,]
	
	# modules are connected components in graphs
	# make sure all genes are included in the network - add self-loops to every gene
	all_coCNA_pairs = as.matrix(unique(all_gene_coCNAs_signf[, c("gene1", "gene2")]))
	genes_with_coCNA = unique(c(all_coCNA_pairs))
	
	all_singletons = setdiff(all_genes, genes_with_coCNA)
	all_singleton_pairs = as.matrix(unique(cbind(all_singletons, all_singletons)))
	all_gene_pairs = rbind(all_coCNA_pairs, all_singleton_pairs)
	colnames(all_gene_pairs) = c("to", "from")
	cat("\n", annot, "n_pairs:", nrow(all_gene_pairs), "\n")
	
	coAMP_graph = igraph::graph_from_data_frame(d = all_gene_pairs,  directed = FALSE) 
	coAMP_module_object = components(coAMP_graph)
	modules_from_merged_lists = data.frame(gene = names(components(coAMP_graph)$membership), module_no = components(coAMP_graph)$membership, 
			row.names = NULL, 
			stringsAsFactors = FALSE)
	modules_from_merged_lists = split(modules_from_merged_lists$gene, modules_from_merged_lists$module_no)
	
	# format modules into a table
	cat("\n", annot, "n_modules:", length(modules_from_merged_lists), "\n")
	CNA_modules = do.call(rbind, lapply(modules_from_merged_lists, 
			collapse_module_focalCNAs, gene_enriched_signf, annot, cancer_genes, prostate_cancer_genes))
	CNA_modules = CNA_modules[order(CNA_modules$pval),]
	CNA_modules
}


load(pff("gene_CNA_enriched.rsav"))
load(pff("gene_cna_annots.rsav"))

load(file = pff("gr_prepared_elements.rsav"))
gr_genes = gr_prepared_elements[grep("^CDSgene", gr_prepared_elements$id)]


modules_high_gain = get_coCNA_modules("high_gain", "high_gain", gene_CNA_enriched, gene_cna_annots, gr_genes)
modules_full_loss = get_coCNA_modules("full_loss", "full_loss", gene_CNA_enriched, gene_cna_annots, gr_genes)
modules_all_gain = get_coCNA_modules("all_gain", c("gain", "high_gain"), gene_CNA_enriched, gene_cna_annots, gr_genes)
modules_all_loss = get_coCNA_modules("all_loss", c("loss", "full_loss"), gene_CNA_enriched, gene_cna_annots, gr_genes)

CNA_modules = rbind(modules_high_gain, modules_full_loss, modules_all_gain, modules_all_loss)
save(CNA_modules, file = pff("CNA_modules.rsav"))

