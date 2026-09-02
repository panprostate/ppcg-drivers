source("/PATH_TO_WORKING_DIR/data")
date_tag = "TODAY_DATE"

suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))


# use the same covariates throughout the study
MFS_COVARIATES = c("early_onset", "T_stage", "GleasonGroup_3plus")
MFS_ANALYSIS_TAG = "PPCG__Age_Stage_Grade"


# use the same cluster number and distance/method throughout
SELECT_CLUSTER_METHOD = "wardD2__euclidean.dist"
SELECT_N_CLUSTERS = 4

# assign cluster colors to be the same for ggplot, R, APW
CLUSTER_COLORS = c("forestgreen", "#FFDB00", "#4ed119", "#c20a0a", "#33CCFF", "lightpink1", "red4", "yellow")



pff = function(x = "") {
	paste0("./", date_tag, "/", paste0(x, collapse=""), sep="")
}

opff = function(x) {
	system(paste("open", pff(x)))
}

this_timestamp = gsub("/|\\./", "", pff())



file_open_call = function(fname) {
	fname1 = tail(strsplit(fname, '/')[[1]], 1)
	cat(paste0("\n\nscp -C $hn:", getwd(), "/", fname, " ./ && open ", fname1, "\n\n"))
}

file_open_call2 = function(fname) {
	system(paste("cp ", fname, "~/"))
	fname1 = tail(strsplit(fname, '/')[[1]], 1)
	cat(paste0("\n\n scp -C $cw\\:", fname1, " ./ && open ", fname1, "\n\n"))
}



plot_theme = function(...) {
	theme_bw() + 
	theme(	
		plot.title = element_text(size = 22),
		plot.caption = element_text(size = 12),
		plot.subtitle = element_text(size = 16),
		axis.title = element_text(size = 18),
		axis.text.x = element_text(size = 12, 
				angle = 90, hjust = 1, vjust=0.5, color = "black"),
		axis.text.y = element_text(size = 12, color = "black"),
		legend.title = element_text(size = 16),
		legend.text = element_text(size = 14),
		...
	)
}




# https://gist.github.com/armish/564a65ab874a770e2c26
memoSort <- function(M) {
	geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
	scoreCol <- function(x) {
		score <- 0;
		for(i in 1:length(x)) {
			if(x[i]) {
				score <- score + 2^(length(x)-i);
			}
		}
		return(score);
	}
	scores <- apply(M[geneOrder, ], 2, scoreCol);
	sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
	return(M[geneOrder, sampleOrder]);
}








# gr_canvas; gr_genes = gr_this_genes; gr_CDS = gr_this_CDS

generate_gene_plot = function(gr_canvas, gr_genes, gr_CDS) {
	
	ov = findOverlaps(gr_genes, gr_canvas)	
	gr_genes_here = gr_genes[queryHits(ov)]
	gr_CDS_here = gr_CDS[gr_CDS$id %in% gr_genes_here$id]
	
	# assign a number to each gene and its CDS exons 
	gr_genes_here$rownum = as.numeric(factor(gr_genes_here$id))
	gene2rownum = structure(gr_genes_here$rownum, names = gr_genes_here$id)
	gr_CDS_here$rownum = gene2rownum[gr_CDS_here$id]
	
	# cut CDS exons and gene coords that do not fit the canvas
	end(gr_genes_here)[end(gr_genes_here) > end(gr_canvas)] = end(gr_canvas)
	start(gr_genes_here)[start(gr_genes_here) < start(gr_canvas)] = start(gr_canvas)
	
	ov = findOverlaps(gr_CDS_here, gr_canvas)
	gr_CDS_here = gr_CDS_here[queryHits(ov)]
	end(gr_CDS_here)[end(gr_CDS_here) > end(gr_canvas)] = end(gr_canvas)
	start(gr_CDS_here)[start(gr_CDS_here) < start(gr_canvas)] = start(gr_canvas)

	dfr_genes_here = as.data.frame(gr_genes_here, stringsAsFactors = FALSE)
	dfr_CDS_here = as.data.frame(gr_CDS_here, stringsAsFactors = FALSE)
	
	label_pos = (start(gr_canvas) + end(gr_canvas)) / 2
	
	plt = ggplot(dfr_genes_here, aes (xmin = start, xmax = end, ymin = rownum + 0.1, ymax = rownum + 0.9)) +
		geom_rect(fill = "lightgrey", color = "darkgrey") + 
		geom_rect(data = dfr_CDS_here, aes (xmin = start, xmax = end, ymin = rownum + 0.0, ymax = rownum + 1), color = "black") + 
		geom_text(data = dfr_genes_here, 
				aes(x = (start+end)/2, y = rownum + 0.5, label = paste(symbol, gene_strand)), 
				color = "purple", size = 2, vjust = 0.5) +
		coord_cartesian(xlim = c(start(gr_canvas), end(gr_canvas))) + 
		scale_y_continuous(NULL, breaks = NULL, labels = NULL) + 
		plot_theme()
	plt
}




# fname = pff("gobp_reac.gmt")
get_gmt = function(fname) {

	gmt_content = readLines(fname)
	gmt_content = strsplit(gmt_content, split = "\t")
	names(gmt_content) = sapply(gmt_content, function(x) paste(x[1:2], collapse = "__"))
	gmt_content = sapply(gmt_content, '[', -1:-2)
	
}





# gene = "MYC"
# cluster_id = "cluster_6"; gene = "ERG"
test_difexp_transcriptome = function(gene, expr_mat, samples_here, cluster_id, test_alternative ) {
	
	if (runif(1)<0.001) cat("o")
	
	vals = expr_mat[gene,]
	cluster_vals =  vals[samples_here]
	out_vals = vals[setdiff(names(vals), samples_here)]	
	
	tt = t.test(cluster_vals, out_vals, alternative = test_alternative)
	pval = tt$p.value
	stat = tt$statistic

	mean_clust = mean(cluster_vals)
	mean_out = mean(out_vals)
	fc = mean_clust / mean_out

	data.frame(gene, cluster_id, pval, fc, stat, mean_clust, mean_out)
}


difexp_per_cluster = function(cluster_id, patient_cluster_map, expr_mat, test_alternative) {
	cat("\n", cluster_id, "\n")
	samples_here = names(which(patient_cluster_map == cluster_id))
	all_genes = rownames(expr_mat)
	expr_res = do.call(rbind, mclapply(all_genes, test_difexp_transcriptome, expr_mat, samples_here, cluster_id, test_alternative, mc.cores = 8))
	expr_res = expr_res[order(expr_res$pval),]
	expr_res
}




get_gene_cytobands = function(gr_genes, cytobands) {
	
	gr_cytobands = GRanges(cytobands$V1, IRanges(cytobands$V2, cytobands$V3), id = cytobands$id)

	ov = findOverlaps(gr_genes, gr_cytobands)
	gene2cytoband = data.frame(
			gene = gr_genes[queryHits(ov)]$symbol, 
			cytoband = gr_cytobands[subjectHits(ov)]$id, 
			stringsAsFactors = FALSE)
	gene2cytoband = sapply(split(gene2cytoband$cytoband, gene2cytoband$gene), function(x) sort(x)[[1]])
	gene2cytoband = data.frame(gene = names(gene2cytoband), cytoband = gene2cytoband, stringsAsFactors = FALSE)

	gene2cytoband$cytoband = gsub("chr", "", gene2cytoband$cytoband)
	gene2cytoband$cytoband = gsub("(.+)\\.(.+)", "\\1", gene2cytoband$cytoband)
	
	gene2cytoband
}



## function to label variant annotations from dndscv

# re-label indels
# annot_table = dnds_results_annotations
dnds_variant_filter_label = function(annot_table) {
	
	# remove indels with missing nt changes
	annot_table = annot_table[!is.na(annot_table$ntchange),]

	which_to_relabel = which(annot_table$impact == "no-SNV")
	new_labels = gsub("(.+)-(.+)-(.+)", "\\3", annot_table[which_to_relabel, "ntchange"])
	new_labels = paste0("indel_", new_labels)
	
	annot_table[which_to_relabel, "impact"] = new_labels

	# removing silent mutations here	
	annot_table = annot_table[annot_table$impact != "Synonymous",]
	
	annot_table
}


gsub_TCGA_sample_patient = function(x) {
	gsub("^(TCGA)-(.+)-(.+)-(.+)+", "\\1-\\2-\\3", x)
}
