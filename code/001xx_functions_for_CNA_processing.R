# remove 1kbps sections
MIN_SEGMENT_SIZE = 1000
# what percentage of gene sequence needs to be AMP to call the gene AMP
FRACTION_GENE_ALTERED = 0.5

# discretise copy number calls into loss, gain, full deletion, high gain. 
# separately for diploid and WGD; separately for autosomes and X
CNA_cutoffs = list()
CNA_cutoffs[['gain_diploid_cutoff']] = 3
CNA_cutoffs[['gain_WGD_cutoff']] = 6
CNA_cutoffs[['loss_diploid_cutoff']] = 1
CNA_cutoffs[['loss_WGD_cutoff']] = 1
CNA_cutoffs[['high_gain_diploid_cutoff']] = 5
CNA_cutoffs[['high_gain_WGD_cutoff']] = 8
CNA_cutoffs[['gain_sex_diploid_cutoff']] = 3 
CNA_cutoffs[['gain_sex_WGD_cutoff']] = 4
CNA_cutoffs[['loss_sex_diploid_cutoff']] = 0
CNA_cutoffs[['loss_sex_WGD_cutoff']] = 0
CNA_cutoffs[['high_gain_sex_diploid_cutoff']] = 3 
CNA_cutoffs[['high_gain_sex_WGD_cutoff']] = 4


#
# read in all CNA segments
#
# fname = "DATA_USED__2024-04-16/SCNA_with_BD_SVs_03_05_2023/Subclonal_SCNA_with_Avg_CN//PPCG0001a_DNA_vs_PPCG0001b_DNA_Subclonal_SCNA.txt"
get_CNA_calls = function(fname, WGD_tumors, CNA_cutoffs, MIN_SEGMENT_SIZE) {

	cat(".")
	gain_diploid_cutoff = CNA_cutoffs[['gain_diploid_cutoff']]
	gain_WGD_cutoff = CNA_cutoffs[['gain_WGD_cutoff']]
	loss_diploid_cutoff = CNA_cutoffs[['loss_diploid_cutoff']]
	loss_WGD_cutoff = CNA_cutoffs[['loss_WGD_cutoff']]
	high_gain_diploid_cutoff = CNA_cutoffs[['high_gain_diploid_cutoff']]
	high_gain_WGD_cutoff = CNA_cutoffs[['high_gain_WGD_cutoff']]
	gain_sex_diploid_cutoff = CNA_cutoffs[['gain_sex_diploid_cutoff']]
	gain_sex_WGD_cutoff = CNA_cutoffs[['gain_sex_WGD_cutoff']]
	loss_sex_diploid_cutoff = CNA_cutoffs[['loss_sex_diploid_cutoff']]
	loss_sex_WGD_cutoff = CNA_cutoffs[['loss_sex_WGD_cutoff']]
	high_gain_sex_diploid_cutoff = CNA_cutoffs[['high_gain_sex_diploid_cutoff']]
	high_gain_sex_WGD_cutoff = CNA_cutoffs[['high_gain_sex_WGD_cutoff']]
	
	dat = read.delim(fname, stringsAsFactors = FALSE)
	sample_id = gsub("(.+\\/)+(.+)_vs_(.+)", "\\2", fname)
	is_WGD = sample_id %in% WGD_tumors
	dat$WGD = is_WGD
	
	dat$total_cn = dat$nMaj1_A + dat$nMin1_A
	dat$tumor_id = sample_id
	dat = dat[, c("chr", "startpos", "endpos", "total_cn", "WGD", "tumor_id")]
	dat$chr = paste0("chr", dat$chr)
	dat = dat[!is.na(dat$total_cn), ]
	dat$cn_width = dat$endpos - dat$startpos
	
	dat = dat[dat$cn_width >= MIN_SEGMENT_SIZE,]

	# relative copy number - subtract total copy number	
	genomewide_cn_autosome = ifelse(is_WGD, 4, 2)
	genomewide_cn_sex = ifelse(is_WGD, 2, 1)
	dat$rel_cn = dat$total_cn - genomewide_cn_autosome
	dat[dat$chr %in% c("chrX", "chrY"), "rel_cn"] = dat[dat$chr %in% c("chrX", "chrY"), "total_cn"] - genomewide_cn_sex
	
	# fold-change copy number - divide total copy number by expected copy number; for gistic
	dat$foldchange_cn = dat$total_cn / genomewide_cn_autosome
	dat[dat$chr %in% c("chrX", "chrY"), "foldchange_cn"] = dat[dat$chr %in% c("chrX", "chrY"), "total_cn"] / genomewide_cn_sex

	# what value to use to annotate gains and losses
	gain_cutoff = ifelse (is_WGD, gain_WGD_cutoff, gain_diploid_cutoff)
	loss_cutoff = ifelse (is_WGD, loss_WGD_cutoff, loss_diploid_cutoff)
	chrsex_gain_cutoff = ifelse (is_WGD, gain_sex_WGD_cutoff, gain_sex_diploid_cutoff)
	chrsex_loss_cutoff = ifelse (is_WGD, loss_sex_WGD_cutoff, loss_sex_diploid_cutoff)
	high_gain_cutoff = ifelse(is_WGD, high_gain_WGD_cutoff, high_gain_diploid_cutoff)
	high_chrsex_gain_cutoff = ifelse(is_WGD, high_gain_sex_WGD_cutoff, high_gain_sex_diploid_cutoff)
	
	# add interesting annotations
	dat$annot = "bal"
	dat$CNA_id = paste(dat$chr, dat$startpos, dat$endpos, dat$tumor_id, sep = ":")
	
	# annotate segments into four categories separately for autosomes and sex chromosomes
	chr_auto = paste0("chr", 1:22)
	chr_sex = paste0("chr", c("X", "Y"))
	dat[dat$chr %in% chr_auto & dat$total_cn <= loss_cutoff, "annot"] = "loss"
	dat[dat$chr %in% chr_auto & dat$total_cn >= gain_cutoff, "annot"] = "gain"
	dat[dat$chr %in% chr_auto & dat$total_cn >= high_gain_cutoff, "annot"] = "high_gain"
	dat[dat$chr %in% chr_sex & dat$total_cn <= chrsex_loss_cutoff, "annot"] = "loss"
	dat[dat$chr %in% chr_sex & dat$total_cn >= chrsex_gain_cutoff, "annot"] = "gain"
	dat[dat$chr %in% chr_sex & dat$total_cn >= high_chrsex_gain_cutoff, "annot"] = "high_gain"
	dat[dat$total_cn == 0, "annot"] = "full_loss"
	
	dat
}


# compute percent genome altered
compute_PGA = function(patient, prepared_CNAs) {
	cat(".")
	cnas_here = prepared_CNAs[prepared_CNAs$patient == patient,]
	# autosomes chromosomes only for PGA
	cnas_here = cnas_here[cnas_here$chr %in% paste0("chr", 1:22),]

	PGA = as.double(sum(cnas_here[cnas_here$annot != "bal", "cn_width"], na.rm = TRUE) / sum(cnas_here$cn_width, na.rm = TRUE))
	
	n_small = sum(cnas_here$cn_width < 1000 & cnas_here$annot != "bal")
	n_total = sum(cnas_here$annot != "bal")
	pct_small = n_small / n_total
	
	data.frame(patient, PGA, n_small, n_total, pct_small, stringsAsFactors = FALSE)
}


# percent overlaps of genes and CNAs per patient
find_gene_percent_cna_coverage = function(gene, gr_CNAs_here, gr_genes_with_CNAs) {
	
	if (runif(1) < 0.001) cat("o")

	gr_gene = gr_genes_with_CNAs[gr_genes_with_CNAs$symbol == gene]
		
	gr_intersecting = pintersect(gr_CNAs_here, gr_gene)
	gr_intersecting1 = gr_intersecting[gr_intersecting$hit]
	
	gr_intersecting1$seg_width = width(gr_intersecting1)
	gr_intersecting1$percent_affected = width(gr_intersecting1)/width(gr_gene)
	
	patient_percentages = c(by(gr_intersecting1$percent_affected, gr_intersecting1$patient, sum))
	patient_weighted_rel_cn = c(by(gr_intersecting1$percent_affected * gr_intersecting1$rel_cn, gr_intersecting1$patient, sum))
	this_patients = names(patient_percentages)
	
	data.frame(gene, patient = this_patients, 
			pct_gene_affected = patient_percentages[this_patients], 
			rel_cn = patient_weighted_rel_cn[this_patients], 
			stringsAsFactors = FALSE)
}


gene_cna_coverage = function(annots, prepared_CNAs, gr_genes, FRACTION_GENE_ALTERED, tag) {

	CNAs_here = prepared_CNAs[prepared_CNAs$annot %in% annots,]
	gr_CNAs_here = GRanges(CNAs_here$chr, IRanges(CNAs_here$startpos, CNAs_here$endpos), 
			annot = CNAs_here$annot, patient = CNAs_here$patient, rel_cn = CNAs_here$rel_cn)
	gr_genes_with_CNAs = gr_genes[countOverlaps(gr_genes, gr_CNAs_here) > 0]

	gene_CNA_overlaps = do.call(rbind, mclapply(gr_genes_with_CNAs$symbol, find_gene_percent_cna_coverage, gr_CNAs_here, gr_genes_with_CNAs, 
			mc.cores = 6))
	
	gene_CNA_overlaps$annots = paste(annots, collapse = ",")
	gene_CNA_overlaps$tag = tag
	gene_CNA_overlaps$fraction_gene_altered = FRACTION_GENE_ALTERED

	gene_CNA_overlaps = gene_CNA_overlaps[gene_CNA_overlaps$pct_gene_affected >= FRACTION_GENE_ALTERED,]
	gene_CNA_overlaps
}
