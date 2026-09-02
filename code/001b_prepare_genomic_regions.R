source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(ActiveDriverWGS)
library(rtracklayer)
library(GenomicFeatures)

# for non-coding elements once we remove protein-coding exons we also remove very small fragments between CDS
MIN_SIZE_POST_SUBTRACT = 50

# function to load gencode file and keep only the relevant info 
clean_gencode_file = function(gencode_fname, element_type, permitted_gene_types) {
	
	gencode = read.delim(gencode_fname, 
		skip = 5, header = FALSE, stringsAsFactors = FALSE)
	
	# TSS +/- 1000	
	if (element_type == "TSS") {
		gencode = gencode[gencode$V3 == "transcript", ]	
		tss_mid = sapply(1:nrow(gencode), function(i) ifelse(gencode[i, "V7"]	=="+", gencode[i, "V4"], gencode[i, "V5"]) )
		gencode$V4 = tss_mid - 1000
		gencode$V5 = tss_mid + 999
	} else {
		gencode = gencode[gencode$V3 == element_type, ]
	}

	annots = strsplit(gencode$V9, split = ";")
	gene_id = gsub("^gene_id ", "", sapply(annots, '[[', 1))
	gene_name = gsub("^ gene_name ", "", sapply(annots, '[[', 5))
	gene_type = gsub("^ gene_type ", "", sapply(annots, '[[', 3))
	gene_status = gsub("^ gene_status ", "", sapply(annots, '[[', 4))
	combined_id = paste0("gencode::", gene_name, "::", gene_id)
	
	gencode_genes = data.frame(chr = gencode$V1, start = gencode$V4, end = gencode$V5, 
			id = combined_id, gene_type, gene_status, gene_strand = gencode$V7,
			stringsAsFactors = FALSE)
			
	gencode_genes = gencode_genes[gencode_genes$gene_type %in% permitted_gene_types, ]
}


# id = "gencode::A4GALT::ENSG00000128274.11"; split_prepared_list = split_prepared_CDS
make_reduced_gr = function(id, split_prepared_list) {
	
	if (runif(1) < 0.001) cat("o")
	
	this_coords = split_prepared_list[[id]]
	
	gr_coords = reduce(makeGRangesFromDataFrame(this_coords))
	gr_coords$id = id
	gr_coords	
}




# remove CDS from the noncoding sites
# reassign names of encode coords
subtract_CDS_by_gene = function(id, gr_split_elements, gr_CDS_reduced) {
	
	if (runif(1)<0.001) cat("o")
	
	gr_this_el = gr_split_elements[[id]]
	gr_this_el = setdiff(gr_this_el, gr_CDS_reduced, ignore.strand = TRUE)
	if (length(gr_this_el)==0) {
		return(NULL)
	}
	gr_this_el$id = id
	gr_this_el
}

# gr_elements = gr_chrom_sites [gr_chrom_sites$id == "CHRMM::MISC::chr10_80812000_80812199::PrimEnh_Prost"]
# gr_elements = gr_UTR; min_width = 50
subtract_CDS = function(gr_elements, gr_CDS, min_width = 100) {

	# separate logic: if elements have overlapping parts we need to subtract one by one, split by element ID
	gr_CDS_reduced = reduce(gr_CDS)
	if (length(gr_elements) != length(reduce(gr_elements))) {
		
		gr_split_elements = split(gr_elements, gr_elements$id)
		cat(length(gr_split_elements), ": ")
		gr_elements_minus_CDS = mclapply(names(gr_split_elements), subtract_CDS_by_gene, gr_split_elements, gr_CDS_reduced, mc.cores= 8)
		gr_elements_minus_CDS = do.call(c, gr_elements_minus_CDS)
		
	} else {
		
		gr_elements_minus_CDS = setdiff(gr_elements, gr_CDS_reduced, ignore.strand = TRUE)
			
		## names of CRMs to be added back to the setDiff
		ov = findOverlaps(gr_elements, gr_elements_minus_CDS)
	
		# add IDs back based on overlap
		mcols(gr_elements_minus_CDS)[, "id"] = NA
		gr_elements_minus_CDS[subjectHits(ov)]$id = gr_elements[queryHits(ov)]$id
	}

	gr_elements_minus_CDS = gr_elements_minus_CDS[width(gr_elements_minus_CDS) >= min_width]
	gr_elements_minus_CDS
}






#
# prepare genomic elements object for parallel running of ADWGS
#
# CDS first

gencode_fname = "elements/gencode.v19.annotation.gtf.gz"
prepared_CDS = clean_gencode_file(gencode_fname, "CDS", "protein_coding")
prepared_CDS = prepared_CDS[prepared_CDS$gene_status == "KNOWN",]
prepared_CDS = unique(prepared_CDS)

## curate and remove duplicated symbols
sym = gsub("(.+)::(.+)::(.+)", "\\2", prepared_CDS$id)
ens = gsub("(.+)::(.+)::(.+)", "\\3", prepared_CDS$id)
dup_sym_dfr = data.frame(sym, ens, stringsAsFactors = FALSE)
dup_sym_dfr = unique(dup_sym_dfr)

dup_symbols = unique(dup_sym_dfr$sym[duplicated(dup_sym_dfr$sym)])
dup_sym_dfr = dup_sym_dfr[dup_sym_dfr$sym %in% dup_symbols,]
dup_sym_dfr = dup_sym_dfr[order(dup_sym_dfr$sym),]
dup_sym_dfr$to_remove = FALSE

# review duplicated symbols from file, save new file and load it in again
fname = pff("figures/dup_sym_dfr.tsv")
write.table(dup_sym_dfr, file = fname, sep = "\t", row.names = FALSE)
file_open_call2(fname)

# read in the curated list of ENSGs to keep and to remove
dup_sym_dfr_EDIT = read.delim(pff("figures/dup_sym_dfr_EDIT.tsv"), stringsAsFactors = FALSE)
ENSG_to_remove = unique(dup_sym_dfr_EDIT[dup_sym_dfr_EDIT$to_remove == TRUE, "ens"])
prepared_CDS$ENSG = gsub("(.+)::(.+)::(.+)", "\\3", prepared_CDS$id)
prepared_CDS = prepared_CDS[!prepared_CDS$ENSG %in% ENSG_to_remove,]
save(prepared_CDS, file = pff("prepared_CDS.rsav"))


#
# merge protein-coding exons
#
split_prepared_CDS = split(prepared_CDS, prepared_CDS$id)
gr_CDS = do.call(c, mclapply(names(split_prepared_CDS), make_reduced_gr, split_prepared_CDS, mc.cores = 8))
gr_CDS$id2 = gr_CDS$id
gr_CDS$id = paste0("CDS::", gr_CDS$id)
gr_CDS$symbol = gsub("(.+)::(.+)::(.+)::(.+)", "\\3", gr_CDS$id)

# use this list to filter permitted genes later
CDS_ID2_list = unique(gr_CDS$id2)
save(gr_CDS, file = pff("gr_CDS.rsav"))



#
# UTRs: 5' and 3' separately
#

# UTR_input = utr3; utr_tag = "UTR5"
prepare_UTR_by_prime = function(UTR_input, tx2gene, ensg_to_geneid, gr_CDS, utr_tag, MIN_SIZE_POST_SUBTRACT) {

	# map internal gene coords to ENSG coords
	all_UTR = unlist(UTR_input, use.names = FALSE)
	
	mcols(all_UTR)$tx_id = rep(as.integer(names(UTR_input)), lengths(UTR_input))
		
	match_UTR = match(mcols(all_UTR)$tx_id, tx2gene$tx_id)
	mcols(all_UTR) = cbind(mcols(all_UTR), tx2gene[match_UTR, -1L, drop = FALSE])
		
	all_UTR = all_UTR[all_UTR$gene_id %in% names(ensg_to_geneid)]
	all_UTR$id = ensg_to_geneid[all_UTR$gene_id]

	split_prepared_UTR = split(all_UTR, mcols(all_UTR)$id)
		
	# merge UTRs by gene
	gr_UTR = do.call(c, mclapply(names(split_prepared_UTR), make_reduced_gr, split_prepared_UTR, mc.cores = 8))
	gr_UTR$id2 = gr_UTR$id
	gr_UTR$id = paste0(utr_tag, "::", gr_UTR$id)
	gr_UTR$symbol = gsub("(.+)::(.+)::(.+)::(.+)", "\\3", gr_UTR$id)
	
	# remove CDS from UTR just in case	
	gr_UTR = subtract_CDS(gr_UTR, gr_CDS, MIN_SIZE_POST_SUBTRACT)
	gr_UTR
}


## import GTF
gtf = makeTxDbFromGFF(gencode_fname)
gtf = keepStandardChromosomes(gtf, pruning = 'coarse')
tx2gene = mcols(transcripts(gtf, columns=c("tx_id", "tx_name", "gene_id")))
tx2gene$gene_id = as.character(tx2gene$gene_id)
# keep only CDS, map to long::gene::id using this data structure
ensg_to_geneid = structure(unique(gr_CDS$id2), names = gsub("(.+)::(.+)::(.+)", "\\3", unique(gr_CDS$id2)))

## get 5' and 3' utrs separately
utr3 = threeUTRsByTranscript(gtf)
utr5 = fiveUTRsByTranscript(gtf)

gr_UTR3 = prepare_UTR_by_prime(utr3, tx2gene, ensg_to_geneid, gr_CDS, "UTR3", MIN_SIZE_POST_SUBTRACT)
gr_UTR5 = prepare_UTR_by_prime(utr5, tx2gene, ensg_to_geneid, gr_CDS, "UTR5", MIN_SIZE_POST_SUBTRACT)

save(gr_UTR3, file = pff("gr_UTR3.rsav"))
save(gr_UTR5, file = pff("gr_UTR5.rsav"))


#
# take lncRNAs from Gencode; remove duplicated symbols; subtract CDS sequences from lncRNAs and remove RNAs that are too small after the subtract
#
lncRNA_types = c("lincRNA", "antisense", "miRNA", "sense_intronic")
prepared_lncRNA = clean_gencode_file(gencode_fname, "gene", lncRNA_types)

## remove duplicated symbols
all_LNC_ids = unique(prepared_lncRNA$id)
all_LNC_symbols = gsub("(gencode)::(.+)::(.+)", "\\2", all_LNC_ids)
duplicated_symbols = all_LNC_symbols[which(duplicated(all_LNC_symbols))]
ids_of_duplicated_LNC_symbols = all_LNC_ids[all_LNC_symbols %in% duplicated_symbols]
prepared_lncRNA = prepared_lncRNA[!prepared_lncRNA$id %in% ids_of_duplicated_LNC_symbols,]
save(prepared_lncRNA, file = pff("prepared_lncRNA.rsav"))

gr_lncRNA = GRanges(prepared_lncRNA$chr,
			                 IRanges(start = prepared_lncRNA$start, end = prepared_lncRNA$end),
			                 id = paste0("lncRNA::", prepared_lncRNA$id))
gr_lncRNA$id2 = gsub("^lncRNA::", "", prepared_lncRNA$id)
gr_lncRNA$symbol = gsub("(.+)::(.+)::(.+)::(.+)", "\\3", gr_lncRNA$id)
gr_lncRNA$gene_type = prepared_lncRNA$gene_type


# remove lncRNA parts overlapping protein-coding exons of other genes
gr_lncRNA = subtract_CDS(gr_lncRNA, gr_CDS, MIN_SIZE_POST_SUBTRACT)
gr_lncRNA$id2 = gsub("lncRNA::", "", gr_lncRNA$id)
gr_lncRNA$symbol = gsub("(.+)::(.+)::(.+)", "\\2", gr_lncRNA$id)
# use this list to filter permitted genes later
LNCRNA_ID2_list = unique(gr_lncRNA$id2)
save(gr_lncRNA, file = pff("gr_lncRNA.rsav"))



##
## GENCODE gene coordinates - use those to annotate non-coding elements later
###

# all genes need to be either protein-coding or lncRNA genes we selected previously
previously_selected_ids = unique(c(CDS_ID2_list, LNCRNA_ID2_list))

# select cds and nc genes of interest
gencode_genes = clean_gencode_file(gencode_fname, "gene", c("protein_coding", "lincRNA", "antisense", "miRNA", "sense_intronic"))
gencode_genes = gencode_genes[gencode_genes$id %in% previously_selected_ids,]
gencode_genes$id2 = gencode_genes$id
gencode_genes$id = paste0("gene::", gencode_genes$id)

gr_gene_coords = GRanges(gencode_genes$chr,
			                 IRanges(start = gencode_genes$start, end = gencode_genes$end),
			                 id = gencode_genes$id, 
			                 id2 = gencode_genes$id2,
			                 strand = "*",
			                 gene_strand = gencode_genes$gene_strand, 
			                 symbol = gsub("(gencode)::(.+)::(.+)", "\\2", gencode_genes$id2))

save(gr_gene_coords, file = pff("gr_gene_coords.rsav"))
save(gencode_genes, file = pff("gencode_genes.rsav"))




##
# select cdsGenes for SV analysis
##
gr_CDSgene_coords = gr_gene_coords[gr_gene_coords$id2 %in% CDS_ID2_list]
gencode_CDSgenes = gencode_genes[gencode_genes$id2 %in% CDS_ID2_list,]

gr_CDSgene_coords$id = gsub("^gene::", "CDSgene::", gr_CDSgene_coords$id)
gencode_CDSgenes$id = gsub("^gene::", "CDSgene::", gencode_CDSgenes$id)

save(gencode_CDSgenes, file = pff("gencode_CDSgenes.rsav"))
save(gr_CDSgene_coords, file = pff("gr_CDSgene_coords.rsav"))



##
# create transcription start site dataset to annotate loops and genes
##
prepared_TSS = clean_gencode_file(gencode_fname, "TSS", c("protein_coding", "lincRNA", "antisense", "miRNA", "sense_intronic"))
prepared_TSS = prepared_TSS[prepared_TSS$id %in% previously_selected_ids,]
prepared_TSS = unique(prepared_TSS)
save(prepared_TSS, file = pff("prepared_TSS.rsav"))

split_prepared_TSS = split(prepared_TSS, prepared_TSS$id)
gr_TSS = do.call(c, mclapply(names(split_prepared_TSS), make_reduced_gr, split_prepared_TSS, mc.cores = 8))
gr_TSS$id2 = gr_TSS$id
gr_TSS$id = paste0("TSS::", gr_TSS$id)
gr_TSS$symbol = gsub("(.+)::(.+)::(.+)::(.+)", "\\3", gr_TSS$id)
save(gr_TSS, file = pff("gr_TSS.rsav"))


##
# ATAC seq in hg38, need to convert to hg19
##
ATAC38 = read.delim(
		paste0("DATA_USED__", this_timestamp, "/genomic_elements_Kevin_2022-09-27/PRAD_peakCalls.txt"), 
		stringsAsFactors = FALSE)
gr_ATAC38 = GRanges(ATAC38[,1], IRanges(ATAC38[,2], ATAC38[,3]))

chain_hg38 = import.chain(
		paste0("DATA_USED__", this_timestamp, "/genomic_elements_Kevin_2022-09-27/hg38ToHg19.over.chain"))
gr_ATAC = liftOver(gr_ATAC38, chain_hg38)
gr_ATAC = gr_ATAC[which(sapply(start(gr_ATAC), length) == 1)]
gr_ATAC = unlist(gr_ATAC)
gr_ATAC = reduce(gr_ATAC)
gr_ATAC$id = paste0("ATAC::TCGA::", seqnames(gr_ATAC), "_", start(gr_ATAC), "_", end(gr_ATAC), "::ATAC")

# small atac peaks excluded; CDS excluded
gr_ATAC = gr_ATAC[width(gr_ATAC) >= 100]
gr_ATAC = subtract_CDS(gr_ATAC, gr_CDS, MIN_SIZE_POST_SUBTRACT)


##
# normal AR binding sites
##
narbs_sites = prepare_elements_from_BED4(
		paste0("DATA_USED__", this_timestamp, "/genomic_elements_Kevin_2022-09-27/NARBS_with_ID.bed"))
gr_narbs_sites = GRanges(narbs_sites[,1], IRanges(narbs_sites[,2], narbs_sites[,3]))
# make sure the sites are unique
gr_narbs_sites = reduce(gr_narbs_sites)

gr_narbs_sites$id = 
		paste0("NARBS::MISC::", seqnames(gr_narbs_sites), "_", 
			start(gr_narbs_sites), "_", end(gr_narbs_sites), "::NARBS")
# CDS excluded
gr_narbs_sites = subtract_CDS(gr_narbs_sites, gr_CDS, min_width = MIN_SIZE_POST_SUBTRACT)


##
#tumor AR binding sites
##
tarbs_sites = prepare_elements_from_BED4(
		paste0("DATA_USED__", this_timestamp, "/genomic_elements_Kevin_2022-09-27/TARBS_with_ID.bed"))
gr_tarbs_sites = GRanges(tarbs_sites[,1], IRanges(tarbs_sites[,2], tarbs_sites[,3]))
gr_tarbs_sites = reduce(gr_tarbs_sites)

gr_tarbs_sites$id = 
		paste0("TARBS::MISC::", seqnames(gr_tarbs_sites), "_", 
			start(gr_tarbs_sites), "_", end(gr_tarbs_sites), "::TARBS")

# CDS excluded
gr_tarbs_sites = subtract_CDS(gr_tarbs_sites, gr_CDS, min_width = MIN_SIZE_POST_SUBTRACT)


##
# original ChromHMM tracks
##
# long elements excluded
MAX_ELEMENT_LENGTH = 25000

chrom_sites = read.delim(
		paste0("DATA_USED__", this_timestamp, "/pomerantz_chromatin_states_2024-04-16/combined_10_states_colour_final_mod.bed"),
		stringsAsFactors = FALSE, header = FALSE, sep = "\t", skip = 1)
chrom_sites = chrom_sites[,1:4]
colnames(chrom_sites) = c("chr", "start", "end", "element_type")

# reduce end coordinate to make sure neighbours stop overlapping
chrom_sites$end = chrom_sites$end - 2

grgr = GRanges(chrom_sites$chr, IRanges(chrom_sites$start, chrom_sites$end))

# this is the list of elements we consider - enhancers, promoters, etc.
element_types = read.delim(
		paste0("DATA_USED__", this_timestamp, "/pomerantz_chromatin_states_2024-04-16/Pomerantz_states_names.txt"),
		stringsAsFactors = FALSE, sep = "\t", header = FALSE)
colnames(element_types) = c("short_name", "long_name")

chrom_sites = merge(chrom_sites, element_types, by.x = "element_type", by.y = "long_name")
chrom_sites = chrom_sites[,-1]
colnames(chrom_sites) = c("chr", "start", "end", "element_type")


states_exclude = c("HetChrom", "RepChrom")
chrom_sites = chrom_sites[!chrom_sites$element_type %in% states_exclude,] 
chrom_sites$id = 
		paste0("CHRMM::MISC::", chrom_sites$chr, "_", 
			chrom_sites$start, "_", chrom_sites$end, "::", chrom_sites$element_type)
chrom_sites = chrom_sites[, c("chr", "start", "end", "id")]

# chrom HMM tracks
gr_chrom_sites = GRanges(chrom_sites[,1], IRanges(chrom_sites[,2], chrom_sites[,3]), id = chrom_sites$id)
gr_chrom_sites = gr_chrom_sites[width(gr_chrom_sites) < MAX_ELEMENT_LENGTH]

# exclude overlaps with CDS
gr_chrom_sites = subtract_CDS(gr_chrom_sites, gr_CDS, min_width = MIN_SIZE_POST_SUBTRACT)


#
# chromatin loops - both for linking elements and as driver search space
# 
loops = read.delim(
		paste0("DATA_USED__", this_timestamp, "/prostate_HiC_Joachim_2022-09-29/20220929_GSE164347_loops_hg19/GSE164347.loops.merge.hg19.bedpe"), 
		stringsAsFactors = FALSE, header = FALSE)
colnames(loops) = c("chr1", "start1", "end1", "chr2", "start2", "end2", "id_combined", "id1", "id2", 
		"score1", "score2", "score3", "score4", "samples")

loops = loops[which(apply(loops, 1, function(x) !any(is.na(x)))),]
loops = loops[loops$chr1 == loops$chr2,]

loops$mid1 = (loops$start1 + loops$end1 ) / 2
loops$mid2 = (loops$start2 + loops$end2 ) / 2
loops$tag = paste(loops$chr1, loops$start1, loops$end1, loops$chr2, loops$start2, loops$end2, sep = "__")

save(loops, file = pff("loops.rsav"))


gr_anchors = c(
		GRanges(loops$chr1, IRanges(loops$start1, loops$end1)),
		GRanges(loops$chr2, IRanges(loops$start2, loops$end2))
)
gr_anchors = reduce(gr_anchors)

gr_anchors$id = paste0("LOOPS::NA::", seqnames(gr_anchors), ":", start(gr_anchors), "-", end(gr_anchors), "::LOOPS")
gr_anchors = subtract_CDS(gr_anchors, gr_CDS, min_width = MIN_SIZE_POST_SUBTRACT)

#
# put all elements all together into a single object
#
gr_prepared_elements = c(gr_CDS, gr_UTR5, gr_UTR3, gr_ATAC, gr_lncRNA, gr_narbs_sites, gr_tarbs_sites, gr_chrom_sites, gr_anchors, gr_CDSgene_coords)
save(gr_prepared_elements, file = pff("gr_prepared_elements.rsav"))
