source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")



##
# annotate non-coding elements with target genes
##

flank_window = 3000

load(file = pff("gr_prepared_elements.rsav"))
gr_cds = gr_prepared_elements[grep("^CDS::", gr_prepared_elements$id)]
gr_lnc = gr_prepared_elements[grep("^lncRNA::", gr_prepared_elements$id)]
gr_cdsGene = gr_prepared_elements[grep("^CDSgene::", gr_prepared_elements$id)]

gr_genes = c(gr_lnc, gr_cdsGene)
gr_noncoding = gr_prepared_elements[grep("^CDS::|^lncRNA::|^CDSgene::", gr_prepared_elements$id, invert = TRUE)]


create_noncoding_annotations = function(gr_genes, gr_noncoding, bps_flank = 0) {

	gr_genes_win = gr_genes
	
	start(gr_genes_win) = start(gr_genes_win) - bps_flank
	end(gr_genes_win) = end(gr_genes_win) + bps_flank
	
	ov = findOverlaps(gr_genes_win, gr_noncoding)
	
	gene_noncoding_map = unique(data.frame(
			nc_id = gr_noncoding[subjectHits(ov)]$id, 
			gene_id = gr_genes_win[queryHits(ov)]$id,
			stringsAsFactors = FALSE
	))
	gene_noncoding_map
}



annot_flank_orig = create_noncoding_annotations(gr_genes, gr_noncoding, flank_window)
save(annot_flank_orig, file = pff("annot_flank_orig.rsav"))

# add 1-1 mapping between CDS regions for completeness
annot_cds = unique(cbind(nc_id = gr_cds$id, gene_id = gr_cds$id))
annot_lnc = unique(cbind(nc_id = gr_lnc$id, gene_id = gr_lnc$id))
annot_cdsGene = unique(cbind(nc_id = gr_cdsGene$id, gene_id = gr_cdsGene$id))

# each CDS element flanks the same elements as CDSgene
annot_flank_orig1 = annot_flank_orig[, c("nc_id", "gene_id")]
annot_flank_orig2 = annot_flank_orig1[grep("^CDSgene::", annot_flank_orig1$gene),]
annot_flank_orig2$gene_id = gsub("^CDSgene::", "CDS::", annot_flank_orig2$gene_id)

annot_flank = rbind(annot_flank_orig1, annot_flank_orig2, annot_cds, annot_lnc, annot_cdsGene)
rm(annot_flank_orig, annot_flank_orig1, annot_flank_orig2)


annot_flank$gene2element = paste(annot_flank$gene_id, annot_flank$nc_id, sep = "__")
save(annot_flank, file = pff("annot_flank.rsav"))


gene_to_nc_flank = split(annot_flank[,"nc_id"], annot_flank[,"gene_id"])
nc_to_gene_flank = split(annot_flank[,"gene_id"], annot_flank[,"nc_id"])

save(gene_to_nc_flank, file = pff("gene_to_nc_flank.rsav"))
save(nc_to_gene_flank, file = pff("nc_to_gene_flank.rsav"))



###
## loops to elements to genes
###
map_elements_via_loops = function(gr_anchor1, gr_anchor2, gr_nc_elements, gr_TSS) {

	# first, get subset of anchors that have an overlap with distal elements
	gr_elements_combined = c(gr_nc_elements, gr_TSS)
	
	ov_anchor1_element = findOverlaps(gr_anchor1, gr_elements_combined)
	anchor1_to_element = data.frame(
			loop_id = gr_anchor1$id[queryHits(ov_anchor1_element)],
			element1_id = gr_elements_combined$id[subjectHits(ov_anchor1_element)], 
			stringsAsFactors = FALSE)
	
	ov_anchor2_element = findOverlaps(gr_anchor2, gr_elements_combined)
	anchor2_to_element = data.frame(
			loop_id = gr_anchor2$id[queryHits(ov_anchor2_element)],
			element2_id = gr_elements_combined$id[subjectHits(ov_anchor2_element)], 
			stringsAsFactors = FALSE)

	# one end of loop needs to be TSS, and the other end a non-coding element			
	dfr = merge(anchor1_to_element, anchor2_to_element, by = "loop_id")		
	dfr1 = dfr[grep("^TSS::", dfr$element1_id), c("element1_id", "element2_id", "loop_id"),]
	dfr2 = dfr[grep("^TSS::", dfr$element2_id), c("element2_id", "element1_id", "loop_id"),]
	colnames(dfr1) = colnames(dfr2) = c("gene_id", "nc_id", "loop_id")

	# keep loops that connect TSS with other elements
	dfr = unique(rbind(dfr1, dfr2))
	
	# exclude self-loops
	gene1 = gsub("(.+)::(.+)::(.+)::(.+)", "\\2::\\3::\\4", dfr$gene_id)
	gene2 = gsub("(.+)::(.+)::(.+)::(.+)", "\\2::\\3::\\4", dfr$nc_id)
	dfr = dfr[gene1 != gene2,]
	
	element_to_loop_to_element = dfr
	
	colnames(anchor1_to_element) = colnames(anchor2_to_element) = c("loop_id", "element_id")
	anchor1_to_element$anchor = "1"
	anchor2_to_element$anchor = "2"
	anchor_to_element = rbind(anchor1_to_element, anchor2_to_element)
	
	list(element_to_loop_to_element, anchor_to_element)
}

#
# loops initially assigned to TSSs; take these and reassign them to coding or lncRNA genes; 
# coding genes twice: for CDS for SNV analysis and for gene bodies for SV analysis
#

# dfr_element_loops = dfr_loops_element_to_TSS
assign_genes_to_element_loops = function(dfr_element_loops, gr_lnc, gr_cdsGene) {
	
	lnc_ids = gsub("^lncRNA::", "", gr_lnc$id)
	cds_ids = gsub("^CDSgene::", "", gr_cdsGene$id)
	
	# replace gene ids with CDSgene or lncRNA
	gene_ids = gsub("^TSS::", "", dfr_element_loops$gene_id)
	gene_ids[gene_ids %in% lnc_ids] = paste("lncRNA::", gene_ids[gene_ids %in% lnc_ids], sep = "")
	gene_ids[gene_ids %in% cds_ids] = paste("CDSgene::", gene_ids[gene_ids %in% cds_ids], sep = "")
	dfr_element_loops$gene_id = gene_ids

	# replace non-coding element ids with CDSgene or lncRNA	
	nc_ids = gsub("^TSS::", "", dfr_element_loops$nc_id)
	nc_ids[nc_ids %in% lnc_ids] = paste("lncRNA::", nc_ids[nc_ids %in% lnc_ids], sep = "")
	nc_ids[nc_ids %in% cds_ids] = paste("CDSgene::", nc_ids[nc_ids %in% cds_ids], sep = "")
	dfr_element_loops$nc_id = nc_ids
	
	# both gene IDs and nc IDs will be doubled - CDS and CDSgene have it
	dfr_element_loops_CDS = dfr_element_loops[grep("^CDSgene::", dfr_element_loops$gene_id),]
	dfr_element_loops_CDS$gene_id = gsub("^CDSgene::", "CDS::", dfr_element_loops_CDS$gene_id)
	
	dfr_element_loops_CDS2 = dfr_element_loops[grep("^CDSgene::", dfr_element_loops$nc_id),]
	dfr_element_loops_CDS2$nc_id = gsub("^CDSgene::", "CDS::", dfr_element_loops_CDS2$nc_id)
	
	dfr_element_loops = rbind(dfr_element_loops, dfr_element_loops_CDS, dfr_element_loops_CDS2)
	
	# tag each interaction to filter later vs locals
	dfr_element_loops$gene2element = paste(dfr_element_loops$gene_id, dfr_element_loops$nc_id, sep = "__")
	dfr_element_loops	
}


load(file = pff("gr_prepared_elements.rsav"))
load(file = pff("loops.rsav"))
load(file = pff("gr_TSS.rsav"))
load(file = pff("gene_to_nc_flank.rsav"))

gr_anchor1 = GRanges(loops$chr1, IRanges(loops$start1, loops$end1), id = loops$tag)
gr_anchor2 = GRanges(loops$chr2, IRanges(loops$start2, loops$end2), id = loops$tag)

# same noncoding elements as for flanking annotations
list_elements_via_loops = map_elements_via_loops(gr_anchor1, gr_anchor2, gr_noncoding, gr_TSS)
dfr_loops_element_to_TSS = list_elements_via_loops[[1]]
dfr_loops_to_element = list_elements_via_loops[[2]]
# convert TSSs to either CDS, CDSgene, or lncRNA
dfr_loops_element_to_gene = assign_genes_to_element_loops(dfr_loops_element_to_TSS, gr_lnc, gr_cdsGene)
# remove direct interactions of genes and nc elements
dfr_loops_element_to_gene = dfr_loops_element_to_gene[!dfr_loops_element_to_gene$gene2element %in% annot_flank$gene2element,]
save(dfr_loops_element_to_gene, file = pff("dfr_loops_element_to_gene.rsav"))
save(dfr_loops_to_element, file = pff("dfr_loops_to_element.rsav"))


# create simple lists to get gene/nc interactions
gene_to_nc_loops = split(dfr_loops_element_to_gene[,"nc_id"], dfr_loops_element_to_gene[,"gene_id"])
nc_to_gene_loops = split(dfr_loops_element_to_gene[,"gene_id"], dfr_loops_element_to_gene[,"nc_id"])

save(gene_to_nc_loops, file = pff("gene_to_nc_loops.rsav"))
save(nc_to_gene_loops, file = pff("nc_to_gene_loops.rsav"))
