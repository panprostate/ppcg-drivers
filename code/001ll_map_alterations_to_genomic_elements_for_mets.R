source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")

# create mapping of elements and SNVs/indels in mets
load(file = pff("gr_prepared_elements.rsav"))
gr_elements_for_SNV = gr_prepared_elements[!grepl("^CDSgene::", gr_prepared_elements$id)]
gr_elements_for_SV = gr_prepared_elements[!grepl("^CDS::", gr_prepared_elements$id)]

load(file = pff("variants_mets__SNV_indel.rsav"))
gr_variants_mets_SNV_indel = GRanges(variants_mets__SNV_indel$chr, 
		IRanges(variants_mets__SNV_indel$pos1, variants_mets__SNV_indel$pos2))
ov = findOverlaps(gr_elements_for_SNV, gr_variants_mets_SNV_indel)
variants_to_elements_mets_snv = cbind(
		element_id = gr_elements_for_SNV$id[queryHits(ov)], 
		variants_mets__SNV_indel[subjectHits(ov),],
		stringsAsFactors = FALSE
)
variants_to_elements_mets_snv = split(variants_to_elements_mets_snv, variants_to_elements_mets_snv$element_id)
save(variants_to_elements_mets_snv, file = pff("variants_to_elements_mets_snv.rsav"))


# create mapping of elements and SV BPs

load(file = pff("variants_mets__SV.rsav"))
gr_variants_mets_SV = GRanges(variants_mets__SV$chr, 
		IRanges(variants_mets__SV$start, variants_mets__SV$end))
ov = findOverlaps(gr_elements_for_SV, gr_variants_mets_SV)

variants_to_elements_sv = cbind(
		element_id = gr_elements_for_SV$id[queryHits(ov)], 
		variants_mets__SV[subjectHits(ov),],
		stringsAsFactors = FALSE
)

variants_to_elements_mets_sv = split(variants_to_elements_sv, variants_to_elements_sv$element_id)
save(variants_to_elements_mets_sv, file = pff("variants_to_elements_mets_sv.rsav"))
