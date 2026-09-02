source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")

##
## create mapping of elements and SNVs/indels
##
variants_to_elements = list()

print(load(file = pff("gr_prepared_elements.rsav")))
gr_elements_for_SNV = gr_prepared_elements[!grepl("^CDSgene::", gr_prepared_elements$id)]
gr_elements_for_SV = gr_prepared_elements[!grepl("^CDS::", gr_prepared_elements$id)]

print(load(file = pff("prepared_variants__SNV_indel.rsav")))
print(load(file = pff("gr_prepared_variants__SNV_indel.rsav")))

ov = findOverlaps(gr_elements_for_SNV, gr_prepared_variants)
variants_to_elements_snv = cbind(
		element_id = gr_elements_for_SNV$id[queryHits(ov)], 
		prepared_variants[subjectHits(ov),],
		stringsAsFactors = FALSE
)

##
# load splicing alterations from DnDS-cv; add these again to CDS variant dataset for consistency
## 
load(file = pff("dnds_results_annotations_all.rsav"))
splice_annotations = dnds_results_annotations_all[dnds_results_annotations_all$impact == "Essential_Splice",]

variants_splice_tag = paste(
		splice_annotations$sampleID,
		paste0("chr", splice_annotations$chr),
		splice_annotations$pos,
		splice_annotations$ref,
		splice_annotations$mut, 
		sep = "__")
all_variants_tag = paste(
		prepared_variants$patient, 
		prepared_variants$chr, 
		prepared_variants$pos1, 
		prepared_variants$ref, 
		prepared_variants$alt, 
		sep = "__")
all_element_variants_tag = paste(
		variants_to_elements_snv$patient, 
		variants_to_elements_snv$chr, 
		variants_to_elements_snv$pos1, 
		variants_to_elements_snv$ref, 
		variants_to_elements_snv$alt, 
		sep = "__")
		
splice_variants = prepared_variants [all_variants_tag %in% variants_splice_tag,]
splice_variants$Func.refGene = splice_variants$ExonicFunc.refGene = "EssentialSplicing_dndscv"

# map symbols to full CDS IDs
CDS_full_ids = unique(grep("^CDS::", gr_prepared_elements$id, value = T))
symbol_map = structure(CDS_full_ids, names = gsub("CDS::gencode::(.+)::(.+)", "\\1", CDS_full_ids))
splice_variants = cbind(element_id = symbol_map[splice_variants$Gene.refGene], splice_variants)
variants_to_elements_snv_minus_previous_splice = variants_to_elements_snv[!all_element_variants_tag %in% variants_splice_tag,]
# create new variants_to_elements_snv object by adding splice-annotated variants from DNDScv
variants_to_elements_snv = rbind(variants_to_elements_snv_minus_previous_splice, splice_variants)

# exclude silent mutations from protein-coding sequnce
silent_annotations = c("", "unknown", "synonymous SNV")
which_silent = which(grepl("CDS::", variants_to_elements_snv$element_id) & variants_to_elements_snv$ExonicFunc.refGene %in% silent_annotations)
variants_to_elements_snv = variants_to_elements_snv[-which_silent,]

variants_to_elements[['SNV_indel']] = 
		split(variants_to_elements_snv, variants_to_elements_snv$element_id)



##
## create mapping of elements and SV BPs
##
load(file = pff("prepared_variants_SV.rsav"))
load(file = pff("gr_prepared_variants__SV.rsav"))

ov = findOverlaps(gr_elements_for_SV, gr_prepared_variants)
variants_to_elements_sv = cbind(
		element_id = gr_elements_for_SV$id[queryHits(ov)], 
		prepared_variants_SV[subjectHits(ov),],
		stringsAsFactors = FALSE
)
variants_to_elements[['SV']] = split(variants_to_elements_sv, variants_to_elements_sv$element_id)
save(variants_to_elements, file = pff("variants_to_elements.rsav"))




#
# create first to-do list that contains all elements to be tested
#

# for SVs
SV_elements = unique(names(variants_to_elements[['SV']]))

# for SNVs
# remove CDS mutations from here - are analysed using dndscv
SNV_elements = unique(names(variants_to_elements[['SNV_indel']]))
SNV_elements = grep("^CDS::", SNV_elements, invert = TRUE, value = TRUE)


todo_list = rbind(
		cbind(muts_dataset = "SNV_indel", element = SNV_elements),
		cbind(muts_dataset = "SV", element = SV_elements))
todo_list = data.frame(id = 1:nrow(todo_list), todo_list, win_size = 50000, stringsAsFactors = FALSE)
colnames(todo_list) = c("id", "muts_dataset", "element", "win_size")

# SVs need a larger window
todo_list[todo_list$muts_dataset == "SV", "win_size"] = 500000

system(paste("mkdir -p ", pff("ADWGS_tmp/")))
save(todo_list, file = pff("ADWGS_tmp/todo_list.0.rsav"))


print(load(file = pff("gr_prepared_elements.rsav")))
gr_prepared_elements = gr_prepared_elements[!grepl("^CDS::", gr_prepared_elements$id)]
gr_elements_for_SNV = gr_prepared_elements[!grepl("^CDSgene::", gr_prepared_elements$id)]
gr_elements_for_SV = gr_prepared_elements

## make silent todo list with no mutations
SNV_elements_silent = setdiff(gr_elements_for_SNV$id, SNV_elements)
SV_elements_silent = setdiff(gr_elements_for_SV$id, SV_elements)

todo_list_silent = rbind(
		cbind(muts_dataset = "SNV_indel", element = SNV_elements_silent),
		cbind(muts_dataset = "SV", element = SV_elements_silent))

todo_list_silent = data.frame(id = 1:nrow(todo_list_silent), todo_list_silent, win_size = 50000, 
		stringsAsFactors = FALSE)
colnames(todo_list_silent) = c("id", "muts_dataset", "element", "win_size")

todo_list_silent[todo_list_silent$muts_dataset == "SV", "win_size"] = 500000

save(todo_list_silent, file = pff("ADWGS_tmp/todo_list_silent.rsav"))
