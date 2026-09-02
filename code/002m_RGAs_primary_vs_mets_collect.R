source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
load(file = pff("patient_sets_for_drivers.rsav"))



# protein-coding drivers separately
CDS_drivers =  gsub("__SNV_CDS", "", grep("SNV_CDS", names(patient_sets_for_drivers), value = T))
load(pff("dnds_results_mets_annotations_all.rsav"))
dnds_results_mets_annotations_all = dnds_results_mets_annotations_all[dnds_results_mets_annotations_all$gene %in% CDS_drivers,] 
SNV_CDS_drivers_mets = lapply(split(dnds_results_mets_annotations_all$sampleID, dnds_results_mets_annotations_all$gene), unique)
SNV_CDS_drivers_mets = lapply(SNV_CDS_drivers_mets, function(x)  gsub("(.+)._DNA", "\\1", x) )
names(SNV_CDS_drivers_mets) = paste0(names(SNV_CDS_drivers_mets), "__SNV_CDS")

# CNA drivers
load(pff("genes_with_full_CNAs_mets.rsav"))
todo = rbind(
	c("__loss", "all_loss"),
	c("__gain", "all_gain"),
	c("__hAMP", "high_gain"))

get_mets_cna_events = function(i, todo, patient_sets_for_drivers, genes_with_full_CNAs_mets) {
	driver_list_tag = todo[i, 1]
	cna_table_tag = todo[i, 2]
	
	event_drivers =  gsub(driver_list_tag, "", grep(driver_list_tag, names(patient_sets_for_drivers), value = T))
	
	event_drivers_mets = genes_with_full_CNAs_mets[genes_with_full_CNAs_mets$gene %in% event_drivers & 
			genes_with_full_CNAs_mets$tag == cna_table_tag,]
	
	event_drivers_mets = sapply(split(event_drivers_mets$patient, event_drivers_mets$gene), unique)

	names(event_drivers_mets) = paste0(names(event_drivers_mets), driver_list_tag)
	event_drivers_mets
}


CNA_drivers_mets = lapply(1:nrow(todo), get_mets_cna_events, todo, patient_sets_for_drivers, genes_with_full_CNAs_mets)
names(CNA_drivers_mets) = todo[,1]

gain_drivers_mets = CNA_drivers_mets[["__gain"]]
loss_drivers_mets = CNA_drivers_mets[["__loss"]]
hAMP_drivers_mets = CNA_drivers_mets[["__hAMP"]]


# nc-SNV drivers
ncSNV_drivers =  gsub("__SNV_NC", "", grep("SNV_NC", names(patient_sets_for_drivers), value = T))
load(file = pff("results_signf_merged_annot.rsav"))
load(file = pff("variants_to_elements_mets_snv.rsav"))

get_patients_with_nc_muts = function(mut_tag, results_signf_merged_annot, nc_drivers, variants_to_elements_mets) {

	nc_driver_elements = results_signf_merged_annot[results_signf_merged_annot$mut_type == mut_tag,]
	nc_driver_elements = split(nc_driver_elements$element_ids, nc_driver_elements$annots_MAIN)
	nc_driver_elements = sapply(nc_driver_elements, function(x) unique(unlist(strsplit(x, split = ","))))
	nc_driver_elements = nc_driver_elements[nc_drivers]

	nc_drivers_mets = lapply(names(nc_driver_elements), function(driver_id)  
			do.call(rbind, variants_to_elements_mets_snv[nc_driver_elements[[driver_id]]] ))
	names(nc_drivers_mets) = names(nc_driver_elements)
	nc_drivers_mets = sapply(sapply(nc_drivers_mets, '[[', 'patient'), unique)
	nc_drivers_mets	
}


SNV_NC_drivers_mets = get_patients_with_nc_muts("SNV_NC", results_signf_merged_annot, ncSNV_drivers, variants_to_elements_mets_snv)
names(SNV_NC_drivers_mets) = paste0(names(SNV_NC_drivers_mets), "__SNV_NC")

# SV drivers
load(file = pff("variants_to_elements_mets_sv.rsav"))
SV_drivers =  gsub("__SV", "", grep("__SV", names(patient_sets_for_drivers), value = T))
SV_drivers_mets = get_patients_with_nc_muts("SV", results_signf_merged_annot, SV_drivers, variants_to_elements_mets_sv)
names(SV_drivers_mets) = paste0(names(SV_drivers_mets), "__SV")
# add ETS separately
load(pff("patients_mets_ETS.rsav"))
SV_drivers_mets[['ETS__SV']] = patients_mets_ETS

# BI drivers
load(pff("BI_patient_sets_mets.rsav"))
BI_patient_sets_mets = lapply(BI_patient_sets_mets, function(x) gsub("(.+)._DNA", "\\1", x)  )
names(BI_patient_sets_mets) = paste0(names(BI_patient_sets_mets), "__BI")

# combine all together
patient_sets_for_drivers_mets = c(
			SV_drivers_mets, SNV_NC_drivers_mets, SNV_CDS_drivers_mets, 
			gain_drivers_mets, loss_drivers_mets, hAMP_drivers_mets,
			BI_patient_sets_mets)
save(patient_sets_for_drivers_mets, file = pff("patient_sets_for_drivers_mets.rsav"))
