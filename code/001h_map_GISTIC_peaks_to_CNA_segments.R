source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")


#
# annotate significant cytobands from GISTIC2
#

FDR_cutoff = 0.05
add_flank_fraction = 0.2	# percentage of genomic space around GISTIC peak calls, to map key genes in poorly defined peaks

significant_CNAs_GISTIC = read.delim(pff("GISTIC2_run1/GISTIC2_results/PPCG_primaryOnly.all_lesions.conf_90.txt"), stringsAsFactors = FALSE)
significant_CNAs_GISTIC = significant_CNAs_GISTIC[grep("CN values", significant_CNAs_GISTIC$Unique.Name, invert = TRUE),]
significant_CNAs_GISTIC = significant_CNAs_GISTIC[
		significant_CNAs_GISTIC$Residual.q.values.after.removing.segments.shared.with.higher.peaks < FDR_cutoff, ]

# x = significant_CNAs_GISTIC$Peak.Limits[1]
get_coords = function(x, add_flank_fraction) {
	chr = gsub("(.+):(.+)-(.+)\\((.+)", "\\1", x)
	start = as.numeric(gsub("(.+):(.+)-(.+)\\((.+)", "\\2", x))
	end = as.numeric(gsub("(.+):(.+)-(.+)\\((.+)", "\\3", x))
	
	add_flank_bps = round((end - start) * (add_flank_fraction))
	start = start - add_flank_bps
	end = end + add_flank_bps
	
	start[start < 1] = 1
	
	data.frame(chr, start, end, stringsAsFactors = FALSE)
}



peak_coords = do.call(rbind, lapply(significant_CNAs_GISTIC$Wide.Peak.Limits, get_coords, add_flank_fraction))
significant_CNAs_GISTIC = cbind(peak_coords, significant_CNAs_GISTIC)

# get samples annotated to the peak
PPCG_colnames = grep("^PPCG", colnames(significant_CNAs_GISTIC), value = TRUE)

# i = 1
get_patients = function(i, significant_CNAs_GISTIC, PPCG_colnames) {
	cat(i, " ")
	tumor_ids = names(which(significant_CNAs_GISTIC[i, PPCG_colnames, drop = TRUE] != 0))
	tumor_ids = unique(tumor_ids)
	patient_ids = unique(gsub("._DNA$", "", tumor_ids))
	tumor_ids_concat = paste(tumor_ids, collapse = ",")
	patient_ids_concat = paste(patient_ids, collapse = ",")
	n_patients = length(patient_ids)
	data.frame(n_patients, tumors = tumor_ids_concat, patients = patient_ids_concat, stringsAsFactors = FALSE)
}

patients_info = do.call(rbind, lapply(1:nrow(significant_CNAs_GISTIC), get_patients, significant_CNAs_GISTIC, PPCG_colnames))
significant_CNAs_GISTIC = cbind(significant_CNAs_GISTIC, patients_info)

significant_CNAs_GISTIC$fdr = significant_CNAs_GISTIC$q.values
significant_CNAs_GISTIC$fdr_residual = significant_CNAs_GISTIC[, "Residual.q.values.after.removing.segments.shared.with.higher.peaks"]
significant_CNAs_GISTIC$CNA_type = substr(significant_CNAs_GISTIC$Unique.Name, 1, 3)

significant_CNAs_GISTIC$id = gsub(" ", "", paste0("chr", significant_CNAs_GISTIC$Descriptor, "_", significant_CNAs_GISTIC$CNA_type))

cols_to_select = c("id", "chr", "CNA_type", "start", "end", "fdr", "fdr_residual", "n_patients", "tumors", "patients")
significant_CNAs_GISTIC = significant_CNAs_GISTIC[, cols_to_select]
save(significant_CNAs_GISTIC, file = pff("significant_CNAs_GISTIC.rsav"))

load(file = pff("prepared_CNAs.rsav"))
gr_prepared_CNAs = GRanges(prepared_CNAs$chr, IRanges(prepared_CNAs$startpos, prepared_CNAs$endpos), 
		patient = prepared_CNAs$patient, annot = prepared_CNAs$annot, cn_width = prepared_CNAs$cn_width)

# id = "chr1q21.3_Amp"
get_gistic_patient_ids = function(id, significant_CNAs_GISTIC, gr_prepared_CNAs) {
	cat(id, " ")
	this_peak = significant_CNAs_GISTIC[significant_CNAs_GISTIC$id == id,]
	CNA_type = this_peak$CNA_type
	
	annot_to_select = NA
	if (CNA_type == "Amp") {
		annot_to_select = c("gain", "high_gain")		
	}
	if (CNA_type == "Del") {
		annot_to_select = c("loss", "full_loss")
	}	
	gr_this_peak = GRanges(this_peak$chr, IRanges(this_peak$start, this_peak$end))

	gr_CNAs = gr_prepared_CNAs[gr_prepared_CNAs$annot %in% annot_to_select, ]
	gr_CNAs = gr_CNAs[subjectHits(findOverlaps(gr_this_peak, gr_CNAs))]
	
	patients_from_segs = unique(gr_CNAs$patient)
	n_patients_from_segs = length(patients_from_segs)
	
	this_peak$n_patients_segs = n_patients_from_segs
	this_peak$patients_segs = paste(patients_from_segs, collapse = ",")
	this_peak
}

significant_CNAs_GISTIC_with_segs = do.call(rbind, 
		lapply(significant_CNAs_GISTIC$id, get_gistic_patient_ids, significant_CNAs_GISTIC, gr_prepared_CNAs))

save(significant_CNAs_GISTIC_with_segs, file = pff("significant_CNAs_GISTIC_with_segs.rsav"))
