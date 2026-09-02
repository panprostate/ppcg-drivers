#' 300 series scripts: Motif analysis with FIMO
library(R.utils)


#' Get nucleotide sequence flanking a mutation 
get_mut_seq <- function(mut_gr, flank_bp) {
	# Region of interest
	# mut_gr = GenomicRanges::reduce(mut_gr)
	mut_gr = mut_gr[1]
	expand_gr = mut_gr + flank_bp
	flank_gr = setdiff(expand_gr, mut_gr)
	# Get ref and alt sequence
	ref_seq = paste0(get_hg19_seq(flank_gr[1]), mut_gr$ref, get_hg19_seq(flank_gr[2]))
	alt_seq = paste0(get_hg19_seq(flank_gr[1]), mut_gr$alt, get_hg19_seq(flank_gr[2]))
	return(list(ref_seq, alt_seq))
}


#' HOCOMOCO file locations
fimo.motif_file <- function(core_motif_only = TRUE, annotation = FALSE) {
	# Motif file location
	motif_dir = paste0(gsub('/data$', '', lab_dir), '/pkg/meme_motif')
	if (core_motif_only) {
		motif_file = paste0(motif_dir, 
							'/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme')
	} else {
		motif_file = paste0(motif_dir, 
							'/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme')
	}	
	# 
	if (annotation) {
		motif_file = paste0(motif_dir, 
							'/motif_databases/HUMAN/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv')
	}
	return(motif_file)
}


#' Run FIMO on a mutation
run_fimo <- function(mut_gr, flank_bp = NULL, verbosity = 1, save_to = '', 
						motif_file = fimo.motif_file(core_motif_only = TRUE), 
						make_sequence = TRUE, p_threshold = NULL) {
	if (make_sequence) {
		if (is.null(flank_bp)) {
			stop('run_fimo(): flank_bp arg required for make_sequence = TRUE option')
		}
		# Make flanking seqeunce for FIMO
		seq_list = get_mut_seq(mut_gr, flank_bp)
		ref_seq = seq_list[[1]]; alt_seq = seq_list[[2]]
	} else {
		# Use provided sequence
		ref_seq = mut_gr$ref
		alt_seq = mut_gr$alt
	}
	# Generate unique filename
	file_ok = FALSE
	while (!file_ok) {
		time_str = as.character(as.numeric(Sys.time()))
		time_str = gsub('\\.', '_', time_str)
		seq_file = paste0(lab_dir, '/temp/fimo_seq_', time_str, '.txt')
		if (file.exists(seq_file)) {
			Sys.sleep(0.001)
		} else {
			system(paste0('touch ', seq_file))
			file_ok = TRUE
		}
	}
	# Write sequence to fasta file
	save_txt(c('>ref', ref_seq, '>alt', alt_seq), seq_file)
	# Run FIMO
	output_dir = paste0(lab_dir, '/temp/fimo_out_', time_str)
	res = run_fimo_fasta(fasta_file = seq_file, motif_file = motif_file,
                            output_dir = output_dir, p_threshold = p_threshold, 
                            motif_name = NULL, max_stored_scores = NULL, 
                            verbosity = verbosity, no_run = FALSE)
	# Optionally copy tsv file to another location
	if (save_to != '') {
		unlink(save_to)
		file.copy(paste0(output_dir, '/fimo.tsv'), save_to)
		unlink(output_dir, recursive = TRUE)
	}
	# Clean up temporary files
	unlink(seq_file)
	return(res)
}


#' Compare ref and alt FIMO results to get motif gain/lost
fimo.diff <- function(save_path, overlap_mut_only = FALSE, 
						mut_range = NULL, get_sequence = FALSE) {
	# Load FIMO results
	fimo_df = read.table(save_path, sep = '\t', header = TRUE)
	fimo_df$motif_id = gsub('_HUMAN.H11MO', '', fimo_df$motif_id)
	if (nrow(fimo_df) == 0) return(data.frame())
	# Add placeholder lines in case there are no ref or alt lines
	fimo_df = rbind(
		fimo_df, 
		data.frame(motif_id = 'placeholder', motif_alt_id = 'placeholder', 
					sequence_name = c('ref', 'alt'), start = 0, stop = 0, 
					strand = c('+', '-'), score = 0, p.value = 0, q.value = 0, 
					matched_sequence = 'placeholder')
	)
	# Make unique IDs for dcast
	fimo_df$id = paste(fimo_df$motif_id, fimo_df$start, fimo_df$stop, 
						fimo_df$strand, sep = '_')
	# Make dataframe comparing ref/alt motifs
	dcast_str = 'id ~ sequence_name'
	# dcast_str = 'id + matched_sequence ~ sequence_name'
	diff_df = reshape2::dcast(fimo_df, dcast_str, value.var = 'q.value')
								# fun.aggregate = mean)
	names(diff_df)[2:3] = paste0(names(diff_df)[2:3], '_q')
	#
	for (col_name in c('score', 'p.value')) {
		ref_alt_df = reshape2::dcast(fimo_df, dcast_str, value.var = col_name)
		names(ref_alt_df)[2:3] = paste0(names(ref_alt_df)[2:3], '_', col_name)
		diff_df = dplyr::left_join(diff_df, ref_alt_df, by = 'id')
	}
	if (get_sequence) {
		# seq_df = data.frame(id = '', ref = '', alt = '')
		seq_df = reshape2::dcast(fimo_df, dcast_str, 
									value.var = 'matched_sequence')
		names(seq_df)[2:3] = paste0(names(seq_df)[2:3], '_seq')
		diff_df = dplyr::left_join(diff_df, seq_df, by = 'id')
	}
	# Remove placeholder
	diff_df = diff_df[!grepl('placeholder', diff_df$id), ]
	# Formatting for plot
	diff_df$motif_start = sapply(
		diff_df$id, 
		function(x) as.numeric(unlist(strsplit(x, '_'))[2]), 
		USE.NAMES = FALSE
	)
	diff_df$motif_stop = sapply(
		diff_df$id, 
		function(x) as.numeric(unlist(strsplit(x, '_'))[3]), 
		USE.NAMES = FALSE
	)
	diff_df$motif_id = sapply(
		diff_df$id, 
		function(x) unlist(strsplit(x, '_'))[1], 
		USE.NAMES = FALSE
	)
	diff_df$strand = sapply(
		diff_df$id, 
		function(x) unlist(strsplit(x, '_'))[4], 
		USE.NAMES = FALSE
	)

	diff_df$change = 'Conserved'
	diff_df$change[is.na(diff_df$ref_q)] = 'Gain'
	diff_df$change[is.na(diff_df$alt_q)] = 'Lost'
	# Optionally screen out motifs not overlapping mutations
	if (overlap_mut_only) {
		if (length(mut_range)!= 2) {
			stop(paste0('fimo.diff(): Missing mut_range for option overlap_mut_only = TRUE'))
		}
		diff_df$in_range = !(diff_df$motif_stop < mut_range[1] | diff_df$motif_start > mut_range[2]) 
		diff_df = diff_df[diff_df$in_range, ]
	}
	return(diff_df)
}


#' Overall motif change for multiple occurences overlapping same mutation
fimo_diff.motif_change <- function(motif_key, diff_df) {
	motif_df = diff_df[grepl(motif_key, diff_df$motif_id), ]
	# Categorize change
	if (nrow(motif_df) == 0) {
		output_str = 'None'
	} else {
		change_str = paste(sort(unique(motif_df$change)), collapse = '_')
		if (change_str == 'Conserved') {
			output_str = 'Conserved'
		} else if (change_str == 'Conserved_Gain') {
			output_str = 'Partial_Gain'
		} else if (change_str == 'Conserved_Lost') {
			output_str = 'Partial_Loss'
		} else if (change_str == 'Gain') {
			output_str = 'Gain'
		} else if (change_str == 'Lost') {
			output_str = 'Lost'
		} else if (grepl('Gain', change_str) && grepl('Lost', change_str)) {
			output_str = 'Mixed'
		} else {
			stop(paste0('fimo_diff.motif_change(): unrecognized change_str - ', change_str))
		}
	}
	return(output_str)
}


#' Change Conserved and Altered to No_effect
format_change <- function(df, change_col) {
	df[[change_col]] = gsub('Conserved|Altered|Partial_Loss|Partial_Gain|Mixed', 
							'No_effect', df[[change_col]])
	return(df)
}

#' Plot colors for motif change
motif_change.plot_colors = c(
	# 'Gain' = 'palegreen', 
	'Gain' = 'green', 
	'Lost' = 'red', 				
	'Altered' = 'yellow', 
	'Conserved' = 'gray', 
	'No_effect' = 'gray'
)

#' Get reverse complement for a DNA sequenc
reverse_complement <- function(input_str) {
	output_str = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(input_str)))
	return(output_str)
}

#' Load HOCOMOCO PWMs
load.hocomoco_motifs <- function() {
	fimo_obj = universalmotif::read_meme(fimo.motif_file(core_motif_only = TRUE))
	# Extract positional weight matrix to list
	pwm_list = list()
	for (fimo_item in fimo_obj) {
		mname = gsub('_HUMAN.H11MO', '', fimo_item@name)
		pwm_list[[mname]] = fimo_item@motif
	}
	return(pwm_list)
}


#' Write sequences to a FASTA file
#' This function writes a set of sequences and their corresponding identifiers to a FASTA file.
#' @param seq_vec A string vector containing the sequences to be written to the FASTA file.
#' @param id_vec A string vector containing the identifiers for each sequence. The length of `id_vec` should be the same as `seq_vec`.
#' @param file_path Output FASTA file path.
#' @return None. This function writes the sequences to the specified file.
#' @examples
#' seqs <- c("ATGCGT", "CGTACG", "TTAGGC")
#' ids <- c("seq1", "seq2", "seq3")
#' write_fasta(seqs, ids, "output.fasta")
write_fasta <- function(seq_vec, id_vec, file_path) {
	# Check number of sequences and identifiers equal
    if (length(seq_vec) != length(id_vec)) {
        stop("The length of seq_vec must be equal to the length of id_vec")
    }
	# Put together FASTA lines
    fasta_lines <- mapply(function(seq, id) {
        paste0(">", id, '\n', seq) # Second '\n' is added by writeLines in save_txt
    }, seq_vec, id_vec)
	# Write to file
    save_txt(fasta_lines, file_path)
}


#' Run FIMO on a fasta file - base wrapper of shell FIMO
#' @param fasta_file Input FASTA file path
#' @param motif_file Input motif file path. If NULL, use fimo.motif_file() function default.
#' @param output_dir A character string specifying the path to the output directory.
#' @param p_threshold Optionally override FIMO default p-value threshold. 
#' @param motif_name Optionally specify a single motif name to run on.
#' @param no_run Optionally do not run FIMO and only output the command string.
#' @return The function does not return a value but generates output files as specified by FIMO.
run_fimo_fasta <- function(fasta_file, motif_file = NULL,
                            output_dir = '', 
                            p_threshold = NULL, 
                            motif_name = NULL, 
                            max_stored_scores = NULL, 
                            verbosity = 1, no_run = FALSE) {
    unlink(output_dir, recursive = TRUE)
    # --oc allows overwriting if results folder already exists
    cmd_str = paste('fimo', '-oc', output_dir, '--verbosity', verbosity)
    if (!is.null(p_threshold)) { 
        cmd_str = paste(cmd_str, '--thresh', 
                        format(p_threshold, scientific = FALSE))
    }
    if (!is.null(motif_name)) {
        cmd_str = paste(cmd_str, '--motif', motif_name) 
    }
    if (!is.null(max_stored_scores)) {
        cmd_str = paste(cmd_str, '--max-stored-scores', 
                        format(max_stored_scores, scientific = FALSE))
    }
    if (is.null(motif_file)) {
        motif_file = fimo.motif_file(core_motif_only = TRUE) # Default motif catalog
    }
    cmd_str = paste(cmd_str, motif_file, fasta_file)
    # Optionally skip running command and output only command string
    res = ifelse(no_run, -1, system(cmd_str))
    return(list(cmd_str = cmd_str, sys.return = res))
}
