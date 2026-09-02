#' FIMO analysis on selected recurrent mutations
source('FigureS3_utils.R')
source('FigureS3_fimo_utils.R')
library(R.utils)
library(ggseqlogo)
library(cowplot)


# --- Functions ---
#' Set all non-matching PWM values to 0 according to a sequence
modify_pwm <- function(pwm, seq, identity_mat = FALSE) {
	chars <- unlist(strsplit(seq, ''))
	new_pwm <- pwm
	for (i in 1:ncol(new_pwm)) {
		for (j in 1:nrow(new_pwm)) {
			if (rownames(new_pwm)[j] != chars[i]) {
				new_pwm[j, i] <- 0
			}
			# Set matching nucleotide to 1 - for plotting sequence
			if (identity_mat & rownames(new_pwm)[j] == chars[i]) {
				new_pwm[j, i] <- 1
			}
		}
	}
	return(new_pwm)
}


# --- Plot functions/objects ---
#' Theme for logo plots
custom_logo_theme <- function(p, nbase) {
	return(p + scale_x_discrete(limits = c(1:nbase)) + 
			theme_bw()) 
}

#' Add global labels to grid plots
plot_grid.add_labels <- function(plot_p, xlab, ylab, text_size = 12) {
	x_grob = grid::textGrob(xlab, gp = grid::gpar(col = 'black', fontsize = text_size))
	y_grob = grid::textGrob(ylab, gp = grid::gpar(col = 'black', fontsize = text_size), 
							rot = 90)
	plot_p = gridExtra::grid.arrange(gridExtra::arrangeGrob(plot_p, left = y_grob, bottom = x_grob))
	return(plot_p)
}

#' Rotate y-axis title
rotate_y_title = theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))

#' Remove y-axis grid lines
rm_grids = theme(panel.grid.major.y = element_blank(), 
					panel.grid.minor.y = element_blank())

#' Remove all grid lines
rm_all_grids = theme(panel.grid.major.y = element_blank(), 
						panel.grid.minor.y = element_blank(), 
						panel.grid.major.x = element_blank(), 
						panel.grid.minor.x = element_blank())

#' Remove x-axis ticks
rm_ticks = theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
					axis.text.y = element_blank(), axis.ticks.y = element_blank())


# PARAMS
flank_bp = 25
overlap_mut_only = TRUE
plots_path = new.dir(paste0(plot_dir, '/ppcg/drivers_MS_repo'))
save_folder = new.dir(paste0(plots_path, '/fimo_run'))
no_filter_save_folder = new.dir(paste0(plots_path, '/fimo_run_no_filter'))


# --- Selected recurrent mutations ---
#' chr15 57178146 C > T, 16 samples on ZNF280D
#' chr1 202447842 202447842 T > C, 8 samples on PPP1R12B
mut_df = data.frame(
	chr   = c('chr15', 'chr1'),
	start = c(57178146, 202447842),
	end   = c(57178146, 202447842),
	ref   = c('C', 'T'),
	alt   = c('T', 'C'),
	stringsAsFactors = FALSE
)
# Make GRanges
mut_df$mut_id = mut_index(mut_df, id_col = NULL)
mut_gr = mut_df2gr(mut_df, add_col = TRUE)

# HOCOMOCOv11 core motifs	
motif_file = fimo.motif_file(core_motif_only = TRUE)

# Load PWMs
pwm_list = load.hocomoco_motifs()


# --- Run FIMO on each mutation ---
id_list = unique(mut_df$mut_id)
i = 0; skipped_n = 0
cat(paste0(length(id_list), ' mutations to process\n'))
for (mut_id in id_list) {
	save_path = paste0(save_folder, '/', mut_id, '_fimo.tsv')
	this_mut_gr = mut_gr[mut_gr$mut_id == mut_id]
	# Run FIMO
	fimo_output = run_fimo(mut_gr = this_mut_gr, flank_bp = flank_bp, 
							motif_file = motif_file, save_to = save_path, 
							verbosity = 1)
	# Run a no-p-filter version to get ref sequence motif scores 
	fimo_output = run_fimo(mut_gr = this_mut_gr, flank_bp = flank_bp, 
							motif_file = motif_file, 
							save_to = paste0(no_filter_save_folder, '/', mut_id, '_fimo.tsv'), 
							p_threshold = 1, verbosity = 1)
	i = i + 1; cat(paste0('\rProcessed ', i, ' mutations, skipped ', skipped_n))
}
cat('\n')


# --- Combine FIMO results ---
results_df = data.frame()
missing_fimo = c()
pb = progress_bar_init(length(id_list), msg = 'Motif analysis')
for (mut_id in id_list) {
	fimo_file = paste0(save_folder, '/', mut_id, '_fimo.tsv')
	# Skip if motif file missing
	if (!file.exists(fimo_file)) {
		missing_fimo = c(missing_fimo, fimo_file)
		pb$tick()
		next
	}
	# No motif case 
	line_count = system(paste0("grep -v '#' ", fimo_file, "| wc -l"), intern = TRUE)
	if (line_count <= 1) {
		pb$tick()
		next
	}
	# Mutation index in FIMO output
	this_mut_df = mut_df[mut_df$mut_id == mut_id, ]
	ref_len = nchar(this_mut_df$ref)
	alt_len = nchar(this_mut_df$alt)
	mut_range = c(flank_bp + 1, flank_bp + max(ref_len, alt_len))
	# Ref-alt motif difference
	diff_df = fimo.diff(fimo_file, overlap_mut_only = overlap_mut_only, 
						mut_range = mut_range)
	diff_df$mut_id = mut_id
	results_df = rbind(results_df, diff_df)
	pb$tick()
}
pb$terminate()
cat(paste0(length(missing_fimo), ' motif files missing\n'))
# Save
write.csv(results_df, file = paste0(plots_path, '/fimo_results.csv'), 
			row.names = FALSE)


# --- Logo plots ---
# Find gene names corresponding to motifs
hcmc_annot_df = read.table(fimo.motif_file(annotation = TRUE), sep = '\t', header = TRUE)
hcmc_annot_df$formatted_motif_id = gsub('_HUMAN.H11MO', '', hcmc_annot_df$Model)
results_df$gene_name = hcmc_annot_df$Transcription.factor[match(results_df$motif_id, hcmc_annot_df$formatted_motif_id)]
for (i in 1:nrow(results_df)) {
	plot_df1 = results_df[i, ]
	gene_name = plot_df1$gene_name
	motif_name = plot_df1$motif_id
	mut_id = plot_df1$mut_id
	pwm = pwm_list[[motif_name]]
	motif_length = ncol(pwm)
	# Calculate site position
	plot_df1$mut_start = 26; plot_df1$mut_end = 26
	plot_df1$mut_start_in_motif = plot_df1$mut_start - plot_df1$motif_start + 1
	plot_df1$mut_end_in_motif = plot_df1$mut_end - plot_df1$motif_start + 1
	plot_df1$site_pos = with(plot_df1, ifelse(strand == '-', motif_length - mut_start_in_motif + 1, mut_start_in_motif))
	plot_df1$mut_pos = plot_df1$site_pos
	# Get ref and alt sequence
	fimo_df = read.table(paste0(no_filter_save_folder, '/', mut_id, '_fimo.tsv'), 
							sep = '\t', header = TRUE)
	fimo_df$motif_id = gsub('_HUMAN.H11MO', '', fimo_df$motif_id)
	fimo_df$id = paste(fimo_df$motif_id, fimo_df$start, fimo_df$stop, 
						fimo_df$strand, sep = '_')
	ref_seq = fimo_df$matched_sequence[(fimo_df$id == plot_df1$id & fimo_df$sequence_name == 'ref')]
	alt_seq = fimo_df$matched_sequence[(fimo_df$id == plot_df1$id & fimo_df$sequence_name == 'alt')]
	# Get a 'sequence only pwm'
	ref_pwm <- modify_pwm(pwm, ref_seq, identity_mat = TRUE)
	alt_pwm <- modify_pwm(pwm, alt_seq, identity_mat = TRUE)

	# --- Plot PWM ---
	y_max = 2
	plot_title = paste0(mut_id, ' (', gene_name, ' motif ', plot_df1$change, ')')
	# Motif PWM
	p1 = ggseqlogo(pwm) + 
		labs(x = NULL, y = paste0(gene_name, '\nmotif')) + 
		ggtitle(plot_title) + 
		ylim(0, y_max) + 
		theme(plot.title = element_text(hjust = 0.5, size = 10)) # centre title
	p1 = custom_logo_theme(p1, motif_length) + 
		rotate_y_title + 
		theme()
	# Ref PWM
	p2 = ggseqlogo(ref_pwm, method = 'custom') + 
		labs(x = NULL, y = 'Reference', title = NULL) + 
		theme()
	p2 = custom_logo_theme(p2, motif_length) + 
		rm_ticks + 
		rm_all_grids + 
		rotate_y_title + 
		theme()
	# Alt PWM
	p3 = ggseqlogo(alt_pwm, method = 'custom') + 
		labs(x = NULL, y = 'Mutated', title = NULL) + 
		theme()
	p3 = custom_logo_theme(p3, motif_length) + 
		rm_ticks + 
		rm_all_grids + 
		rotate_y_title + 
		theme()
	# Highlight site position
	highlight_site_pos <- geom_rect(aes(xmin = unique(plot_df1$site_pos) - 0.5, 
										xmax = unique(plot_df1$site_pos) + 0.5, 
										ymin = -Inf, ymax = Inf),
									fill = 'salmon', alpha = 0.3, inherit.aes = FALSE)
	p1 <- p1 + highlight_site_pos
	p2 <- p2 + highlight_site_pos
	p3 <- p3 + highlight_site_pos
	# Stitch and label
	plot_p = plot_grid(p1, p2 + rm_legend, p3 + rm_legend, ncol = 1, align = 'v', 
						rel_heights = c(4, 1, 1))
	plot_p = plot_grid.add_labels(plot_p, xlab = 'Motif position', ylab = NULL, 
									text_size = 12)
	pdf_file = paste0(plots_path, '/', mut_id, '_', plot_df1$id, '.pdf')
	pdf(pdf_file, width = (ncol(pwm) + 3) * 0.3, height = 3.5)
	grid::grid.draw(plot_p)
	dev.off()
}


