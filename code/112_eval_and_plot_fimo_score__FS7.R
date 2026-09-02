#' Plot FIMO scores and q-values
source('FigureS3_utils.R')
source('FigureS3_fimo_utils.R')
library(R.utils)


#
plots_path = new.dir(paste0(plot_dir, '/ppcg/drivers_MS_repo'))
no_filter_save_folder = new.dir(paste0(plots_path, '/fimo_run_no_filter'))
flank_bp = 25
overlap_mut_only = TRUE


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


# --- Get FIMO scores for both ref and alt ---
id_list = unique(mut_df$mut_id)
results_df = data.frame()
missing_fimo = c()
pb = progress_bar_init(length(id_list), msg = 'Motif analysis')
for (mut_id in id_list) {
	fimo_file = paste0(no_filter_save_folder, '/', mut_id, '_fimo.tsv')
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


# --- Score bar plots ---
# Find gene names corresponding to motifs
hcmc_annot_df = read.table(fimo.motif_file(annotation = TRUE), sep = '\t', header = TRUE)
hcmc_annot_df$formatted_motif_id = gsub('_HUMAN.H11MO', '', hcmc_annot_df$Model)
results_df$gene_name = hcmc_annot_df$Transcription.factor[match(results_df$motif_id, hcmc_annot_df$formatted_motif_id)]
# Load results from default param FIMO run
def_results_df = read.csv(paste0(plots_path, '/fimo_results.csv'))
def_results_df$ref_q[is.na(def_results_df$ref_q)] <- 1
# Plot FIMO score
pdf(paste0(plots_path, '/score_bar.pdf'), width = 2, height = 2.5)
for (id in c('NR1I2.0.C_15_33_-', 'NR1I3.0.C_15_32_-')) {
	sub_df <- results_df[results_df$id == id, ]
	score_df = reshape2::melt(sub_df[c('id', 'ref_score', 'alt_score')], id = 'id')
	score_df$variable <- factor(score_df$variable, 
								levels = c('ref_score', 'alt_score'),
								labels = c('Reference', 'Mutated'))
	# Add default params run q-values
	score_df$q_value = ifelse(
		score_df$variable == 'Reference',
		def_results_df$ref_q[match(sub_df$id, def_results_df$id)],
		def_results_df$alt_q[match(sub_df$id, def_results_df$id)]
	)
	score_df$asterisk = pval2asterisk(score_df$q_value, ns = 'n.s.') # Asterisks for plot
	p <- ggplot(score_df, aes(x = variable, y = value, fill = variable)) +
		geom_bar(stat = 'identity') +
		geom_text(aes(label = asterisk, y = value + 0.1), vjust = 0, size = 3) +
		labs(title = paste0(sub_df$gene_name, ' motif'), x = NULL, y = 'FIMO score') +
		scale_fill_manual(values = c('Reference' = 'grey', 'Mutated' = 'orange')) +
		ylim(0, 12.5) + 
		theme_bw()
	print(p + rm_legend)
}
dev.off()
