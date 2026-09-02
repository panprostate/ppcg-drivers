#' Check ARBS neighboring SLC family genes
source('FigureS3_utils.R')


plots_path = new.dir(paste0(plot_dir, '/ppcg/drivers_MS_repo'))
dist_threshold = 1 * 10^3
cat(paste0('Distance threshold = ', dist_threshold, ' bp\n'))


# --- Load AR binding sites ---
# Load TARBS/NARBS (tumor-specific or normal tissue-specific binding sites)
narbs_df = read_sites(paste0(sites_dir, '/ppcg/NARBS.bed'), 
                mode = 'nohead_open_chrom')
narbs_gr = mut_df2gr(narbs_df, no_seq = TRUE, add_col = TRUE, use_ref = TRUE)
tarbs_df = read_sites(paste0(sites_dir, '/ppcg/TARBS.bed'), 
                mode = 'nohead_open_chrom')
tarbs_gr = mut_df2gr(tarbs_df, no_seq = TRUE, add_col = TRUE, use_ref = TRUE)
# Load Chip-seq data (contains also tumor and normal common binding sites)
file_list = list.files(paste0(sites_dir, '/ppcg/tarbs_narbs'), pattern = 'DF', full.names = TRUE)
cat('Loading ChIPseq data...')
chip_df = data.frame()
for (file_path in file_list) {
    this_df = read_sites(file_path, mode = 'nohead_open_chrom', all_info = TRUE)
    names(this_df)[4:5] = c('site', 'score')
    this_df$sample = basename(file_path)
    chip_df = rbind(chip_df, this_df)
}
chip_df$sample = gsub('_peaks.bed.txt', '', chip_df$sample)
chip_df$sample_type = ifelse(grepl('tumor', chip_df$sample), 'tumor', 'normal')
chip_df$coded_sample = sapply(chip_df$site, function(x) unlist(strsplit(x, '_'))[1])
chip_df$rep = sapply(chip_df$site, function(x) stringr::str_match(x, 'rep[0-9+]'))
chip_gr = mut_df2gr(chip_df, no_seq = TRUE, use_ref = TRUE, add_col = TRUE)
cat('Done\n')


# --- Combine all available data ---
merged_gr = GenomicRanges::reduce(c(chip_gr, narbs_gr, tarbs_gr))
# Annotate NARBS and TARBS
merged_gr$annotation = 'both'
overlap_narbs = findOverlaps(merged_gr, narbs_gr)
overlap_tarbs = findOverlaps(merged_gr, tarbs_gr)
merged_gr$annotation[queryHits(overlap_narbs)] = 'narbs'
merged_gr$annotation[queryHits(overlap_tarbs)] = 'tarbs'


# --- Annotate binding scores in merged_gr ---
total_tumor_sample = length(unique(chip_df$sample[grepl('tumor', chip_df$sample)]))
total_normal_sample = length(unique(chip_df$sample[grepl('normal', chip_df$sample)]))

combined_metrics_df = data.frame(peak_id = 1:length(merged_gr), 
                                    n_overlap = 0, n_tumor = 0, n_normal = 0, 
                                    mean_tumor_score = 0, mean_normal_score = 0, 
                                    mean_overall_score = 0)
pb = progress_bar_init(length(merged_gr), msg = 'Annotating merged sites')
for (i in 1:length(merged_gr)) {
    overlap_gr = subsetByOverlaps(chip_gr, merged_gr[i])
    tumor_num = sum(grepl('tumor', unique(overlap_gr$sample)))
    normal_num = sum(grepl('normal', unique(overlap_gr$sample)))
    # This line will average multiple scores from the same sample and split normal and tumor scores by column
    if (length(overlap_gr) == 0) {
        mean_tumor_score = 0; mean_normal_score = 0; mean_overall_score = 0
    } else {
        # Tumor score and normal score
        score_df = reshape2::dcast(data.frame(mcols(overlap_gr)), sample ~ sample_type, value.var = 'score', 
                                    fun.aggregate = mean, fill = 0)
        mean_tumor_score = sum(score_df$tumor, na.rm = TRUE) / total_tumor_sample
        mean_normal_score = sum(score_df$normal, na.rm = TRUE) / total_normal_sample
        # Overall score
        score_df2 = reshape2::dcast(data.frame(mcols(overlap_gr)), sample ~ rep, value.var = 'score', 
                                fun.aggregate = mean, fill = 0)
        mean_overall_score = sum(score_df2$rep1) / (total_tumor_sample + total_normal_sample)
    }
    combined_metrics_df[i, ] = c(i, length(overlap_gr), tumor_num, normal_num, 
                                    mean_tumor_score, mean_normal_score, mean_overall_score)
    pb$tick()
}
pb$terminate()

# Save
saveRDS(combined_metrics_df, file = paste0(plots_path, '/ARBS_combined_metrics_df.rds'))
cat('Site binding score annotation done\n')

combined_metrics_df = readRDS(paste0(plots_path, '/ARBS_combined_metrics_df.rds'))
# Add combined metrics from combined_metrics_df as metadata columns to merged_gr
mcols(merged_gr) <- cbind(mcols(merged_gr), combined_metrics_df[, -1])
saveRDS(merged_gr, file = paste0(plots_path, '/ARBS_merged_gr.rds'))
arbs_gr = merged_gr


# --- Load SV breakpoints ---
cat('Loading SVs...')
sv_combined = read.csv(paste0(ppcg_dir, '/raw/SV/SV_driver_MS.csv'), stringsAsFactors = FALSE)
sv_gr = mut_df2gr(sv_combined, no_seq = TRUE, add_col = TRUE)
cat('Done\n')


# --- Transcription start sites of genes ---
cat('Loading genome annotations...')
gene_gr = get_annotations(src = 'grch37')
gene_gr = gene_gr[gene_gr$type == 'gene'] # Select main transcript
mcols(gene_gr)[, 'gene_name'] = gsub('^gene-', '', gene_gr$site)
tss_gr = IRanges::promoters(gene_gr, upstream = dist_threshold, downstream = dist_threshold)
cat('Done\n')


# --- Annotate all genes ---
# Initialize additional columns in annot_gr
annot_gr = tss_gr
new_cols <- list(
    sv_count               = rep(0, length(annot_gr)),
    sv_samples             = rep(0, length(annot_gr)),
    sv_class               = rep('', length(annot_gr)),
    n_sites_in_range       = rep(0, length(annot_gr)),
    n_normal               = rep(0, length(annot_gr)),
    n_tumor                = rep(0, length(annot_gr)),
    n_both                 = rep(0, length(annot_gr)),
    mean_tumor_score       = rep(0, length(annot_gr)),
    mean_normal_score      = rep(0, length(annot_gr)),
    mean_overall_score     = rep(0, length(annot_gr)),
    max_tumor_score        = rep(0, length(annot_gr)),
    max_normal_score       = rep(0, length(annot_gr)),
    max_overall_score      = rep(0, length(annot_gr))
)
mcols(annot_gr) <- cbind(mcols(annot_gr), new_cols)
# Subset ranges of interest for more efficient annotation
overlapped_sv_gr = subsetByOverlaps(sv_gr, annot_gr)
overlapped_arbs_gr = subsetByOverlaps(arbs_gr, annot_gr)
#
annot_gr$elem_index = 1:length(annot_gr)
gr_with_svs = subsetByOverlaps(annot_gr, sv_gr)
gr_with_arbs = subsetByOverlaps(annot_gr, arbs_gr)
index_to_annotate = unique(c(gr_with_svs$elem_index, gr_with_arbs$elem_index))
#
pb = progress_bar_init(length(index_to_annotate), msg = 'Annotating gene flanks')
for (i in index_to_annotate) {
    this_gr = annot_gr[i]
    # Find overlapping SVs
    sv_count = countOverlaps(this_gr, overlapped_sv_gr)
    if (sv_count > 0) {
        mcols(annot_gr)[i, 'sv_count'] = sv_count
        # Number of samples with SVs
        sv_hits <- findOverlaps(this_gr, overlapped_sv_gr)
        sv_samples <- length(unique(overlapped_sv_gr$samples[subjectHits(sv_hits)]))
        mcols(annot_gr)[i, 'sv_samples'] <- sv_samples
        sv_table <- table(overlapped_sv_gr$svclass[subjectHits(sv_hits)])
        sv_table_str <- paste(sapply(names(sv_table), function(x) paste(x, sv_table[x], sep=':')), collapse=',')
        mcols(annot_gr)[i, 'sv_class'] <- sv_table_str
    }
    # Find overlapping ARBS
    arbs_count = countOverlaps(this_gr, overlapped_arbs_gr)
    if (arbs_count > 0) {
        mcols(annot_gr)[i, 'n_sites_in_range'] = arbs_count
        arbs_in_range = subsetByOverlaps(overlapped_arbs_gr, this_gr)
        mcols(annot_gr)[i, 'n_normal'] = sum(arbs_in_range$annotation == 'narbs')
        mcols(annot_gr)[i, 'n_tumor'] = sum(arbs_in_range$annotation == 'tarbs')
        mcols(annot_gr)[i, 'n_both'] = sum(arbs_in_range$annotation == 'both')
        mcols(annot_gr)[i, 'mean_tumor_score'] = mean(arbs_in_range$mean_tumor_score)
        mcols(annot_gr)[i, 'mean_normal_score'] = mean(arbs_in_range$mean_normal_score)
        mcols(annot_gr)[i, 'mean_overall_score'] = mean(arbs_in_range$mean_overall_score)
        mcols(annot_gr)[i, 'max_tumor_score'] = max(arbs_in_range$mean_tumor_score)
        mcols(annot_gr)[i, 'max_normal_score'] = max(arbs_in_range$mean_normal_score)
        mcols(annot_gr)[i, 'max_overall_score'] = max(arbs_in_range$mean_overall_score)
    }
    pb$tick()
}
pb$terminate()
saveRDS(annot_gr, file = paste0(plots_path, '/ARBS_annotated_TSS_gr.rds'))
cat('Gene binding annotation done\n')

