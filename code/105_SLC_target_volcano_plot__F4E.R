# Figure 4E
# Load libraries
library(dplyr)
library(data.table)
library(ggpubr)
library(ggrepel)
library(tidyverse)

# Set working directory
setwd("/CESAM_PPCG/")

# Load files
# Load overlap file between TADs and SVs
overlap <- fread("TADS_overlap_SVS_peirtobed_2025_02_10.bed", header = TRUE, sep = "\t") %>%
  unique() %>%
  mutate(TAD = paste0("chr", chromTAD, ":", startTAD, "-", endTAD))

# Load CESAM2 results
permutations <- fread("results_2025_02_10/permutation_intraTAD_2025_02_10.txt", header = FALSE, sep = " ") %>%
  filter(!is.na(V20)) %>%
  arrange(V20)

# Load genotype matrix
genotype <- fread("all_matrix.final_new2_sorted_abc_0kb.vcf.gz")

# Load PCGC cohort (only samples that we do have SVs)
cohort <- fread("WGS_master_tracking_sheet.tsv", na.strings = c("", "NA")) %>%
  filter(selected_one_sample_per_donor == TRUE & SV == TRUE & SampleType == "TUM")

# Filter overlap data to include only samples in the cohort
overlap <- overlap %>%
  filter(samples %in% cohort$WGS_AssayID)

# Convert genotype data to a data.table
setDT(genotype)

# Exclude non-genotype columns to count genotypes across sample columns
genotype_cols <- setdiff(names(genotype), c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))

# Efficiently count occurrences of each genotype (0/0, 0/1, 1/1) for each ID
genotype_counts <- genotype[, .(
  count_0_0 = sum(unlist(.SD) == "0/0"),
  count_0_1 = sum(unlist(.SD) == "0/1"),
  count_1_1 = sum(unlist(.SD) == "1/1")
), by = ID, .SDcols = genotype_cols]

# Load gene expression data
gene_expr <- fread("PPCG_RNAseqData_Tier1.ruv.Var20Filtered_sorted.bed.gz")
setDT(gene_expr)

# Find common samples between gene expression and genotype data
common_samples <- intersect(names(gene_expr), names(genotype))

# Keep only the relevant genotype columns that exist in both data frames
genotype_cols <- common_samples

# Efficiently count occurrences of each genotype (0/0, 0/1, 1/1) for each ID in matched samples
genotype_counts_matched <- genotype[, .(
  count_0_0_WGS_RNA = sum(unlist(.SD) == "0/0"),
  count_0_1_WGS_RNA = sum(unlist(.SD) == "0/1"),
  count_1_1_WGS_RNA = sum(unlist(.SD) == "1/1")
), by = ID, .SDcols = genotype_cols]

# Merge overlap data with permutations
overlap <- merge(
  overlap,
  permutations[, c("V1", "V7", "V8", "V18", "V20")],
  by.x = "TAD",
  by.y = "V8",
  all.x = TRUE,
  allow.cartesian = TRUE  # Allows multiple matches per key
)
colnames(overlap)[colnames(overlap) %in% c("V1", "V7", "V18", "V20")] <- c("CESAM_gene", "CESAM_distance", "CESAM_SizeEffect", "CESAM_p.adj")

# Extract SLCs
SLCs <- overlap %>%
  filter(grepl("^SLC", name_gene2) | grepl("^SLC", name_gene1)) %>%
  select(TAD, name_gene1, name_gene2) %>%
  mutate(SLC = ifelse(grepl("^SLC", name_gene1), name_gene1, ifelse(grepl("^SLC", name_gene2), name_gene2, NA))) %>%
  select(TAD, SLC) %>%
  distinct()

# Merge permutations with SLCs
permutations_SLC <- merge(permutations, SLCs, by.x = "V8", by.y = "TAD", all.x = TRUE) %>%
  select(V8, V1, V7, V18, V20, SLC) %>%
  dplyr::rename(TAD = V8, gene_name = V1, distance = V7, EffectSize = V18, p.adj = V20)

# Merge TADs info with gene expression data
gene_expr_TADs <- merge(gene_expr, permutations[, c("V1", "V8")], by.x = "symbol", by.y = "V1", all.x = TRUE)

# Melt gene expression data for merging
gene_expr_TADs_m <- reshape2::melt(gene_expr_TADs[, c(1, 508, 7:507)])

# Melt genotype data directly using data.table's melt function
sample_cols <- intersect(colnames(gene_expr_TADs), colnames(genotype))
selected_cols <- c("ID", sample_cols)

genotype_selected <- melt(genotype[, ..selected_cols], id.vars = "ID", variable.name = "Sample", value.name = "Genotype")

# Merge genotype and gene expression data
merged_data <- merge(gene_expr_TADs_m, genotype_selected, by.x = c("V8", "variable"), by.y = c("ID", "Sample")) %>%
  mutate(Genotype2 = ifelse(Genotype == "0/1", "1/1", ifelse(Genotype == "0/0", "0/0", "1/1")))

# Calculate fold change
fold_change <- merged_data %>%
  group_by(symbol, V8, Genotype2) %>%
  summarise(mean = mean(value))

fold_change_dcast <- reshape2::dcast(fold_change, symbol + V8 ~ Genotype2, value.var = "mean") %>%
  mutate(log2_fc = `1/1` - `0/0`)

# Combine logFC with permutation results
permutations_SLC_logFC <- merge(permutations_SLC, fold_change_dcast, by.x = c("TAD", "gene_name"), by.y = c("V8", "symbol"), all.x = TRUE)

# Combine with overlap
overlap <- merge(overlap, fold_change_dcast[, c("symbol", "log2_fc")], by.x = "CESAM_gene", by.y = "symbol", allow.cartesian = TRUE)

# Merge with cohort data
overlap <- merge(overlap, cohort[, c("WGS_AssayID", "selected_one_sample_per_donor", "SampleType")], by.x = "samples", by.y = "WGS_AssayID", all.x = TRUE)

# Calculate size and filter data
overlap <- overlap %>%
  mutate(Size = abs(start2 - start1)) %>%
  filter((chrom1 != chrom2 | Size > 1000000) & (abs(dist_gene1) < 200000 & abs(dist_gene2) < 200000)) %>%
  distinct(samples, SampleType, selected_one_sample_per_donor, name_gene1, name_gene2) %>%
  filter(name_gene1 != name_gene2) %>%
  distinct()

# Define genes of interest
genes_of_interest <- c("ERG", "ETV1", "ETV4", "ETV5", "FLI1", "MYC", "AR", "ETS1")

# Filter data for genes of interest
filtered_data <- overlap %>%
  filter(selected_one_sample_per_donor == TRUE & SampleType == "TUM") %>%
  filter(name_gene1 %in% genes_of_interest | name_gene2 %in% genes_of_interest)

# Volcano plot
permutations_SLC_logFC_filtered <- permutations_SLC_logFC %>%
  filter(!is.na(log2_fc)) %>%
  mutate(
    log10P = -log10(p.adj),
    SLC_present = ifelse(is.na(SLC), 0, 1),
    color_group = ifelse(log10P > 2 & log2_fc >= 1 & SLC_present == 1, "darkblue", "grey")
  ) %>%
  distinct()

# Volcano plot only SLCs
permutations_SLC_logFC_filtered_2 <- permutations_SLC_logFC_filtered %>%
  group_by(gene_name) %>%
  mutate(SLC_count = n_distinct(SLC)) %>%  # Count unique SLCs for each gene_name
  ungroup() %>%
  mutate(color_group = ifelse(log10P > 2 & log2_fc >= 0.5 & SLC_present == 1, "darkblue", "grey"))

label_data <- permutations_SLC_logFC_filtered_2 %>%
  filter(log10P >= 10 & log2_fc >= 0.5 & SLC_present == 1) %>%
  distinct(gene_name, .keep_all = TRUE)

# Define min-max scaling function
scale_log10P <- function(x, new_min = 0, new_max = 20) {
  old_min <- min(x, na.rm = TRUE)
  old_max <- max(x, na.rm = TRUE)
  scaled_x <- new_min + (x - old_min) * (new_max - new_min) / (old_max - old_min)
  return(scaled_x)
}

# Apply scaling
permutations_SLC_logFC_filtered_2 <- permutations_SLC_logFC_filtered_2 %>%
  mutate(scaled_log10P = scale_log10P(log10P, 0, 20))

permutations_SLC_logFC_filtered_3 <- permutations_SLC_logFC_filtered_2 %>%
  mutate(
    capped_log10P = pmin(log10P, 60),
    capped_size_scaled = scales::rescale(capped_log10P, to = c(2, 10))
  ) %>%
  distinct(gene_name, log2_fc, capped_log10P, SLC_present, color_group)

label_data <- permutations_SLC_logFC_filtered_3 %>%
  filter(capped_log10P >= 10 & log2_fc >= 0.5 & SLC_present == 1) %>%
  distinct(gene_name, .keep_all = TRUE)

# Plot volcano plot for SLCs
volcano <- ggplot(unique(subset(permutations_SLC_logFC_filtered_3, SLC_present == 1))) +
  geom_point(aes(x = log2_fc, y = capped_log10P), color = "grey", alpha = 0.5) +
  geom_point(data = subset(permutations_SLC_logFC_filtered_3, color_group == "darkblue"), aes(x = log2_fc, y = capped_log10P), color = "darkblue") +
  theme_pubr() +
  geom_text_repel(
    data = label_data,
    aes(x = log2_fc, y = capped_log10P, label = gene_name),
    show.legend = FALSE,
    min.segment.length = 0,
    max.overlaps = 15,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = "darkgrey",
    size = 7 / .pt,
    family = "Arial"  # Set font to Arial
  ) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "grey40") +
  labs(
    title = "SLC target genes",
    subtitle = "Blue dots: Genes with at least one SLC with -log10P > 2 & log(fc) ≥ 0.5",
    x = expression(paste(Log, " fold gene expression change (SV+/SV-)", sep = "")),
    y = expression(paste(-Log[10], italic("P"), sep = ""))
  ) +
  theme(
    text = element_text(family = "Arial", size = 7),  # Set global text to Arial and size 7
    plot.title = element_text(size = 7, face = "bold"),
    plot.subtitle = element_text(size = 7),
    plot.caption = element_text(size = 7, face = "italic"),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linetype = "dashed")
  ) +
  xlim(-3.5, 4.5)

# Print the volcano plot
print(volcano)
