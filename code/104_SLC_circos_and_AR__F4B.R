# Figure 4B
# Load libraries
library(dplyr)
library(data.table)
library(circlize)

# Set working directory
setwd("/CESAM2_PPCG/")

# Load files
# Load overlap file between TADs and SVs
overlap <- fread("TADS_overlap_SVS_pairtobed_2025_02_10.bed", header = TRUE, sep = "\t") %>%
  unique() %>%
  mutate(TAD = paste0("chr", chromTAD, ":", startTAD, "-", endTAD))

# Load PCGC cohort data (only samples that have SVs)
cohort <- fread("WGS_master_tracking_sheet.tsv", na.strings = c("", "NA")) %>%
  filter(selected_one_sample_per_donor == TRUE & SV == TRUE & SampleType == "TUM")

# Filter overlap data to include only samples in the cohort
overlap <- overlap %>%
  filter(samples %in% cohort$WGS_AssayID)

# Load ARB data
arb <- read.table("slc_table.csv", header = TRUE, sep = ",")

# Define oncogenes of interest
oncogenes <- c("ERG", "ETV4", "ETV1", "FLI1", "SKIL", "MYC", "ETV5", "AR", "ETS1", "ELK4")

# Extract SLCs from overlap data
SLCs <- overlap %>%
  filter(grepl("^SLC", name_gene1) | grepl("^SLC", name_gene2)) %>%
  mutate(SLC = ifelse(grepl("^SLC", name_gene1), name_gene1, name_gene2)) %>%
  filter(!is.na(SLC))

# Extract SLCs that are associated with oncogenes
SLCs_oncogene <- SLCs %>%
  filter(name_gene1 %in% oncogenes | name_gene2 %in% oncogenes) %>%
  mutate(oncogene = ifelse(name_gene1 %in% oncogenes, name_gene1, name_gene2)) %>%
  select(chrom1, start1, end1, chrom2, start2, end2, SLC, oncogene, samples) %>%
  mutate(Index = row_number())

# Summarize the number of oncogene pairs
oncogene_pairs <- SLCs_oncogene %>%
  group_by(oncogene) %>%
  summarise(n = n())

# Define genome version
genome <- "hg19"

# Ensure chromosome names have 'chr' prefix
SLCs_oncogene <- SLCs_oncogene %>%
  mutate(chrom1 = paste0("chr", chrom1), chrom2 = paste0("chr", chrom2))

# Convert positions to numeric
SLCs_oncogene <- SLCs_oncogene %>%
  mutate(across(c(start1, start2), as.numeric))

# Get unique SLCs for dots
SLCs_oncogene_unique <- SLCs_oncogene %>%
  distinct(SLC, .keep_all = TRUE)

# Prepare oncogenes bed data
oncogenes_bed <- SLCs %>%
  filter(name_gene1 %in% oncogenes | name_gene2 %in% oncogenes) %>%
  mutate(
    chrom = ifelse(name_gene1 %in% oncogenes, chrom_gene1, chrom_gene2),
    start = ifelse(name_gene1 %in% oncogenes, start_gene1, start_gene2),
    end = ifelse(name_gene1 %in% oncogenes, end_gene1, end_gene2),
    gene = ifelse(name_gene1 %in% oncogenes, name_gene1, name_gene2)
  ) %>%
  select(chrom, start, end, gene) %>%
  distinct()

# Prepare SLC bed data
SLC_bed <- SLCs %>%
  filter(name_gene1 %in% SLCs_oncogene_unique$SLC | name_gene2 %in% SLCs_oncogene_unique$SLC) %>%
  mutate(
    chrom = ifelse(name_gene1 %in% SLCs_oncogene_unique$SLC, chrom_gene1, chrom_gene2),
    start = ifelse(name_gene1 %in% SLCs_oncogene_unique$SLC, start_gene1, start_gene2),
    end = ifelse(name_gene1 %in% SLCs_oncogene_unique$SLC, end_gene1, end_gene2),
    gene = ifelse(name_gene1 %in% SLCs_oncogene_unique$SLC, name_gene1, name_gene2)
  ) %>%
  select(chrom, start, end, gene) %>%
  distinct()

# Annotate ARB and Type for SLC bed data
SLC_bed <- SLC_bed %>%
  mutate(
    ARB = ifelse(gene %in% arb[arb$ar_regulated == TRUE,]$gene_name, "ARB", NA),
    Type = "SLC"
  )

# Annotate ARB and Type for oncogenes bed data
oncogenes_bed <- oncogenes_bed %>%
  mutate(
    ARB = ifelse(gene %in% arb[arb$ar_regulated == TRUE,]$gene_name, "ARB", NA),
    Type = "oncogene"
  )

# Combine SLC and oncogenes bed data
all <- bind_rows(oncogenes_bed, SLC_bed) %>%
  mutate(chrom = paste0("chr", chrom)) %>%
  data.frame()

# Figure 4B
# Initialize circos plot
circos.clear()
circos.initializeWithIdeogram(species = genome, chromosome.index = unique(all$chrom), plotType = NULL)

# Add genomic labels to circos plot
circos.genomicLabels(all, labels.column = 4, side = "outside", niceFacing = TRUE, cex = 0.6, col = as.numeric(factor(all$Type)))

# Add genomic ideogram to circos plot
circos.genomicIdeogram(species = genome)

# Define annotation colors
annotation_colors <- c("ARB" = "#ffbb44")

# Add midpoint to all data
all <- all %>%
  mutate(midpoint = (start + end) / 2)

# Add genomic track for ARB to circos plot
circos.genomicTrack(all[all$ARB == "ARB", ], panel.fun = function(region, value, ...) {
  circos.genomicPoints(
    region = region,
    value = data.frame(y = rep(0.25, nrow(region))),
    col = annotation_colors["ARB"],
    pch = 16,
    cex = 0.8,
    ...
  )
}, ylim = c(0, 0.5), track.height = 0.05, bg.border = NA)

# Add links between SLCs and oncogenes to circos plot
for (i in seq_len(nrow(SLCs_oncogene))) {
  circos.link(SLCs_oncogene$chrom1[i], SLCs_oncogene$start1[i], 
              SLCs_oncogene$chrom2[i], SLCs_oncogene$start2[i])
}

# Load gene expression data
gene_expression <- fread("PPCG_RNAseqData_Tier1.ruv.Var20Filtered_sorted.bed.gz")
#setDT(gene_expression)

# Filter gene expression data for oncogenes
gene_expression <- gene_expression %>%
  filter(symbol %in% SLCs_oncogene$oncogene)

# Melt gene expression data for merging
gene_expression_m <- melt(gene_expression[, c(5, 7:507)])
gene_expression_m <- unique(merge(gene_expression_m, SLCs_oncogene, by.x = c("symbol", "variable"), by.y = c("oncogene", "samples"), all.x = TRUE))
gene_expression_m$SV <- ifelse(is.na(gene_expression_m$chrom1), "SV-", "SV+")

# Calculate median values for SV- samples
median_values <- gene_expression_m %>%
  filter(SV == "SV-") %>%
  group_by(symbol) %>%
  summarise(median_value = median(value, na.rm = TRUE), .groups = "drop")

# Calculate fold change for SV+ samples
df <- gene_expression_m %>%
  left_join(median_values, by = "symbol") %>%
  mutate(fold_change = ifelse(SV == "SV+", value / median_value, NA))

# Supplementary table
# Filter and select relevant columns for final data frame
supp_tab_figure4B <- df %>%
  filter(SV == "SV+") %>%
  select(symbol, SLC, variable, chrom1, start1, end1, chrom2, start2, end2, fold_change)%>%
  dplyr::rename(oncogene = symbol, Sample = variable)%>%
  unique()

