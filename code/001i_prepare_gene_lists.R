source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")


#
# other gene lists for interpretation
#
load(file = pff("gencode_genes.rsav"))
gencode_symbols = unique(gsub("(.+)::(.+)::(.+)::(.+)", "\\3", gencode_genes$id))

## cancer gene census
cgc2024 = read.delim(
		paste0("DATA_USED__", this_timestamp, "/cancer_gene_lists_Juri_2024-09-18/cgc2024/Census_allWed_Sep_18_18_01_40_2024.tsv"), 
		stringsAsFactors = FALSE)
cgc2024 = unique(toupper(cgc2024$Gene.Symbol))

missing_cgc_genes = cgc2024[!cgc2024 %in% gencode_symbols]
cgc2024[cgc2024 %in% "AFDN"] = "MLLT4"
cgc2024[cgc2024 %in% "LHFPL6"] = "LHFP"
cgc2024[cgc2024 %in% "MRTFA"] = "MKL1"
cgc2024[cgc2024 %in% "NSD2"] = "WHSC1"
cgc2024[cgc2024 %in% "NSD3"] = "WHSC1L1"
cgc2024[cgc2024 %in% "SHTN1"] = "KIAA1598"
cgc2024[cgc2024 %in% "TENT5C"] = "FAM46C"
cgc2024[cgc2024 %in% "WDCP"] = "C2orf44"

cgc2024 = setdiff(cgc2024, missing_cgc_genes)
save(cgc2024, file = pff("cgc2024.rsav"))



##
# prostate cancer genes
##
wedge_drivers = read.csv(
		paste0("DATA_USED__", this_timestamp, "/cancer_gene_lists_Juri_2024-09-18/Wedge2018/Wedge2018_genes.csv"), 
		stringsAsFactors = FALSE)
wedge_drivers = setdiff(unique(wedge_drivers$gene), "")

missing_wedge_genes = setdiff(wedge_drivers, gencode_symbols)
wedge_drivers[wedge_drivers == "MLL2"] = "KMT2D"
wedge_drivers[wedge_drivers == "MLL3"] = "KMT2C"
wedge_drivers[wedge_drivers == "MYST3"] = "KAT6A"
save(wedge_drivers, file = pff("wedge_drivers.rsav"))

armenia_drivers = read.csv(
		paste0("DATA_USED__", this_timestamp, "/cancer_gene_lists_Juri_2024-09-18/Armenia2018/Armenia2018_genes.csv"), 
		stringsAsFactors = FALSE)
armenia_drivers = setdiff(unique(armenia_drivers[armenia_drivers$final_list.n.97.total. == 1,1]), "")
missing_armenia_genes = setdiff(armenia_drivers, gencode_symbols)
armenia_drivers [armenia_drivers == "JADE2"] = "PHF15"
save(armenia_drivers, file = pff("armenia_drivers.rsav"))

prostate_cancer_genes = unique(c(wedge_drivers, armenia_drivers))
save(prostate_cancer_genes, file = pff("prostate_cancer_genes.rsav"))

##
# GMT file for pathway analyses
##
gobp = readLines(
		paste0("DATA_USED__", this_timestamp, 
		"/cancer_gene_lists_Juri_2024-09-18/gene_sets/gene_sets_gprofiler_2024-09-17/gprofiler_hsapiens.name/hsapiens.GO:BP.name.gmt"))
reac = readLines(
		paste0("DATA_USED__", this_timestamp, 
		"/cancer_gene_lists_Juri_2024-09-18/gene_sets/gene_sets_gprofiler_2024-09-17/gprofiler_hsapiens.name/hsapiens.REAC.name.gmt"))

gobp_reac_gmt = c(gobp, reac)
writeLines(gobp_reac_gmt, pff("gobp_reac.gmt"))
writeLines(reac, pff("reac.gmt"))