#' Check corrleation between EPHA1-AS1/EPHA1 and NEAT1/MALAT1
source('FigureS3_utils.R')

#
plots_path = new.dir(paste0(plot_dir, '/ppcg/drivers_MS_repo'))

# --- Load RNA data ---
load(paste0(ppcg_dir, '/raw/RNA/expr_mat_20250106.rsav')) 
rownames(expr_mat) = gsub('\\-', '_', rownames(expr_mat))
rna_df = as.data.frame(t(expr_mat))


# --- Plot gene correlation ---
pair_list = list(
    c('NEAT1', 'MALAT1'), 
    c('EPHA1_AS1', 'EPHA1')
)
#
pdf(paste0(plots_path, '/gene_correlations.pdf'), width = 4, height = 4)
for (gene_pair in pair_list) {
    gene_x = gene_pair[1]; gene_y = gene_pair[2]
    p = corr_scat_plot(rna_df, x = gene_x, y = gene_y) 
    print(p)
}
dev.off()


