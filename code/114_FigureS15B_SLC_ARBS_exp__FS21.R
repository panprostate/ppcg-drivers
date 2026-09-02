#' Associating SLC genes ARBS with expression 
source('FigureS3_utils.R')
library(ggplot2)
library(ggrepel)
library(dplyr)
library(scales)


plots_path = paste0(plot_dir, '/ppcg/drivers_MS_repo')


# --- Load annotated gene_gr ---
gene_gr <- readRDS(paste0(plots_path, '/ARBS_annotated_TSS_gr.rds'))
gene_df <- mut_gr2df(gene_gr, no_seq = TRUE, add_col = TRUE)
gene_df <- gene_df[order(gene_df$n_sites_in_range, decreasing = TRUE), ]
gene_df$has_arbs = (gene_df$n_sites_in_range > 0)


# --- Add mean PPCG expresssion to gene_df ---
load(paste0(ppcg_dir, '/raw/RNA/expr_mat_20250106.rsav')) 
rna_df = as.data.frame(t(expr_mat))
mean_expr = rowMeans(expr_mat)
gene_df$mean_expr = mean_expr[gene_df$gene_name]
gene_df = gene_df[!is.na(gene_df$mean_expr), ]


# --- Compare expression between SLC gene w/wo ARBS ---
# Empirical p-value with ECDF
slc_gene_df = gene_df[grepl('SLC', gene_df$gene_name), ]
med_expr = median(slc_gene_df$mean_expr[slc_gene_df$n_sites_in_range == 0])
ecdf_expr = ecdf(slc_gene_df$mean_expr[slc_gene_df$n_sites_in_range == 0])
p_values = sapply(slc_gene_df$mean_expr, function(g_value) {
    1 - ecdf_expr(g_value)
})
slc_gene_df$p_value = p_values
# Test on genes with TSS ARBS only 
slc_gene_df$p_value[!slc_gene_df$has_arbs] = NA
slc_gene_df$p_0.05[!slc_gene_df$has_arbs] = NA
# Define AR-regulated genes
slc_gene_df$p_0.05 = (slc_gene_df$p_value < 0.05)
slc_gene_df$ar_regulated = (slc_gene_df$has_arbs & slc_gene_df$p_0.05)
# Save SLC table
slc_table_df = slc_gene_df[c('gene_name', 'n_sites_in_range', 'mean_expr', 'p_value', 
                                'has_arbs', 'sv_count')]
slc_table_df = slc_table_df[order(-slc_table_df$has_arbs, -slc_table_df$mean_expr, slc_table_df$p_value), ]
write.csv(slc_table_df, file = paste0(plots_path, '/slc_table.csv'), row.names = FALSE)


# --- Plot expression comparison ---
plot_df1 = slc_gene_df
# Highlight putative AR-regulated SLC genes
ar_regulated_slc = slc_table_df$gene_name[slc_table_df$ar_regulated]
plot_df1$highlight_color = ifelse(plot_df1$gene_name %in% ar_regulated_slc, 'highlight', 'normal')
# Label highly recurrent SV SLC genes
label_genes = c('SLC30A4', 'SLC45A3')
plot_df1$is_label = plot_df1$gene_name %in% label_genes
plot_df1$label_name = ifelse(plot_df1$is_label, as.character(plot_df1$gene_name), NA)
# Group Wilcox test
x = 'has_arbs'; y = 'mean_expr'
wilcox_result = wilcox.test(as.formula(paste0(y, ' ~ ', x)), data = plot_df1)
#
pdf(paste0(plots_path, '/SLC_ARBS_boxplot.pdf'), width = 3, height = 4)
p = ggplot(plot_df1, aes(x = .data[[x]], y = .data[[y]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = highlight_color, shape = is_label), width = 0.2, alpha = 0.6, show.legend = FALSE) +
    geom_text_repel(aes(label = label_name), size = 2, 
                    point.padding = unit(0.5, 'lines')) + 
    scale_color_manual(values = c('highlight' = 'red', 'normal' = 'black')) +
    labs(title = paste0('SLC genes\n(Wilcox p = ', 
                    dec_format(wilcox_result$p.value), ')'),
            x = var_format(x), y = var_format(y)) +
    theme_bw()
print(p)
dev.off()


