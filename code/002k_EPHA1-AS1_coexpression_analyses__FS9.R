source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(ggrepel)

CUTOFF = 0.05
N_PAIRS = 1e6



set.seed(12345)


load(pff("cgc2024.rsav"))
load(pff("prostate_cancer_genes.rsav"))
cancer_genes = c(cgc2024, prostate_cancer_genes)

load(file = pff("expr_mat_with_LNCs.rsav"))


# generate random gene pairs
all_genes = rownames(expr_mat)
g1 = sample(all_genes, N_PAIRS, replace = T)
g2 = sample(all_genes, N_PAIRS, replace = T)

g_pairs = cbind(g1, g2)
g_pairs = g_pairs[g_pairs[,1] != g_pairs[,2],]

# i = 1 
correlate_pairs = function(i, g_pairs, expr_mat, type = "bg") {
	
	if (i %% 1000 == 0) cat("O")
	
	g1 = g_pairs[i,1]
	g2 = g_pairs[i,2]
	
	g1_val = expr_mat[g1,, drop = T]
	g2_val = expr_mat[g2,, drop = T]
	
	
	cor_test = cor.test(g1_val, g2_val, method = "spearman")
	cor_val = cor_test$estimate[[1]]
	pval = cor_test$p.value
	res = cor_val
	
	if (type != "bg") {
		res = data.frame(g1, g2, cor = res, pval, stringsAsFactors = FALSE)
	}
	return(res)
}


corrs_bg = unlist(mclapply(1:nrow(g_pairs), correlate_pairs, g_pairs, expr_mat, mc.cores = 8))


lower_quantile = CUTOFF / 2
higher_quantile = 1 - CUTOFF / 2

lower_cut = quantile(corrs_bg, lower_quantile)
higher_cut = quantile(corrs_bg, higher_quantile)

cancer_genes_here = intersect(cancer_genes, rownames(expr_mat))
pairs_to_test = cbind("EPHA1-AS1", cancer_genes_here)

corrs_tests = do.call(rbind, mclapply(1:nrow(pairs_to_test), correlate_pairs, pairs_to_test, expr_mat, type = "test", mc.cores = 8))
corrs_tests$select = "no"
corrs_tests[corrs_tests$cor < lower_cut, "select"] = "lower"
corrs_tests[corrs_tests$cor >= higher_cut, "select"]  = "higher"

corrs_tests_select = corrs_tests[corrs_tests$select %in% c("lower", "higher"),]
corrs_tests_select$g2 = factor(corrs_tests_select$g2, levels = corrs_tests_select$g2[order(corrs_tests_select$cor)])

corrs_tests_select$is_prca_gene = c("CGC", "PrCa")[1 + as.character(corrs_tests_select$g2) %in% prostate_cancer_genes]
corrs_tests_select$label = NA
corrs_tests_select$label[corrs_tests_select$is_prca_gene == "PrCa"] = as.character(corrs_tests_select$g2)[corrs_tests_select$is_prca_gene == "PrCa"]

titl = paste(nrow(corrs_tests_select), " cancer genes, cor with EPHA1-AS1\nPermutation P < ", CUTOFF)

plt = ggplot(corrs_tests_select, aes (g2, cor, fill = is_prca_gene, label = label)) +
		geom_bar(stat = "identity")  +
		geom_text_repel(y = 0, angle = 45, size = 4, max.overlaps = nrow(corrs_tests_select), min.segment.length = 0) + 
		scale_fill_manual(values = c("PrCa" = "darkred", "CGC" = "darkorange")) + 
		scale_y_continuous("Spearman rho") +
		plot_theme() + 
		theme(axis.text.x = element_text(size = 5), legend.position = "bottom") + 
		ggtitle(NULL, titl)

fname = pff("figures/EPHA1AS1_significant_correlations.pdf")
ggsave(plt, file = fname, width = 8, height = 5)
file_open_call2(fname)

