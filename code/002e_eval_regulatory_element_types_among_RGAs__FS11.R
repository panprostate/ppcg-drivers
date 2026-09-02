source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library("gtools")


# sample elements directly from ADWGS results, measure how many elements of each type captured 
total_times_sampled = 100000

get_element_type = function(res_dfr, elements_having_ENSG) {
	element_type2 = gsub("(.+)::(.+)::(.+)::(.+)", "\\1__\\4", res_dfr$id)
	which_elements_to_remove_ensg = which(res_dfr$element_type %in% elements_having_ENSG)
	element_type2[which_elements_to_remove_ensg] = res_dfr[which_elements_to_remove_ensg, "element_type"]
	element_type2		
}

get_expected_element_types = function(results_all, n_iters, n_to_sample, all_element_types, mc.cores) {
	res = do.call(rbind, mclapply(1:n_iters, get_sample, results_all, all_element_types, n_to_sample, mc.cores = mc.cores))
	colnames(res) = all_element_types
	res[is.na(res)] = 0
	res
}

get_sample = function(i, results_all, all_element_types, n_to_sample) {
	if (runif(1) < 0.001) cat("O")
	res = c(table(sample(results_all$element_type2, n_to_sample))[all_element_types])
	res
}

get_empirical_p = function(element_type, observed_element_stats, sampled_element_types) {
	obs_val = observed_element_stats[observed_element_stats$element_type == element_type, "median_val"]
	p_higher = (1 + sum(sampled_element_types[,element_type] >= obs_val)) / (1 + nrow(sampled_element_types))
	p_lower = (1 + sum(sampled_element_types[,element_type] < obs_val)) / (1 + nrow(sampled_element_types))
	data.frame(element_type, p_higher, p_lower, stringsAsFactors = FALSE)
}

load(pff("results_signf.rsav"))
elements_having_ENSG = c("UTR3", "UTR5", "lncRNA")
results_signf = results_signf[!results_signf$element_type %in% c("CDS", "CDSgene"),]
results_signf$element_type2 = get_element_type(results_signf, elements_having_ENSG)

load(pff("results_ADWGS.rsav"))
results_all = results[!results$element_type %in% c("CDS", "CDSgene"),]
results_all$element_type2 = get_element_type(results_all, elements_having_ENSG)
n_to_sample = nrow(results_signf)
all_element_types = unique(results_all$element_type2)

# get observed / expected values from sampling 
sampled_element_types = get_expected_element_types(results_all, total_times_sampled, n_to_sample, all_element_types, mc.cores = 8)
expected_element_stats = data.frame(
		element_type = colnames(sampled_element_types),
		median_val = apply(sampled_element_types, 2, median), 
		bottom_val = apply(sampled_element_types, 2, quantile, 0.025, na.rm = T), 
		top_val = apply(sampled_element_types, 2, quantile, 0.975, na.rm = T), 
		value_type = "expected",
		stringsAsFactors = FALSE
)
obs_vals = c(table(results_signf$element_type2)[expected_element_stats$element_type])
obs_vals[is.na(obs_vals)] = 0
names(obs_vals) = expected_element_stats$element_type
observed_element_stats = data.frame(
		element_type = expected_element_stats$element_type,
		median_val = obs_vals, 
		bottom_val = NA,
		top_val = NA,
		value_type = "observed",
		stringsAsFactors = FALSE
)

combined_element_stats = rbind(expected_element_stats, observed_element_stats)
sampled_pvals_stats = do.call(rbind, lapply(all_element_types, get_empirical_p, observed_element_stats, sampled_element_types))
sampled_pvals_stats = melt(sampled_pvals_stats)
sampled_pvals_stats$fdr = p.adjust(sampled_pvals_stats$value, method = "fdr")
sampled_pvals_stats$fdr_stars = stars.pval(sampled_pvals_stats$fdr)

# print pvals as stars, order by median count of elements
element_type_to_fdr_stars = c(by(sampled_pvals_stats$fdr_stars, sampled_pvals_stats$element_type, paste, collapse = "/"))
combined_element_stats$fdr_label = element_type_to_fdr_stars[combined_element_stats$element_type]
combined_element_stats[combined_element_stats$value_type == "expected", "fdr_label"] = NA
element_order = combined_element_stats[combined_element_stats$value_type == "observed", ]
element_order = rev(element_order[order(element_order$median_val), "element_type"])
combined_element_stats$element_type = factor(combined_element_stats$element_type, levels = element_order)
save(combined_element_stats, file = pff("combined_element_stats.rsav"))

plt_title = paste0("non-coding element types, permutation test, n_iter=", total_times_sampled)
plt = ggplot(combined_element_stats, aes(element_type, median_val, ymin = bottom_val, ymax = top_val, fill = value_type, label = fdr_label)) + 
		geom_bar(stat = "identity", position = position_dodge()) +
		geom_text(size = 3) + 
		scale_fill_manual(values = c("expected" = "darkgrey", "observed" = "steelblue4")) +
		scale_y_continuous ("Number of genomic elements") + 
		geom_errorbar(position = position_dodge(width = 0.9), width = 0.001) + 
		plot_theme() + 
		ggtitle(plt_title)

fname = pff("figures/noncoding_elements_altered_barplot_obs_exp.pdf")
ggsave(plt, file = fname)
file_open_call2(fname)

