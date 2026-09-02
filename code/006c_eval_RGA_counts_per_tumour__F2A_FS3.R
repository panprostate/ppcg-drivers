source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(ggrepel)
library(scales)
library(patchwork)
library(RColorBrewer)

load(pff("patient_sets_for_drivers.rsav"))
load(file = pff("all_patients.rsav"))

drivers_stack = data.frame(stack(patient_sets_for_drivers), stringsAsFactors = FALSE)
drivers_stack = cbind(drivers_stack[,1], do.call(rbind, strsplit(as.character(drivers_stack[,2]), split = "__")))
drivers_stack = data.frame(drivers_stack, stringsAsFactors = FALSE)
colnames(drivers_stack) = c("patient", "gene", "mut_type")

# how many drivers per gene
compute_percentage = function(mut_tag, all_mut_types, select_muts, drivers_stack, all_patients) {

	select_mut_types = all_mut_types[[mut_tag]]
	drivers_stack_here = drivers_stack[drivers_stack$mut_type %in% select_mut_types, ]
	# stack only contains observed muts; zeroes not here
	obs_count_per_patient = split(drivers_stack_here$gene, drivers_stack_here$patient)
	obs_count_per_patient = sapply(obs_count_per_patient, function(x) length(unique(x)))

	count_per_patient = structure(names = all_patients, rep(0, length(all_patients)))
	count_per_patient[names(obs_count_per_patient)] = obs_count_per_patient
	count_mean = mean(count_per_patient)
	count_median = median(count_per_patient)
	count_lower = quantile(count_per_patient, 0.025)
	count_upper = quantile(count_per_patient, 0.975)
	n_some_drivers = sum(count_per_patient > 0)
	pct_some_drivers = 100 * n_some_drivers / length(all_patients)
	
	dfr_stats = data.frame(mut_tag, count_mean, count_lower, count_upper, count_median, 
			pct_some_drivers, select_mut_types = paste(select_mut_types, collapse = ","), 
			stringsAsFactors = FALSE)
	dfr_counts = data.frame(mut_tag, patient = names(count_per_patient), count_RGAs = count_per_patient, stringsAsFactors = FALSE)
	
	list(dfr_stats, dfr_counts)
}

all_mut_types = unique(drivers_stack$mut_type)
all_mut_types = setdiff(all_mut_types, c("gain", "loss"))  # these are combined as the CNA category
all_mut_types = lapply(all_mut_types, function(x) x)
names(all_mut_types) = all_mut_types
all_mut_types$CNA = c("gain", "loss", "hAMP")
all_mut_types$all_minus_CNA = setdiff(unique(drivers_stack$mut_type), c("gain", "loss", "hAMP"))
all_mut_types$ALL = unique(drivers_stack$mut_type)

stats_counts_by_mut_type = lapply(names(all_mut_types), function(x) 
		compute_percentage(x, all_mut_types, select_muts, drivers_stack, all_patients))
		
stats_by_mut_type = do.call(rbind, lapply(stats_counts_by_mut_type, '[[', 1))
counts_by_mut_type = do.call(rbind, lapply(stats_counts_by_mut_type, '[[', 2))	
stats_by_mut_type$tag = factor(stats_by_mut_type$mut_tag, levels = stats_by_mut_type$mut_tag[order(stats_by_mut_type$count_mean)])


plt_pct_some_drivers = ggplot(stats_by_mut_type, aes(tag, pct_some_drivers, label = signif(pct_some_drivers, 2))) + 
		geom_bar(stat = "identity") + 
		geom_text(color = "darkred", vjust = -0.2, size = 2) + 
		coord_flip() + 
		scale_x_discrete(NULL) +
		scale_y_continuous("Percent") +
		plot_theme() + 
		ggtitle(NULL, "Samples with drivers")
		
plt_counts = ggplot(stats_by_mut_type, aes(tag, count_mean, ymin = count_lower, ymax = count_upper, label = signif(count_mean, 2))) + 
		geom_point(size = 3) + 
		geom_text(color = "darkred", vjust = -0.2, size = 2) + 
		geom_errorbar(width = 0.2) +  
		coord_flip() + 
		scale_x_discrete(NULL) +
		scale_y_continuous("Mean (+/- 95% CI)") +
		plot_theme() + 
		ggtitle(NULL, "Drivers per sample")

plt_combined = plt_counts / plt_pct_some_drivers

fname = pff("figures/plot_driver_counts_per_sample.pdf")
ggsave(plt_combined, file = fname, width = 4, height = 4)
file_open_call2(fname)

stats_drivers_per_sample = stats_by_mut_type
counts_drivers_per_sample = counts_by_mut_type
save(counts_drivers_per_sample, file = pff("counts_drivers_per_sample.rsav"))
save(stats_drivers_per_sample, file = pff("stats_drivers_per_sample.rsav"))


# histogram drivers per sample
counts_drivers_per_sample$mut_tag = factor(counts_drivers_per_sample$mut_tag, levels = rev(stats_by_mut_type$tag))

plt = ggplot(counts_drivers_per_sample, aes(x = count_RGAs)) + 
		geom_histogram(bins = 8) + 
		facet_wrap(~mut_tag, scale = "free", ncol = 4) +
		scale_x_continuous("Number of tumours") +
		scale_y_continuous("Number of RGAs per tumour") + 
		plot_theme() + 
		ggtitle(NULL, "RGAs per tumour")

fname = pff("figures/histogram_driver_counts_per_sample_v2.pdf")
ggsave(plt, file = fname, width = 8, height = 6)
file_open_call2(fname)

