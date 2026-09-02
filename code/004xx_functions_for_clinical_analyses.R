library(survival)
library(survminer)
library(patchwork)


time_to_event = function(clin1, col_interval, col_event, event_value, min_mut_patients, covariates, clin_var_label) {
	
	if (length(covariates) == 0) {
		h0_formula = "this_surv ~ 1"
		h1_formula = "this_surv ~ is_mutated"
	} else {
		h0_formula = paste("this_surv ~", paste(covariates, collapse = " + "))
		h1_formula = paste(h0_formula, " + is_mutated")
	}
	
	# select the variables we need 
	clin_here = clin1[, c(covariates, col_interval, col_event, "is_mutated", "patient")]
	# try converting to numerics
	clin_here$this_interval = as.numeric(clin_here[, col_interval])
	# convert interval to years for visualisation purposes
	clin_here$this_interval = clin_here$this_interval / 365
	
	clin_here$this_event = clin_here[, col_event] == event_value

	# remove any NAs using complete.cases
	clin_here = clin_here[complete.cases(clin_here),]
	
	patients_mutated = unique(clin_here[clin_here$is_mutated == "mutated", "patient"])
	patients_mutated_concat = paste(patients_mutated, collapse = ",")

	# skip this analysis if after interval=NA filtering few mutated patients are found
	n_mut = sum(clin_here$is_mutated == "mutated")
	if (n_mut < min_mut_patients) {
		plt = ggplot(clin_here) + geom_blank()	
		pval = pval_logrank = HR = effect_lo = effect_hi = 1 
		dfr = data.frame(clin_var = clin_var_label, pval, effect = HR, effect_lo, effect_hi, n_mut, patients_mutated_concat, pval_logrank,
				stringsAsFactors = FALSE)
		return(list(plt, dfr))
	}
	
	clin_here$is_mutated = factor(clin_here$is_mutated, levels = c("not_mutated", "mutated"))

	this_surv = Surv(clin_here$this_interval, clin_here$this_event)
	h0 = coxph(as.formula(h0_formula), data = clin_here)
	h1 = coxph(as.formula(h1_formula), data = clin_here)
	pval = anova(h0, h1, test = "Chisq")[2,4]

	model_summary = summary(h1)
	ci_matrix = model_summary$conf.int
	HR = ci_matrix["is_mutatedmutated", "exp(coef)"]
	HR_lower = ci_matrix["is_mutatedmutated", "lower .95"]
	HR_upper = ci_matrix["is_mutatedmutated", "upper .95"]
	
	this_formula = paste( 
			'Surv(clin_here$this_interval, clin_here$this_event)', 
			'~ is_mutated')
	ggsurv_fit <- survfit(as.formula(this_formula), data = clin_here)
	ggsurv_fit$call$formula = as.formula(this_formula)
	pval_logrank = survdiff(ggsurv_fit$call$formula, data = clin_here)$pval

	plt1 = ggsurvplot(ggsurv_fit, data = clin_here, risk.table = FALSE, pval = signif(pval_logrank, 2),
			legend.labs = c("not mutated", "mutated"), palette = c("grey", "red"))
	plt1 = plt1 + ggtitle(clin_var_label, 
			paste0("p= ", signif(pval, 2), "; n= ", n_mut, "\nHR= ", 
					signif(HR, 2), " [", signif(HR_lower, 2), "-", signif(HR_upper, 2), "]"))
	plt1 = plt1$plot
	
	dfr = data.frame(clin_var = clin_var_label, pval, effect = HR, 
			effect_lo = HR_lower, effect_hi = HR_upper, n_mut, patients_mutated_concat, pval_logrank, 
			stringsAsFactors = FALSE)
	list(plt1, dfr)
}

test_binary_category = function(clin1, colname, col_value, min_mut_patients, covariates, clin_var_label) {

	# avoid circularity	
	if (colname == "early_onset") {
		covariates = setdiff(covariates, "early_onset")
	}
	
	if (length(covariates) == 0) {
		h0_formula = "outcome ~ 1"
		h1_formula = "outcome ~ is_mutated"
	} else {
		h0_formula = paste("outcome ~", paste(covariates, collapse = " + "))
		h1_formula = paste(h0_formula, " + is_mutated")
	}
		
	# select the variables we need 
	clin_here = clin1[, c(covariates, colname, "is_mutated", "patient"), drop = FALSE]
	clin_here$outcome = 0+(clin_here[, colname] == col_value)

	# remove any NAs using complete.cases
	clin_here = clin_here[complete.cases(clin_here),]
	patients_mutated = unique(clin_here[clin_here$is_mutated == "mutated", "patient"])
	patients_mutated_concat = paste(patients_mutated, collapse = ",")

	clin_here$is_mutated = factor(clin_here$is_mutated, levels = c("not_mutated", "mutated"))
	
	# skip this analysis if after interval=NA filtering few mutated patients are found
	n_mut = sum(clin_here$is_mutated == "mutated")
	if (n_mut < min_mut_patients) {
		plt = ggplot(clin_here) + geom_blank()	
		pval = pval_logrank = or = effect_lo = effect_hi = 1
		dfr = data.frame(clin_var = clin_var_label, pval, effect = or, effect_lo, effect_hi, n_mut, patients_mutated_concat, pval_logrank, 
				stringsAsFactors = FALSE)
		return(list(plt, dfr))
	}
	
	h0 = glm(as.formula(h0_formula), data = clin_here, family = binomial)
	h1 = glm(as.formula(h1_formula), data = clin_here, family = binomial)
	pval = anova(h0, h1, test = "Chisq")[2, "Pr(>Chi)"]
	or = exp(coef(h1)[['is_mutatedmutated']])
	
	# Calculate the profile likelihood confidence intervals
	ci_matrix = confint(h1)
	OR_lower = exp(ci_matrix["is_mutatedmutated", "2.5 %"])
	OR_upper = exp(ci_matrix["is_mutatedmutated", "97.5 %"])
	
	clin_here$binary_class = clin_here[,colname] == col_value
	clin_here$is_mutated = factor(clin_here$is_mutated, levels = c("mutated", "not_mutated"))
	colors_here = c("TRUE" = "darkcyan", "FALSE" = "lightgrey")
	
	plt = ggplot(clin_here, aes(factor(1), fill = binary_class)) + 
			geom_bar(color = "black") + 
			scale_x_discrete(NULL, breaks = NULL, labels = NULL) + 
			scale_fill_manual(NULL, values = colors_here) +
			facet_wrap(~is_mutated, scale = "free") + 
			plot_theme() +
			ggtitle(clin_var_label, paste0("p = ", signif(pval, 2), "; n = ", n_mut, 
					"\nOR = ", signif(or, 2), " [", signif(OR_lower, 2), "-", signif(OR_upper, 2), "]")) + 
			theme(legend.position = "bottom", plot.title = element_text(size = 13), , plot.subtitle = element_text(size = 11))

	dfr = data.frame(clin_var = clin_var_label, pval, effect = or, 
			effect_lo = OR_lower, effect_hi = OR_upper, n_mut, patients_mutated_concat, pval_logrank = NA,
			stringsAsFactors = FALSE)
	list(plt, dfr)
}


clinical_analysis = function(i, patient_clinical_sets, clin, do_plot, min_mut_patients, covariates, precomp_clin_stats) {
	
	mut_patients = unique(patient_clinical_sets[[i]])
	patient_set_id = names(patient_clinical_sets)[[i]]
	
	cat(".")
	
	if (length(mut_patients) < min_mut_patients)  {
		return(NULL)
	}
	
	clin$is_mutated = c("not_mutated", "mutated")[1 + clin$patient %in% mut_patients]
	n_mutated_patients = sum(clin$patient %in% mut_patients)
	plot_title = paste0(patient_set_id, "; n=", n_mutated_patients, ";\ncovr=", paste(covariates, collapse = ","))
	
	# if clinical stats are precomputed and added here, we can add the FDR to the plot title
	this_merged_fdr = NA
	if (!is.null(precomp_clin_stats)[[1]] & patient_set_id %in% precomp_clin_stats$id) {
		this_merged_fdr = precomp_clin_stats[precomp_clin_stats$id == patient_set_id, "merged_fdr"]
		plot_title = paste0("; ", plot_title, "; FDR=", signif(this_merged_fdr, 2))
	}

	# time to relapse KM plot; 10-year censored
	stats_time_to_relapse = time_to_event(clin, "relapse_time_10year", "relapse_event_10year", "relapsed", 
			min_mut_patients, covariates, "relapse_time")
	plt_relapse_time = stats_time_to_relapse[[1]]
	dfr_relapse_time = stats_time_to_relapse[[2]]
	
	# time to mets KM plot; 10-year
	stats_time_to_mets = time_to_event(clin, "mets_time_10year", "mets_event_10year", "mets", 
			min_mut_patients, covariates, "mets_time")
	plt_mets_time = stats_time_to_mets[[1]]
	dfr_mets_time = stats_time_to_mets[[2]]
	
#	# metastatic biology	
	colname = "new_metastatic_biology_indicator"
	col_value = "mets_biol"
	stats_metsbio = test_binary_category(clin, colname, col_value, min_mut_patients, covariates, "mets_biol_ind")
	plt_metsbio = stats_metsbio[[1]]
	dfr_metsbio = stats_metsbio[[2]]
	
#	# early onset
	colname = "early_onset"
	col_value = TRUE
	stats_EO = test_binary_category(clin, colname, col_value, min_mut_patients, covariates, "early_onset_55")
	plt_EO = stats_EO[[1]]
	dfr_EO = stats_EO[[2]]

	# combine individual plots above	
	combo_plt = (plt_relapse_time + plt_mets_time) / (plot_spacer() + plt_metsbio + plt_EO)
	combo_plt = combo_plt + plot_annotation(title = plot_title) # + plot_layout(widths = unit(c(1, 1, 0.5, 0.5), 'null'))
	if (do_plot) {
		print(combo_plt)
	}

	pval_dfr = rbind(dfr_relapse_time, dfr_mets_time, dfr_metsbio, dfr_EO)
	dfr = data.frame(id = patient_set_id, n_mut_patients = length(mut_patients), pval_dfr, stringsAsFactors = FALSE)
	rownames(dfr) = NULL
	dfr
}


censor_time_to_event = function(time_col, event_col, event_value, no_event_value, n_years_censor, new_time_col, new_event_col, patient_clinical_data) {
		
	n_days_censor = (365 * n_years_censor) + floor(1 * n_years_censor / 4)
	
	patient_clinical_data1 = patient_clinical_data
	patient_clinical_data1[, new_time_col] = patient_clinical_data1[, time_col]
	patient_clinical_data1[, new_event_col] = patient_clinical_data1[, event_col]

	# everyone's last observation of survival happens at censoring point	
	which_to_censor_time = !is.na(patient_clinical_data1[, time_col]) & patient_clinical_data1[,time_col] > n_days_censor
	which_to_censor_event = which_to_censor_time & !is.na(patient_clinical_data1[, event_col]) & patient_clinical_data1[,event_col] == event_value

	which_to_censor_time = which(which_to_censor_time)
	which_to_censor_event = which(which_to_censor_event)
	
	patient_clinical_data1[which_to_censor_time, new_time_col] = n_days_censor
	patient_clinical_data1[which_to_censor_event, new_event_col] = no_event_value	
	patient_clinical_data1
}


multigroup_time_to_event = function(clin1, col_interval, col_event, event_value, covariates, min_mut_patients, group_colors, plot_tag) {
	
	columns_needed = c(col_event, col_interval, "group_id", covariates)
	clin2 = clin1[,columns_needed]
	clin2$this_interval = as.numeric(clin2[, col_interval])
	clin2 = clin2[complete.cases(clin2), ]
	
	# convert interval to years for visualisation purposes
	clin2$this_interval = clin2$this_interval / 365
	clin2$this_event = 0 + (clin2[, col_event] == event_value)

	# skip this analysis if after interval=NA filtering few mutated patients are found
	if (any(table(clin2$group_id) < min_mut_patients)) {
		cat("not enough patients")
		plt = ggplot(clin2) + geom_blank()	
		pval = HR = 1
		return(list(plt, pval, HR))
	}

	# consistent ordering 
	clin2$group_id = factor(clin2$group_id, levels = sort(unique(clin2$group_id)))
	
	# add covariates as specificed previously
	h0_formula = paste("this_surv ~", paste(covariates, collapse = " + "))
	h1_formula = paste(h0_formula, " + group_id")
	
	this_surv = Surv(clin2$this_interval, clin2$this_event)
	h0 = coxph(as.formula(h0_formula), data = clin2)
	h1 = coxph(as.formula(h1_formula), data = clin2)
	pval = anova(h0, h1, test = "Chisq")[2,4]
	HR = NA

	this_formula = paste('Surv(clin2$this_interval, clin2$this_event)', '~ group_id')
	ggsurv_fit <- survfit(as.formula(this_formula), data = clin2)
	ggsurv_fit$call$formula = as.formula(this_formula)
	pval_logrank = survdiff(ggsurv_fit$call$formula, data = clin2)$pval
		
	plt1 = ggsurvplot(ggsurv_fit, data = clin2, risk.table = FALSE, pval = signif(pval_logrank, 2), 
			legend.labs = levels(clin2$group_id), legend = "right", 
			palette = group_colors[levels(clin2$group_id)])
			
	ggtitl = paste0(plot_tag, "\n", col_interval,", n=", nrow(clin2), "\np=",signif(pval, 2))
	plt1 = plt1 + ggtitle(label = ggtitl)
	
	dfr_stats = data.frame(plot_tag, pval, HR, pval_logrank, col_interval, stringsAsFactors = FALSE)
	
	list(plt1$plot, dfr_stats)
}



get_clinical_cluster_stats_vs_rest_continuous = function(column_to_test, column_for_group, clin1, alternative) {
	res = do.call(rbind, lapply(unique(clin1[, column_for_group]), function(grp) 
			data.frame(grp, pval = 
				wilcox.test(
					clin1[clin1[, column_for_group] == grp, column_to_test], 
					clin1[clin1[, column_for_group] != grp, column_to_test], 
					alt = alternative)$p.val)
				))
	res$pval_stars = stars.pval(res$pval)
	res$variable = column_to_test
	res
}



get_clinical_cluster_stats_vs_rest_discrete = function(column_to_test, value_to_test, column_for_group, clin1, alternative) {
	res = do.call(rbind, lapply(unique(clin1[, column_for_group]), function(grp) 
			data.frame(grp, pval = 
				fisher.test(
					clin1[, column_to_test] == value_to_test,
					clin1[, column_for_group] == grp,
					alt = alternative)$p.val
					)))
	res$pval_stars = stars.pval(res$pval)
	res$variable = column_to_test
	res
}


assign_median_val_label = function(clin1, clin_column, group_column) {
	ave_label = ave(clin1[, clin_column], clin1[, group_column], FUN = function(x) median(x, na.rm = T))	
	ave_label[duplicated(clin1[, group_column])] = NA
	ave_label
}



assign_sum_val_label = function(clin1, clin_column, clin_value, group_column) {
	ave_label = ave(clin1[, clin_column] == clin_value, clin1[, group_column], FUN = function(x) sum(x, na.rm = T))	
	ave_label[duplicated(clin1[, group_column])] = NA
	ave_label
}


get_multi_group_report = function(group_id, input_data_list, COVARIATES, MIN_MUT_PATIENTS) {
	
	cat(group_id, "\n")

	group_colors = input_data_list[[group_id]][["group_colors"]]
	clin1 = input_data_list[[group_id]][["clin1"]]
	
	clin1$T_stage_num = as.numeric(gsub("^T", "", clin1$T_stage))
			
	res_relapse = multigroup_time_to_event(clin1, "relapse_time_10year", "relapse_event_10year", "relapsed", 
			COVARIATES, MIN_MUT_PATIENTS, group_colors, group_id)
	res_mets = multigroup_time_to_event(clin1, "mets_time_10year", "mets_event_10year", "mets", 
			COVARIATES, MIN_MUT_PATIENTS, group_colors, group_id)
			
	pga_stat = get_clinical_cluster_stats_vs_rest_continuous("PGA", "group_id", clin1, "greater")
	age_stat = get_clinical_cluster_stats_vs_rest_continuous("age_at_tumour_collection", "group_id", clin1, "less")
	grade_stat = get_clinical_cluster_stats_vs_rest_continuous("Gleason_group", "group_id", clin1, "greater")
	stage_stat = get_clinical_cluster_stats_vs_rest_continuous("T_stage_num", "group_id", clin1, "greater")
	PSA_stat = get_clinical_cluster_stats_vs_rest_continuous("PSA_log2", "group_id", clin1, "greater")
	TMB_stat = get_clinical_cluster_stats_vs_rest_continuous("GTMB_log1p", "group_id", clin1, "greater")
	RGA_stat = get_clinical_cluster_stats_vs_rest_continuous("n_rgas", "group_id", clin1, "greater")
	
	clin1$cl_pga_median = assign_median_val_label(clin1, "PGA", "group_id")
	clin1$cl_age_median = assign_median_val_label(clin1, "age_at_tumour_collection", "group_id")
	clin1$cl_Gleason_median = assign_median_val_label(clin1, "Gleason_group", "group_id")
	clin1$cl_stage_median = assign_median_val_label(clin1, "T_stage_num", "group_id")
	clin1$cl_PSA_median = assign_median_val_label(clin1, "PSA_log2", "group_id")
	clin1$cl_TMB_median = assign_median_val_label(clin1, "GTMB_log1p", "group_id")
	clin1$cl_RGA_median = assign_median_val_label(clin1, "n_rgas", "group_id")

	ancestry_stat = get_clinical_cluster_stats_vs_rest_discrete("predicted_ancestry", "AFR", "group_id", clin1, "greater")

	plt_pga = ggplot(clin1[!is.na(clin1$PGA),], aes(group_id, PGA, fill = group_id)) + 
			geom_boxplot(notch = FALSE) + 
			geom_text(size = 3, y = 0, color = "darkgrey", aes(group_id, label = signif(cl_pga_median, 2))) + 
			geom_text(data = pga_stat, aes(grp, label = pval_stars, y = max(clin1$PGA, na.rm = TRUE))) + 
			scale_fill_manual("grp", values = group_colors) +
			scale_x_discrete(NULL) + scale_y_continuous("PGA") +
			plot_theme() +
			theme(legend.position = "bottom")

	plt_age = ggplot(clin1[!is.na(clin1$age_at_tumour_collection),], aes(group_id, age_at_tumour_collection, fill = group_id)) + 
			geom_boxplot(notch = FALSE) + 
			geom_text(size = 3, y = 35, color = "darkgrey", aes(group_id, label = signif(cl_age_median, 2))) + 
			geom_text(data = age_stat, aes(grp, label = pval_stars, y = max(clin1$age_at_tumour_collection, na.rm = TRUE))) + 
			scale_fill_manual("grp", values = group_colors) +
			scale_x_discrete(NULL) + scale_y_continuous("Age") +
			plot_theme() +
			theme(legend.position = "bottom")
			
	plt_grade = ggplot(clin1[!is.na(clin1$Gleason_group),]) + 
			geom_bar(position = "fill", aes(group_id, fill = factor(Gleason_group))) + 
			geom_text(size = 3, y = 0, color = "darkgrey", aes(x = group_id, label = signif(cl_Gleason_median, 2))) + 
			geom_text(data = grade_stat, aes(grp, label = pval_stars, y = 1)) + 
			scale_fill_brewer("grade", palette = "OrRd") + 
			scale_x_discrete(NULL) + scale_y_continuous("Grade (GGG)") +
			plot_theme() +
			theme(legend.position = "bottom")

	plt_stage = ggplot(clin1[!is.na(clin1$T_stage),]) + 
			geom_bar(position = "fill", aes(group_id, fill = factor(T_stage))) + 
			geom_text(size = 3, y = 0, color = "darkgrey", aes(group_id, label = signif(cl_stage_median, 2))) +
			geom_text(data = stage_stat, aes(grp, label = pval_stars, y = 1)) + 
			scale_fill_brewer("stage", palette = "PuRd") + 
			scale_x_discrete(NULL) + scale_y_continuous("Stage") +
			plot_theme() +
			theme(legend.position = "bottom")

	plt_PSA = ggplot(clin1[!is.na(clin1$PSA_log2),], aes(group_id, PSA_log2, fill = group_id)) + 
			geom_boxplot(notch = FALSE) + 
			geom_text(size = 3, y = 0, color = "darkgrey", aes(group_id, label = signif(cl_PSA_median, 2))) + 
			scale_x_discrete(NULL) + scale_y_continuous("PSA,log2") +
			scale_fill_manual("grp", values = group_colors) +
			geom_text(data = PSA_stat, aes(grp, label = pval_stars, y = max(clin1$PSA_log2, na.rm = TRUE))) + 
			plot_theme() +
			theme(legend.position = "bottom")
			
	plt_TMB = ggplot(clin1[!is.na(clin1$GTMB_log1p),], aes(group_id, GTMB_log1p, fill = group_id)) + 
			geom_boxplot(notch = FALSE) + 
			geom_text(size = 3, y = 0, color = "darkgrey", aes(group_id, label = signif(cl_TMB_median, 2))) + 
			scale_x_discrete(NULL) + scale_y_continuous("gTMB, log1p") +
			scale_fill_manual("grp", values = group_colors) +
			geom_text(data = TMB_stat, aes(grp, label = pval_stars, y = max(clin1$GTMB_log1p, na.rm = TRUE))) + 
			plot_theme() +
			theme(legend.position = "bottom")
	
	plt_RGA = ggplot(clin1[!is.na(clin1$n_rgas),], aes(group_id, n_rgas, fill = group_id)) + 
			geom_boxplot(notch = FALSE) + 
			geom_text(size = 3, y = 0, color = "darkgrey", aes(group_id, label = signif(cl_RGA_median, 2))) + 
			scale_x_discrete(NULL) + scale_y_continuous("n_RGAs") +
			scale_fill_manual("grp", values = group_colors) +
			geom_text(data = RGA_stat, aes(grp, label = pval_stars, y = max(clin1$n_rgas, na.rm = TRUE))) + 
			plot_theme() +
			theme(legend.position = "bottom")

	clin1$is_AFR = clin1$predicted_ancestry == "AFR"
	pct_cluster_AFR = stack(by(clin1$is_AFR, clin1$group_id, function(x) 100 * sum(x, na.rm = T) / sum(!is.na(x))))
	count_cluster_AFR = stack(by(clin1$is_AFR, clin1$group_id, function(x) sum(x, na.rm = T) ))
	stat_cluster_afr = merge(count_cluster_AFR, pct_cluster_AFR, by = "ind")
	colnames(stat_cluster_afr) = c("group_id", "n_samples", "pct_samples")
			
	plt_ancestry = ggplot(stat_cluster_afr, aes(group_id, pct_samples, label = n_samples)) + 
			geom_bar(stat = "identity") +
			geom_text(size = 3, y = 0, color = "darkgrey", aes(group_id, label = n_samples)) +
			geom_text(data = ancestry_stat, aes(grp, label = pval_stars, y = max(stat_cluster_afr$pct_samples))) + 
			scale_x_discrete(NULL) + scale_y_continuous("% cluster AFR") +
			plot_theme() +
			theme(legend.position = "bottom")
			
	plt_n_patient = ggplot(clin1, aes(group_id, after_stat(count), fill = group_id)) + 
			geom_bar() + 
			geom_text(stat = "count", aes(y = after_stat(count)-0.5, label=after_stat(count))) + 
			scale_fill_manual("grp", values = group_colors) +
			plot_theme() +
			scale_x_discrete(NULL) + scale_y_continuous("n_pat") +
			theme(legend.position = "bottom")
	
	combo_plt = 
			( res_mets[[1]] | res_relapse[[1]] ) / 
			( plt_n_patient | plt_grade | plt_stage | plt_age | plt_PSA) /
			( plot_spacer() | plt_ancestry  | plt_pga | plt_TMB | plt_RGA ) +
			plot_annotation(title = paste0(group_id, "\n", paste(COVARIATES, collapse = ", ")))

	res = list()	
	res[[1]] = combo_plt
	res[[2]] = rbind(res_mets[[2]], res_relapse[[2]])
	res[[3]] = rbind(pga_stat, age_stat, grade_stat, stage_stat, PSA_stat, TMB_stat, RGA_stat, ancestry_stat)
	res	
}


