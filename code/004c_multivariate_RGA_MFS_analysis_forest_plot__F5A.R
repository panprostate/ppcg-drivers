source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(gplots)
library(survival)
library(survminer)
library(gtools)

PVAL_PRESELECT_FROM_UNIVAR = 0.1
PVAL_FILTER_IN_MULTIVAR = 0.05
NULL_ANALYSIS_TAG = "PPCG__No_covar"

load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
N_MUT_PATIENTS_FOR_CLIN_ANALYSIS = MIN_DRIVER_FREQ * length(all_patients)

predictor_type = "mets_time"
colname_time = "mets_time_10year"
colname_event = "mets_event_10year"
event_value = "mets"
	
# load univariate clinical stats from 004b to filter the multivariate scenario
load(file = pff("patient_clinical_data.rsav"))
load(file = pff(c("driver_clinical_stats__", NULL_ANALYSIS_TAG, ".rsav")))
load(file = pff("patient_sets_for_drivers.rsav"))
all_patients = patient_clinical_data$patient

# create driver matrix, make sure all drivers are included
drivers_stacked = data.frame(as.matrix(stack(patient_sets_for_drivers)), stringsAsFactors = FALSE)
colnames(drivers_stacked) = c("patient", "driver")
drivers_stacked$patient = factor(drivers_stacked$patient, levels = all_patients)

driver_matrix = dcast(patient ~ driver, data = drivers_stacked, fun.aggregate = length, value.var = "driver", drop = FALSE) 

univar_predictors = driver_clinical_stats[driver_clinical_stats$clin_var == predictor_type,]
univar_predictors = univar_predictors[univar_predictors$pval < PVAL_PRESELECT_FROM_UNIVAR, ]
univar_predictors = univar_predictors[univar_predictors$n_mut >= N_MUT_PATIENTS_FOR_CLIN_ANALYSIS, ]
univar_predictors = unique(univar_predictors$id)

# take drivers that survived univariate stats
driver_matrix_prefilter = driver_matrix[,colnames(driver_matrix) %in% c("patient", univar_predictors)] 

clinical_matrix = patient_clinical_data[, c("patient", colname_time, colname_event, MFS_COVARIATES), drop = FALSE]
select_complete_cases = complete.cases(clinical_matrix)
clinical_matrix = clinical_matrix[select_complete_cases,, drop = FALSE]
drivers_with_clin = merge(driver_matrix_prefilter, clinical_matrix, by = "patient")

time_to_response = drivers_with_clin[, c(colname_time, colname_event)]
this_surv = Surv(time_to_response[, colname_time], time_to_response[, colname_event] == event_value)

# two predictor matrices: one with drivers+clin, and one control as clin-only
driver_colnames = setdiff(colnames(driver_matrix_prefilter), "patient")
predictors_clin_and_drivers = drivers_with_clin[, c(MFS_COVARIATES, driver_colnames)]
predictors_clin = drivers_with_clin[, MFS_COVARIATES]

fit_main = coxph(this_surv ~ ., data = predictors_clin_and_drivers)
fit_ctrl = coxph(this_surv ~ ., data = predictors_clin)	
p_anova =  anova(fit_main, fit_ctrl)[2, "Pr(>|Chi|)"]

# data for forest plot as coefficients, Pvals, errors
forest_input_data = data.frame(coef(summary(fit_main)), stringsAsFactors = FALSE)
forest_input_data$driver = rownames(forest_input_data)
colnames(forest_input_data) = c("coef", "HR", "se_coef", "z", "pval", "driver")

#ordering and confidence intervals for HR
forest_input_data$HR_2e_lower = exp(forest_input_data$coef - 1.96 * (forest_input_data$se_coef))
forest_input_data$HR_2e_upper = exp(forest_input_data$coef + 1.96 * (forest_input_data$se_coef))
driver_order = forest_input_data[order(forest_input_data$HR), "driver"]
forest_input_data$driver = factor(forest_input_data$driver, levels = driver_order)
forest_input_data$label = signif(forest_input_data$HR, 3)
forest_input_data$predictor_type = predictor_type
forest_input_data$analysis_tag_univar = NULL_ANALYSIS_TAG
forest_input_data$multivar_covariates = paste(MFS_COVARIATES, collapse = ",")

# visualise only with Pval cutoff
forest_input_data_for_vis = forest_input_data[forest_input_data$pval < PVAL_FILTER_IN_MULTIVAR,]

multivariate_driver_clinical_stats = forest_input_data
save(multivariate_driver_clinical_stats, file = pff("multivariate_driver_clinical_stats.rsav"))
	
plot_title = paste0(predictor_type, 
		"\nanalysis_tag_univar = ", NULL_ANALYSIS_TAG,
		"\nmultivar_covariates = ", paste(MFS_COVARIATES, collapse = ","), 
		"\np_preselect < ", PVAL_PRESELECT_FROM_UNIVAR, 
		", Wald p < ", PVAL_FILTER_IN_MULTIVAR, "\nP_anova = ", signif(p_anova, 2))

plt = ggplot(forest_input_data_for_vis, 
	aes(driver, HR, ymin = HR_2e_lower, ymax = HR_2e_upper, label = label)) + 
		geom_point() +
		geom_hline(yintercept = 0, linetype = "dotted") +
		geom_linerange() + 
		geom_text (vjust = -1, size = 2) + 
		geom_text (data = forest_input_data_for_vis, y = 0, aes(driver, label = stars.pval(pval)), color = "darkred") + 
		coord_flip() + 
		plot_theme() + 
		theme() + 
		theme(axis.text.x = element_text(angle = 0)) +
		ggtitle(NULL, plot_title)

fname = pff(c("figures/Forest_plot_NEW__", predictor_type, "__", NULL_ANALYSIS_TAG, ".pdf"))
ggsave(plt, file = fname, height = 7, width = 7)
file_open_call2(fname)	

