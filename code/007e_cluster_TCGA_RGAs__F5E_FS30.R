source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(survminer)
library(survival)
library(caret)
library(randomForest)
library(ranger)
library(patchwork)

clust_to_try = SELECT_CLUSTER_METHOD
n_clusters = SELECT_N_CLUSTERS

# perform this on most frequent RGAs only
load(pff("all_patients.rsav"))
MIN_DRIVER_FREQ = 0.01
MIN_N_PATIENTS = MIN_DRIVER_FREQ * length(all_patients)

load(file = pff("TCGA_validation/TCGA_clinical_data.rsav"))
load(file = pff("TCGA_validation/PPCG_driver_matrix.rsav"))
load(file = pff("TCGA_validation/TCGA_driver_matrix.rsav"))
load(pff("patient_cluster_membership.rsav"))
PPCG_patient_cluster = structure(names = patient_cluster_membership$patient, patient_cluster_membership$group_id)

select_RGAs = names(which(apply(PPCG_driver_matrix, 2, sum) >= MIN_N_PATIENTS))
TCGA_driver_matrix = TCGA_driver_matrix[,select_RGAs]
PPCG_driver_matrix = PPCG_driver_matrix[,select_RGAs]
PPCG_driver_dfr = data.frame(PPCG_driver_matrix[rownames(PPCG_driver_matrix) %in% patient_cluster_membership$patient,])
PPCG_driver_dfr$cluster_label = factor(PPCG_patient_cluster[rownames(PPCG_driver_dfr)])
PPCG_driver_dfr$cluster_label = relevel(PPCG_driver_dfr$cluster_label, ref = "cluster_1")
TCGA_driver_dfr = data.frame(TCGA_driver_matrix, stringsAsFactors = FALSE)

# remove a point formerly a dash from variable names
colnames(TCGA_driver_dfr) = gsub("\\.", "_", colnames(TCGA_driver_dfr))
colnames(PPCG_driver_dfr) = gsub("\\.", "_", colnames(PPCG_driver_dfr))


# random forest: train and tune using balanced sampling; split into train (70%) and test (30%)
set.seed(123)
train_idx = createDataPartition(PPCG_driver_dfr$cluster_label, p = 0.7, list = FALSE)
train_data = PPCG_driver_dfr[train_idx, ]
test_data = PPCG_driver_dfr[-train_idx, ]

# cross-validation control with down-sampling for class imbalance
cv_control = trainControl(
	method = "cv",
	number = 5,
	sampling = "down",
	classProbs = TRUE,
	verboseIter = TRUE
)

tune_grid = expand.grid(
	mtry = c(5, 9, 15, 25),
	splitrule = c("gini", "extratrees"),
	min.node.size = c(1, 3, 5)
)

rf_tuned = train(
	cluster_label ~ .,
	data = train_data,
	method = "ranger",
	trControl = cv_control,
	tuneGrid = tune_grid,
	metric = "Kappa",             
	num.trees = 750,
	importance = "permutation"     
)

fname = pff(c("figures/random_forest_for_PPCG_cluster_label_prediction_ccp__", clust_to_try, "__K", n_clusters, ".pdf"))
pdf(fname)
plot(rf_tuned)
dev.off()
file_open_call2(fname)

# predict on test data using the best model
PPCG_test_predictions = predict(rf_tuned, newdata = test_data)

# evaluate performance
PPCG_test_evaluation = confusionMatrix(PPCG_test_predictions, test_data$cluster_label)
conf_matrix = melt(PPCG_test_evaluation$table)
conf_matrix_proportional <- as.data.frame(prop.table(PPCG_test_evaluation$table, margin = 1))
conf_matrix_both = merge(conf_matrix_proportional, conf_matrix, by = c("Prediction", "Reference"))
performance_stats = paste(paste(names(PPCG_test_evaluation$overall), signif(PPCG_test_evaluation$overall, 2)), collapse = "\n")

plt1 = ggplot(data = conf_matrix_both, aes(x = Prediction, y = Reference, fill = Freq)) +
	  geom_tile(color = "white") +
	  geom_text(aes(label = value), color = "black") +
	  scale_fill_gradient(low = "#e0f2fe", high = "#0284c7") + 
	  theme_minimal() +
	  labs(title = "Confusion Matrix", x = "Predicted", y = "Actual") +
	  ggtitle(NULL, paste0(clust_to_try, "; K=", n_clusters, "\n", performance_stats))
	  
stats_dfr = melt(PPCG_test_evaluation$byClass)
stats_dfr$Var1 = gsub("Class: cluster_(.)", "c\\1", stats_dfr$Var1)

plt2 = ggplot(stats_dfr, aes(Var1, value)) + 
		geom_bar(stat = "identity") + 
		facet_wrap(~Var2, scale = "free") + 
		plot_theme()

plt = plt1 / plt2

fname = pff(c("figures/random_forest_confusion_matrix_for_PPCG_cluster_label_prediction_ccp__", clust_to_try, "__K", n_clusters, ".pdf"))
ggsave(plt, file = fname, height = 10)
file_open_call2(fname)

save(rf_tuned, file = pff("TCGA_validation/rf_tuned.rsav"))
save(PPCG_test_evaluation, file = pff("TCGA_validation/PPCG_test_evaluation.rsav"))

TCGA_predicted_classes_RF = predict(rf_tuned, newdata = TCGA_driver_dfr, type = "raw")
TCGA_sample_cluster = structure(as.character(TCGA_predicted_classes_RF), names = rownames(TCGA_driver_dfr))
dfr_TCGA_cluster_membership = data.frame(TCGA_donor_id = names(TCGA_sample_cluster), oncotype = TCGA_sample_cluster, stringsAsFactors = FALSE)
save(dfr_TCGA_cluster_membership, file = pff("TCGA_validation/dfr_TCGA_cluster_membership.rsav"))


model_setup = "this_surv ~ early_onset"
TCGA_clinical_data$cluster = TCGA_sample_cluster[TCGA_clinical_data$patient]
TCGA_clinical_data$cluster = factor(TCGA_clinical_data$cluster)
TCGA_clinical_data$cluster = relevel(TCGA_clinical_data$cluster, ref = "cluster_1")

this_surv = Surv(TCGA_clinical_data$time_PFI_10year, TCGA_clinical_data$status_PFI_10year)
h0 = coxph(as.formula(model_setup), data = TCGA_clinical_data)
h1 = coxph(as.formula(paste0(model_setup, " + cluster")), data = TCGA_clinical_data)
pval = anova(h0, h1, test = "Chisq")[2,4]

this_formula = paste('Surv(TCGA_clinical_data$time_PFI_10year, TCGA_clinical_data$status_PFI_10year)', '~ cluster')
ggsurv_fit <- survfit(as.formula(this_formula), data = TCGA_clinical_data)
ggsurv_fit$call$formula = as.formula(this_formula)

pval_logrank = survdiff(ggsurv_fit$call$formula, data = TCGA_clinical_data)$pval

plot_title = paste0("TCGA\nP=", signif(pval, 2), "; n=", nrow(TCGA_clinical_data), "; cert=0; ", model_setup)
plot_title = paste0(plot_title, "\n", 
		paste(names(table(TCGA_clinical_data$cluster)), 
		table(TCGA_clinical_data$cluster), collapse = ", ", sep = "="))
	
plt1 = ggsurvplot(ggsurv_fit, data = TCGA_clinical_data, risk.table = TRUE, pval = TRUE, 
		legend.labs = levels(TCGA_clinical_data$cluster), legend = "right", 
		palette = CLUSTER_COLORS[as.numeric(gsub("cluster_", "", levels(TCGA_clinical_data$cluster)))], 
		title = plot_title)
combo_plt = ( plt1[[1]] | plt1[[1]] ) / 
			( plot_spacer()) /
			( plot_spacer()) /
			( plot_spacer() ) +
			plot_annotation(title = "TCGA validation\nKM plot")

fname = pff(c("figures/TCGA_clusters_km_plot__CCP_k4_cert0__", clust_to_try, "__K", n_clusters, ".pdf"))
ggsave(combo_plt, file = fname, width = 12, height = 12)
file_open_call2(fname)


dfr_TCGA_validation_coxph = data.frame(pval_cox = pval, t(exp(coef(h1))), pval_logrank, tag = "TCGA", stringsAsFactors = FALSE)
save(dfr_TCGA_validation_coxph, file = pff("TCGA_validation/dfr_TCGA_validation_coxph.rsav"))


TCGA_cluster_freq = 
		data.frame(fraction_cohort = table(TCGA_sample_cluster) / length(TCGA_sample_cluster), 
				cluster_id = names(table(TCGA_sample_cluster)))
PPCG_cluster_freq = 
		data.frame(fraction_cohort = table(PPCG_driver_dfr$cluster_label) / nrow(PPCG_driver_dfr), 
				cluster_id = names(table(PPCG_driver_dfr$cluster_label)))

dfr1 = data.frame(cohort = "TCGA", cluster = TCGA_sample_cluster, patient = names(TCGA_sample_cluster), stringsAsFactors = FALSE)
dfr2 = data.frame(cohort = "PPCG", cluster = as.character(PPCG_driver_dfr$cluster_label), patient = rownames(PPCG_driver_dfr), stringsAsFactors = FALSE)
dfr_comb = rbind(dfr1, dfr2)

plt = ggplot(dfr_comb, aes(factor(1), fill = cluster, group = cluster)) + 
		geom_bar() + 
		plot_theme() + 
		facet_wrap(~cohort, scale = "free") + 
		geom_text(
			stat = "count", 
			aes(label = after_stat(count)), 
			vjust = -0.5,
			position = position_stack(vjust = 0.5)
		) +
		scale_fill_manual(values = CLUSTER_COLORS[1:length(unique(dfr_comb$cluster))])
		
fname = pff("figures/TCGA_PPCG_cluster_freqs.pdf")
ggsave(plt, file = fname, width = 5)
file_open_call2(fname)

