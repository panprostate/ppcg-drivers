source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(eulerr)

load(pff("results_signf_merged_annot.rsav"))
drivers_by_mut_type = split(results_signf_merged_annot$annots_MAIN, results_signf_merged_annot$mut_type)
all_drivers = setdiff(unique(unlist(drivers_by_mut_type)), "ETS")


# euler diagram for sites, plot a few and select the best arranged visual
dfr_for_venn = data.frame(
		all_drivers, 
		SV = all_drivers %in% drivers_by_mut_type[["SV"]],
		SNV_NC = all_drivers %in% drivers_by_mut_type[["SNV_NC"]],
		SNV_CDS = all_drivers %in% drivers_by_mut_type[["SNV_CDS"]],
		stringsAsFactors = FALSE)
dfr_for_venn$mut_types = apply(dfr_for_venn, 1, function(x) paste(names(x)[-1][as.logical(x[-1])], collapse = ","))

top_drivers_by_mut_type = lapply(split(dfr_for_venn$all_drivers, dfr_for_venn$mut_type), function(x) paste(setdiff(x[1:5], NA), collapse = "\n"))
top_drivers_by_mut_type = data.frame(as.matrix(stack(top_drivers_by_mut_type)), stringsAsFactors = FALSE)
colnames(top_drivers_by_mut_type) = c("top_drivers", "mut_type")

plt_topgenes = ggplot(top_drivers_by_mut_type, aes(mut_type, factor(1), label = top_drivers))  +
		geom_text(size = 2)  +
		plot_theme()


set.seed(12345)

fname = pff("figures/VennEuler_drivers_by_mut_type_top5.pdf")
pdf(fname)

for (i in 1:10) {
	cat(i, " ")
	fit2 = euler(dfr_for_venn[,c("SV", "SNV_NC", "SNV_CDS")], shape = "ellipse")
	print(plot(fit2, quantities = TRUE))
}

print(plt_topgenes)

dev.off()
file_open_call2(fname)



