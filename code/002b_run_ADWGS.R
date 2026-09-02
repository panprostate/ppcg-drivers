source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")
library(ActiveDriverWGS)


args = commandArgs(trailingOnly = TRUE)
if (length(args)>=3) {
	iteration = as.integer(args[1])
	n_jobs_total = as.integer(args[2])
	job_id = as.integer(args[3])
	cat("\n\narg iteration=", iteration, "; n_jobs_total=", n_jobs_total, "; job_id=", job_id, "\n\n")
} else {
	stop("parameters missing")
}


load(pff(c("ADWGS_tmp/todo_list.", iteration, ".rsav")))
load(pff("gr_prepared_elements.rsav"))

this_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
todo_here = todo_list[which(1:nrow(todo_list) %% n_jobs_total == (job_id-1)),]

all_prep_vars = list()
load(file = pff(paste0("gr_prepared_variants__SNV_indel.rsav")))
all_prep_vars[["SNV_indel"]] = gr_prepared_variants
load(file = pff(paste0("gr_prepared_variants__SV.rsav")))
all_prep_vars[["SV"]] = gr_prepared_variants


get_ADWGS = function(i, todo_here, gr_prepared_elements, all_prep_vars, this_genome) {
	el = todo_here[i, "element"]
	ds = todo_here[i, "muts_dataset"]
	win_size = as.numeric(todo_here[i, "win_size"])
	row_id = todo_here[i, "id"]
		
	res = ADWGS_test(el, 
			gr_element_coords = gr_prepared_elements,
			gr_site_coords = GRanges(),
			gr_maf = all_prep_vars[[ds]],
			this_genome = this_genome,
			win_size = win_size)
			
	res = cbind(res, ds, row_id)
	res
}

cat(nrow(todo_here), ": ")
result = do.call(rbind, lapply(1:nrow(todo_here), 
		get_ADWGS, todo_here, gr_prepared_elements, all_prep_vars, this_genome))
fname = pff(c("ADWGS_tmp/result_iter_", iteration, "_job_", job_id, ".rsav"))
save(result, file = fname)