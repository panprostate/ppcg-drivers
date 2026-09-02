source("/PATH_TO_WORKING_DIR/bin/DATE_TAG/000_HEADER.R")

iteration = 0
n_jobs_total = 2500
homedir = paste0('/PATH_TO_WORKING_DIR/bin/', gsub("\\.", "", pff("")))
mem_lmt = "10G"

if (iteration > 0) {
	
	load(file = pff(paste0("ADWGS_tmp/todo_list.", iteration - 1, ".rsav")))
	all_row_ids = todo_list[,1]
	
	fnames = unique(unlist(lapply(0:(iteration-1), function(iter)
		list.files(pff("ADWGS_tmp"), pattern = paste0("result_iter_", iter), full.names = TRUE))))
	done_row_ids = unique(do.call(c, lapply(fnames, function(fname) {
		x <- try(load(fname), silent = TRUE); 
		if (class(x)=="try-error") return(NULL)
		if (nrow(result) == 0) {
			return(NULL)
		}
		result[, "row_id"] 
	})))

	rows_to_do = setdiff(all_row_ids, done_row_ids)
	todo_list = todo_list[todo_list[,1] %in% rows_to_do,]
	save(todo_list, file = pff(c("ADWGS_tmp/todo_list.", iteration, ".rsav")))	
	n_jobs_total = min(n_jobs_total, nrow(todo_list))
	cat("\n\nROWS TO DO: ", nrow(todo_list), "\n\n")
}





cat('\n\n\n\nhomedir=', homedir, '; n_jobs_total=', n_jobs_total, '; iteration=', iteration, ' \n
	for i in `seq 1 ', n_jobs_total ,'`; do \n
	touch ${homedir}/log/002b_run_ADWGS.R.$i.log
	echo "module load r-bl; Rscript ${homedir}/002b_run_ADWGS.R $iteration $n_jobs_total $i > ${homedir}/log/002b_run_ADWGS.R.$i.log 2>&1 " | qsub -V -N ADWGS -l h_vmem=', mem_lmt, ' -o /dev/null -e /dev/null \n
	done\n\n\n\n', sep='')
