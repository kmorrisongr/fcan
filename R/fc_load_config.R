#' Load Fc Analysis Config .R File
#'
#' Loads the contents of the provided config file, performs data imputation, any specified filtering, binning, and prints out sessionInfo().
#'
#' @param config_location A string specifying the file path to the config.R file, relative to the location from which this function is called, or an absolute file path.
#' @param experiment_name A string specifying the name of the current experiment. Will be used to create a separate subdirectory for results.
#' @param fc The Fc Array data frame. NULL means that it will be loaded from the location specified in the config.R file. If an Fc Array data frame is passed, then it will be subjected to the filtering, binning, etc. that you have specified.
#' @param bin_method If flags$do_bin is TRUE, this is the function you will use to assign new groups.
#' @param ref_group If flags$do_differs is TRUE, this is a string specifying the group you want to define as baseline.
#'
#' @return A list containing the Fc Array data frame modified to your specifications, and the new results_dir string for this experiment, etc.
#'
#' @export
fc_load_config = function(config_location, experiment_name, flags, gopts, fc=NULL, only_group=NULL, bin_method=NULL, ref_group=NULL){
	backup = getwd()

	config_dir = strsplit(config_location,'/')[[1]]
	config_name = config_dir[length(config_dir)]
	config_dir = paste(config_dir[1:(length(config_dir)-1)], collapse='/')

	setwd(config_dir)
	source(config_name)

	results_dir = paste(results_dir, experiment_name, '/', sep='')
	dir.create(results_dir)

	if (exists("surv_location")){
		surv = read.csv(file=surv_location, stringsAsFactors=FALSE, header=TRUE)
		colnames(surv) = c("id", "group", "infected")
	}

	if (is.null(fc) && exists("fc_location")){
		fc = read.csv(file=fc_location, stringsAsFactors=FALSE, header=TRUE)
		fc$group = as.factor(fc$group)
	}

	# Even if you don't read an fc, fc=NULL being defined means this check will work
	if (!is.null(fc)){
		fc = fcdc_impute(fc)

		if (!is.null(flags$do_standard) && flags$do_standard){
			fc = fcdc_normalize(fc)
		}

		groups = levels(fc$group)
		times = unique(fc$times)

		filtered = fcff_wrap_filter(fc, results_dir, keep_filter, discard_filter, k_behavior, d_behavior)
		fc = filtered$fc
		results_dir = filtered$results_dir

		# RV144
		if (!is.null(flags$do_only_group) && !is.null(only_group) && flags$do_only_group){
			fc = fcff_only_attr(fc, only_group, groups, "group")
			results_dir = paste(results_dir, "only_", tolower(substr(only_group,1,4)), '/', sep='')
			dir.create(results_dir)
		}

		if (!is.null(flags$do_bin) && flags$do_bin){
			fc = bin_method(fc, samples, surv)
			results_dir = paste(results_dir, "binned/", sep='')
			dir.create(results_dir)
			base_ids = (fc[,"group"] == "yes")
		}

		if (!is.null(flags$do_differs) && flags$do_differs){
			if (!flags$do_bin){
				base_ids = (fc[,"group"] == ref_group)
			}
				
			diff_opts = list()
			diff_opts$base_ids = base_ids
			diff_opts$test = wilcox.test
			diff_opts$alternative = "either"
			diff_opts$adj_method = "fdr"
			filtered = fcff_wrap_filter(fc, results_dir, keep_filter, discard_filter, "differs", "NULL", diff_opts)
			fc = filtered$fc
			results_dir = filtered$results_dir
		}
	}

	stamp = gsub(':', '', Sys.time())
	stamp = gsub('-', '', stamp)
	stamp = gsub(' ', '_', stamp)
	sink(paste(results_dir, "sessionInfo_", stamp, ".txt", sep=''))
	print(Sys.time())
	print(sessionInfo())
	cat("\n=============================\n")
	print("Objects defined prior to config loading")
	cat("=============================\n")
	print(sapply(ls.str(1)[!(ls.str(1) %in% lsf.str(1))], get))
	sink()

	# If experiments are one level down
	if (!is.null(flags$adj_cwd) && flags$adj_cwd){
		results_dir = paste("../", results_dir, sep='')
	}

	setwd(backup)

	if (!is.null(fc)){
		return(list(fc=fc, results_dir=results_dir))

	} else {
		return(list(results_dir=results_dir))
	}
}
