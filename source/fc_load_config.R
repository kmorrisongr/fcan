# Loads fc_config.R
fc_load_config = function(config_location, experiment_name, fc=NULL, only_group=NULL, bin_method=NULL, ref_group=NULL){
	backup = getwd()

	config_dir = strsplit(config_location,'/')[[1]]
	config_name = config_dir[length(config_dir)]
	config_dir = paste(config_dir[1:(length(config_dir)-1)],collapse='/')

	setwd(config_dir)
	source(config_name)

	setwd(source_dir)
	source("fc_utils.R")
	source("custom.R")

	setwd(backup)
	setwd(config_dir)

	if (is.null(fc)){
		fc = read.csv(file=fc_location, stringsAsFactors=FALSE, header=TRUE)
		fc$group = as.factor(fc$group)
	}

	fc = fcu_impute(fc,fcu_first_feat_col(fc))

	surv = read.csv(file=surv_location, stringsAsFactors=FALSE, header=TRUE)
	colnames(surv) = c("id","group","infected")

	groups = levels(fc$group)
	times = unique(fc$times)

	results_dir = paste(results_dir, experiment_name, '/', sep='')
	filtered = fcu_wrap_filter(fc, results_dir, keep_filter, discard_filter, k_behavior, d_behavior)
	fc = filtered$fc
	results_dir = filtered$results_dir

	# RV144
	if (flags$do_only_group){
		fc = fcu_only_attr(fc, only_group, groups, "group")
		results_dir = paste(results_dir, "only_", tolower(substr(only_group,1,4)), '/', sep='')
		dir.create(results_dir)
	}

	if (flags$do_bin){
		fc = bin_method(fc, samples, surv)
		results_dir = paste(results_dir, "binned/", sep='')
		dir.create(results_dir)
		base_ids = (fc[,"group"] == "yes")
	}

	if (flags$do_differs){
		if (!flags$do_bin){
			base_ids = (fc[,"group"] == ref_group)
		}
			
		diff_opts = list()
		diff_opts$base_ids = base_ids
		diff_opts$test = wilcox.test
		diff_opts$alternative = "either"
		diff_opts$adj_method = "fdr"
		filtered = fcu_wrap_filter(fc, results_dir, keep_filter, discard_filter, "differs", "NULL", diff_opts)
		fc = filtered$fc
		results_dir = filtered$results_dir
	}

	sink(paste(results_dir,"sessionInfo.txt",sep=''))
	print(sessionInfo())
	sink()

	# If experiments are one level down
	if (flags$adj_cwd){
		results_dir = paste("../",results_dir,sep='')
	}

	setwd(backup)
	return(list(fc=fc,results_dir=results_dir))
}
