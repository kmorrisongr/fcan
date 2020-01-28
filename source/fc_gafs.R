# Written by Kyle Morrison for analysis of Fc Array Data
# This function builds a random forest classification model using genetic algorithm feature selection

library(doParallel)
library(caret)

fc_gafs = function(fc, groups, stats_dir, folds=10, cor_cutoff=0.9, min_feat_scale=0.25, iters=50, popsize=50, elite=0, pmutation=0.2, pcrossover=0.8, do_permute=FALSE, do_cor=FALSE, do_penalize=FALSE){
	# --------------------------------------------------
	# Environment setup
	# --------------------------------------------------

	feat_sel_dir = paste(stats_dir, "gafs/", sep="")
	dir.create(feat_sel_dir)

	# My laptop has 6 physical cores
	cl = makePSOCKcluster(5)
	registerDoParallel(cl)
	set.seed(1337)


	# --------------------------------------------------
	# Data prep, dealing with correlation and balancing
	# --------------------------------------------------

	local = fcu_impute(fc, fcu_first_feat_col(fc))

	if (do_permute){
		feat_sel_dir = paste(stats_dir, "gafs/perm/", sep='')
		dir.create(feat_sel_dir)

		# Shuffle the class labels for permutation control
		local[,"group"] = sample(local[,"group"])
		cat("Shuffling the class labels\n")
	}

	if (do_cor){
		local = remove_cor(local, cor_cutoff, "sample")
		cat("Removing highly-correlated features (>", cor_cutoff, ")\n")
	}

	fcs = fcu_fcs_combs(local, groups)

	# Do balancing as necessary
	group_counts = vector("numeric")
	for (g in groups){
		group_counts = c(group_counts, sum(local[,"group"] == g))
	}
	max_subjs = min(group_counts)
	# This is not written to support >2 groups right now
	min_idx = which(group_counts == max_subjs)
	other_idx = which(group_counts != max_subjs)


	# --------------------------------------------------
	# Balance and/or holdout data
	# --------------------------------------------------

	if (group_counts[1] != group_counts[2]){
		cat("Balancing the groups\n")

		# Randomly drop the extra members of fcs[[other_idx]]
		# For now, do nothing with the leftover data
		balanced_data = rbind(fcs[[min_idx]], fcs[[other_idx]][sample(1:max_subjs, max_subjs),])

		group = balanced_data[,"group"]
		group = as.factor(group)
		balanced_data = balanced_data[,-(1:2)]

		# Run gafs on balanced data - should not need my own CV
		ga = fcu_wrap_gafs(balanced_data, group, min_feat_scale, iters=iters, popsize=popsize, elite=elite, pmutation=pmutation, pcrossover=pcrossover, do_penalize)
		print(ga)
		final_feats = ga$ga$final
		write(final_feats, file=paste(feat_sel_dir, "final_feats.txt", sep=''))


	# --------------------------------------------------
	# Use everything for building the model
	# --------------------------------------------------

	} else {
		keep = 1:nrow(fcs[[other_idx]])
		train = rbind(fcs[[min_idx]], fcs[[other_idx]][keep,])

		group = train[,"group"]
		group = as.factor(as.character(group))
		train = train[,-1]

		if (do_cor){
			train = fcu_remove_cor(train, cor_cutoff, "sample")
			cat("Removing highly-correlated features (>", cor_cutoff, ")\n")
		}

		# --------------------------------------------------
		# Run gafs on each train, get score

		ga = fcu_wrap_gafs(train, group, min_feat_scale, iters=35, elite=0, pmutation=0.2, pcrossover=0.8, do_penalize)
		print(ga)
		final_feats = ga$ga$final
		write(final_feats, file=paste(feat_sel_dir, "final_feats.txt", sep=''))

		cat("All data used to train - couldn't do any testing\n")
	}
}
