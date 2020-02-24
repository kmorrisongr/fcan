# -------------------------------------------
# -------------------------------------------
# FEATURE FILTERING
# -------------------------------------------
# -------------------------------------------

#' Exclude An Attribute
#'
#' Exclude an attribute from an Fc Array data frame.
#'
#' @param fc The Fc Array data frame.
#' @param attr The value of the thing you want gone.
#' @param clm What column in fc attr will be found in.
#'
#' @return The Fc Array data frame, without the rows that have values in clm matching attr.
#'
#' @examples
#' fc_no_plac = fcff_exclude_attr(fc, "PLACEBO", "group")
#'
#' @export
fcff_exclude_attr = function(fc, attr, clm){
	toss = grepl(attr, fc[,clm], fixed=TRUE)

	return(fc[!toss,])
}

#' Filter Features By Difference From Baseline
#'
#' Filter features from an Fc Array data frame by their difference from baseline.
#'
#' @param fc The Fc Array data frame.
#' @param diff_opts A list containing the following sub-parameters:
#' \itemize{
#' 	\item diff_opts$base_ids Uhe row indices in fc that correspond to the subjects who represent baseline (say, those in the PLACEBO group).
#' 	\item diff_opts$test Which statistical test to use to determine difference from baseline.
#' 	\item diff_opts$alternative Alternative hypothesis (two.sided, either (two.sided/2)).
#' 	\item diff_opts$adj_method The method you want to use to adjust p-values. Use "none" as a passthrough.
#' }
#'
#' @return The p-values that represent the comparison between baseline and the other subjects using the test specified in diff_opts.
#'
#' @examples
#' p_vals = fcff_feats_diff(fc, diff_opts)
#' keep = p_vals < 0.05
#'
#' @export
fcff_feats_diff = function(fc, diff_opts){
	# No categorical variables
	fc = data.frame(fc[,fccu_first_feat_col(fc):ncol(fc)])

	p_vals = as.vector(apply(fc, 2, function(x, base_ids, test, alternative){
			       placebo = x[base_ids]
			       vaccine = x[!base_ids]

				# Is feature different?
				if (alternative == "either"){
					t = test(placebo, vaccine, alternative="two.sided")
					t$p.value = t$p.value / 2

				} else {
					t = test(placebo, vaccine, alternative=alternative)
				}

				return(t$p.value)
		       }, base_ids=diff_opts$base_ids, test=diff_opts$test, alternative=diff_opts$alternative))

	p_vals = p.adjust(p_vals, method=diff_opts$adj_method)

	return(p_vals)
}

#' Filter Features
#'
#' Filter features based on keyword. You can keep those, remove those, and many other things.
#'
#' @param fc The Fc Array data frame.
#' @param feats Vector of strings for the features you want to keep or remove.
#' @param action String specifying whether you want to "keep" or "toss" the features in feats.
#' @param behavior Whether the feature names will be matched by:
#' \itemize{
#' 	\item containing any of the strings in feats ("permissive").
#'	\item containing all of the strings in feats ("strict").
#'	\item Alternatively, "differs" specifies that features will be removed by difference from baseline. This will be done instead of keyword filtering (so to do both, call this function twice).
#' }
#' @param diff_opts list of options if behavior = "differs". see ?fcff_feats_diff.
#'
#' @return The Fc Array data frame filtered based on what you specified.
#'
#' @examples
#' new_fc = fcff_filter_features(fc, c("gp41", "gp120"), "keep")
#' my_diff_opts = list() # fill with stuff
#' diff_fc = fccu_filter_features(fc, NULL, NULL, behavior="differs",
#' 					diff_opts=my_diff_opts)
#'
#' @export
fcff_filter_features = function(fc, feats, action, behavior="permissive", diff_opts=NULL){
	if (!(length(feats) > 0) && is.null(diff_opts)){
		return(fc)
	}
	local = data.frame(fc[,fccu_first_feat_col(fc):ncol(fc)])
	colnames(local) = colnames(fc)[fccu_first_feat_col(fc):ncol(fc)]
	subjects = fc$subject

	keep = vector("logical")

	if ( (behavior == "differs") && (!is.null(diff_opts)) ){
		p_vals = fcff_feats_diff(local, diff_opts)
		keep = p_vals < 0.05

	} else {
		# Written this way to accommodate & and |
		for (i in 1:length(feats)){
			i_keep = grepl(feats[i], colnames(local), fixed=TRUE)

			if (i == 1){
				keep = i_keep
				next
			}

			# ?ifelse to understand rep usage
			if (length(i_keep) > 0){
				keep = ifelse(rep(behavior == "strict", length(keep)), keep & i_keep, keep | i_keep)
			}
		}
	}

	indices = ifelse(rep(action == "toss", length(keep)), !keep, keep)
	output = as.data.frame(local[,indices])
	colnames(output) = colnames(local)[indices]

	if (fccu_first_feat_col(fc) != 1){
		output = cbind(as.data.frame(fc[,1:(fccu_first_feat_col(fc)-1)]), output)
		colnames(output)[1:(fccu_first_feat_col(fc)-1)] = colnames(fc)[1:(fccu_first_feat_col(fc)-1)]
	}

	if (!is.null(output$group)){
		output$group = factor(output$group)
	}

	return(output)
}

#' Keep Only Specific Values Of A Column
#'
#' A wrapper for fcff_exclude_attr that performs the exclusion behavior repeatedly in order to produce a dataset with just attr in attrs.
#'
#' @param fc The Fc Array data frame.
#' @param attr The value of the thing you want gone.
#' @param attrs A vector containing all the values that attr can take on in clm.
#' @param clm What column in fc attr will be found in.
#'
#' @return The Fc Array data frame, with only the rows that have values in clm matching attr.
#'
#' @examples
#' groups = c("PLACEBO", "VACCINE")
#' fc_only_placebo = fcff_only_attr(fc, "PLACEBO", groups, "group")
#'
#' @export
fcff_only_attr = function(fc, attr, attrs, clm){
	output = fc
	for (a in attrs){
		if (a != attr){
			output = fcff_exclude_attr(output, a, clm)
		}
	}

	return(output)
}

#' Remove Correlated Features
#'
#' Remove correlated features from fc based on cutoff.
#'
#' @param fc The Fc Array data frame.
#' @param cutoff The correlation coefficient cutoff above which we will remove features from fc.
#' @param behavior How we'll decide what to keep as a representative feature from those to be removed. "sample" randomly selects one feature. "first" picks the first one in the list.
#'
#' @return The Fc Array data frame filtered based on what you specified.
#'
#' @export
fcff_remove_cor = function(fc, cutoff, behavior="sample"){
	# No categorical variables
	output = data.frame(fc[,fccu_first_feat_col(fc):ncol(fc)])

	cor_mat = cor(output)
	high_cor = findCorrelation(cor_mat, cutoff=cor_cutoff)

	# The features should, by definition, be interchangeable
	if (behavior == "sample"){
		keep = sample(high_cor, 1)
		keep = which(high_cor == keep)

	} else if (behavior == "first"){
		keep = high_cor[1]
	}

	high_cor = high_cor[-keep]
	output = output[,-high_cor]

	if (fccu_first_feat_col(fc) != 1){
		output = cbind(as.data.frame(fc[,1:(fccu_first_feat_col(fc)-1)]), output)
		colnames(output)[1:(fccu_first_feat_col(fc)-1)] = colnames(fc)[1:(fccu_first_feat_col(fc)-1)]
	}

	return(output)
}

#' Higher-Level Interface For Filtering Features
#'
#' Filter features in such a way as to set up experiments. Will filter feautres/perform fcff_feats_diff as necessary, and will create files/directories (as necessary) to facilitate result reproduction and keep different versions of fc analysis sequestered.
#'
#' @param fc The Fc Array data frame
#' @param results_dir A string representing the directory where you want to store the results of your analyses.
#' @param keep_filter A vector of strings of features that you want to retain in fc. If no keeping is to be done, pass an empty vector c().
#' @param discard_filter A vector of strings of features that you want to remove from fc. If no discarding is to be done, pass an empty vector c().
#' @param k_behavior A string "permissive" or "strict" controlling keeping behavior. Can also be set to "differs" to perform fcff_feats_diff. See ?fcff_filter_features for more.
#' @param d_behavior A string "permissive" or "strict" controlling discarding behavior. Can also be set to "differs" to perform fcff_feats_diff. See ?fcff_filter_features for more.
#' @param diff_opts A list of subparameters for fcff_feats_diff. Pass "NULL" if you don't want to do any fcff_feats_diff things. See ?fcff_feats_diff for more.
#'
#' @return A list containing the filtered Fc Array data frame (accessed by $fc), and the new results_dir string based on the filtering performed (accessed by $results_dir). If both keeping and discarding are performed, results_dir will be of the form complex_YYYYMMDD_HHMMSS, and the directory will contain a txt file with the filtering that was performed.
#'
#' @examples
#' results_dir = "/home/me/science/experiments/"
#' filtered = fcff_wrap_filter(fc, results_dir, c("gp41", "gp120"), c(), "permissive", "NULL")
#' fc = filtered$fc
#' # is now "/home/me/science/experiments/only_gp41_gp120"
#' results_dir = filtered$results_dir
#'
#' @export
fcff_wrap_filter = function(fc, results_dir, keep_filter, discard_filter, k_behavior, d_behavior, diff_opts="NULL"){
	cat("Fc Array has", ncol(fc)-(fccu_first_feat_col(fc)-1), "features", '\n')

	if (k_behavior == "differs" | d_behavior == "differs"){
		fc = fcff_filter_features(fc, keep_filter, "keep", "differs", diff_opts)
		results_dir = paste(results_dir, "differs", '/', sep='')

	} else if (length(keep_filter) > 0 & length(discard_filter) > 0){
		stamp = toString(Sys.time())
		stamp = paste(strsplit(stamp,' ')[[1]], collapse='_')
		stamp = gsub('-','',stamp)
		stamp = gsub(':','',stamp)

		results_dir = paste(results_dir, "complex_", stamp, '/', sep='')
		dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

		fc = fcff_filter_features(fc, keep_filter, "keep", k_behavior)
		fc = fcff_filter_features(fc, discard_filter, "toss", d_behavior)

		filter_file = paste(results_dir, "feat_filter.txt", sep='')
		write(paste("only_", keep_filter, sep=''), filter_file)
		write(paste("no_", discard_filter, sep=''), filter_file, append=TRUE)

	} else if (length(keep_filter) > 0){
		fc = fcff_filter_features(fc, keep_filter, "keep", k_behavior)
		results_dir = paste(results_dir, "only_", paste(keep_filter, collapse="-"), '/', sep='')

	} else if (length(discard_filter) > 0){
		fc = fcff_filter_features(fc, discard_filter, "toss", d_behavior)
		results_dir = paste(results_dir, "no_", paste(discard_filter, collapse="-"), '/', sep='')
	}

	dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

	cat("Filtered Fc Array has", ncol(fc)-(fccu_first_feat_col(fc)-1), "features", '\n')
	return(list(fc=fc, results_dir=results_dir))
}
