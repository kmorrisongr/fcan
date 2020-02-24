# -------------------------------------------
# -------------------------------------------
# MULTI-FC UTILITES
# -------------------------------------------
# -------------------------------------------

#' Generate Groupwise Combinations Of Fc Arrays
#'
#' Generate a list of data frames that each consist of only two treatment groups. Makes many group comparisons easier when paired with for loops.
#'
#' @param fc The Fc Array data frame.
#' @param order Determines what order the data frames are in the list
#' \itemize{
#'	\item "groups" (g1v2a, g1v2b, g1v3a, g1v3b, g2v3a, g2v3b, etc...).
#'	\item "times" (g1v2a, g1v3a, g2v3a, g1v2b, g1v3b, g2v3b, etc...).
#' }
#' @param num_groups The number of groups to represent in each entry in the returned list. For example, 2 would mean all pairwise combinations of groups would be represented.
#'
#' @return A list of Fc Array data frames, each with only num_groups represented. Returns NULL if it encountered an error.
#'
#' @export
fcmu_fcs_combs = function(fc, order, num_groups){
	groups = levels(fc$group)
	times = unique(fc$times)
	group_combs = combn(groups, num_groups)
	fcs = list()

	if (order == "groups"){
		if (!is.null(times)){
			# For each combination of groups - pairwise because fold change stuff
			# Won't change anything if only two groups
			for (i in 1:ncol(group_combs)){
				fc_comb = fch_fc_comb(fc, group_combs, group_combs[,i])

				for (t in 1:length(times)){
					idx = t + ((i-1) * length(times))
					fcs[[idx]] = fcff_only_attr(fc_comb, times[i], times, "times")
					names(fcs)[idx] = tolower(paste(group_combs[,i], times[i], collapse="_"))
				}
			}

		} else {
			for (i in 1:ncol(group_combs)){
				fcs[[i]] = fch_fc_comb(fc, group_combs, group_combs[,i])
				names(fcs)[i] = tolower(paste(group_combs[,i], collapse="_"))
			}
		}

	} else if (order == "times"){
		if (is.null(times)){
			print("No fc$times information present. Check your data frame.")
			return(NULL)

		} else {
			for (t in 1:length(times)){
				fc_time = fcff_only_attr(fc, times[i], times, "times")

				for (i in 1:ncol(group_combs)){
					idx = i + ((t-1) * ncol(group_combs))
					fcs[[idx]] = fch_fc_comb(fc_time, group_combs, group_combs[,i])
					names(fcs)[idx] = tolower(paste(group_combs[,i], times[i], collapse="_"))
				}
			}
		}
	}

	return(fcs)
}

#' Find y-limits For Many Data Frames
#'
#' Find the y-limits that will encompass all of the data present in the data frames in the list. Allows plots to be drawn with the same y-limits, making visual inspection and comparison easier and less dubious.
#'
#' @param data_list A list of data frames.
#' @param clm Which column in each data frame you are interested in.
#'
#' @return A vector of y-limits c(lower, upper). If, for some reason, no y-limits were found, you get back NULL and an error is written. If only a y-maximum was found, the y-minimum is set to 0 and you get back c(0, y_max).
#'
#' @export
fcmu_fcs_find_ylim = function(data_list, clm){
	y_min = NULL
	y_max = NULL

	for (d in data_list){
		measures = na.omit(d[,clm])

		if (length(measures) > 0){
			# Initialize
			if (is.null(y_min)){
				y_min = min(measures)

			} else if (min(measures) < y_min){
				y_min = min(measures)
			}

			# Initialize
			if (is.null(y_max)){
				y_max = max(measures)

			} else if (max(measures) > y_max){
				y_max = max(measures)
			}
		}
	}

	# If we didn't find anything
	if (is.null(y_min) & is.null(y_max)){
		ylim = NULL
		write("No limits found - check your data, fcmu_fcs_find_ylim(), and potentially fc_covb.R",file='')

	# If y_max but no y_min
	} else if (is.null(y_min)){
		ylim = c(0, y_max)
	}

	return(ylim)
}

#' Perform fcst_mannwhit On Many Fc Array Datasets
#'
#' For each feature in each Fc Array data frame, perform the Mann-Whitney U (Wilcoxon Rank-Sum) test to test the null hypothesis that a randomly selected value from one population (group) will be less than or greater than a randomly selected value from a second population (group). Here, this is to test whether the two sets of samples, which should be independent (treatment group A should not have an impact on treatment group B), were selected from populations having the same distribution, and there thus being no difference in the groups.
#'
#' @param fcs A list of Fc Array data frames. Each data frame should have just one time point and just two groups. The names of the entries in fcs should be a contraction of the groups you are comparing, separated by an underscore.
#' @param results_dir A string representing the directory where you want to store the results of your analyses.
#' @param adj_method The method you want to use to adjust p-values. Use "none" as a passthrough.
#' @param alternative What alternative hypothesis we're testing. Supports "two.sided", "one.sided", and "either" (result of "two.sided"/2).
#'
#' @return A list containing data frames for each corresponding Fc Array data frame in fcs. Each data frame has p-values for each feature (features are columns). 
#'
#' @export
fcmu_fcs_mannwhit = function(fcs, adj_method, alternative="two.sided"){
	if (!fch_fcs_assert_input(fcs)){
		warning("Returning NULL!")
		return(NULL)
	}

	output = list()
	n = 0

	for (i in 1:length(fcs)){
		fc = fcs[[i]]
		
		p_vals = fcst_mannwhit(fc, adj_method="none", alternative, exact)

		# MHC at the end
		n = n + length(p_vals)

		local = t(data.frame(p_vals=p_vals))
		colnames(local) = colnames(fc)[fccu_first_feat_col(fc):ncol(fc)]
		output[[i]] = local
		names(output[[i]]) = names(fcs)[i]
	}

	# Save time
	if (adj_method != "none"){
		for (i in 1:length(output)){
			output[[i]]["p_vals",] = p.adjust(output[[i]]["p_vals",], adj_method, n)
		}
	}

	return(output)
}

#' Perform fcst_fochs On Many Fc Array Datasets
#'
#' For each feature in each Fc Array data frame, calculate the fold change between the measurements of one group and another using method. Which group is the numerator is determined by the structure of comp_name.
#'
#' @param fcs A list of Fc Array data frames. Each data frame should have just one time point and just two groups. The names of the entries in fcs should be a contraction of the groups you are comparing, separated by an underscore. The first group in this string will be the numerator in the fold change ratios.
#' @param method What method will be used to calculate fold change. The default is the median of each group's values.
#'
#' @return A list containing data frames for each corresponding Fc Array data frame in fcs. Each data frame has fold change values for each feature (features are columns). 
#'
#' @export
fcmu_fcs_fochs = function(fcs, method=median){
	if (!fch_fcs_assert_input(fcs)){
		warning("Returning NULL!")
		return(NULL)
	}

	output = list()

	for (i in 1:length(fcs)){
		fc = fcs[[i]]
		current_name = names(fcs)[i]
		
		fochs = fcst_fochs(fc, current_name, method)
		
		local = t(data.frame(fochs=fochs))
		colnames(local) = colnames(fc)[first_feat_col(fc):ncol(fc)]
		output[[i]] = local
		names(output[[i]]) = names(fcs)[i]
	}

	return(output)
}
