# -------------------------------------------
# -------------------------------------------
# STATISTICAL TESTING
# -------------------------------------------
# -------------------------------------------

# TODO: Check for presence/absence of time points

#' Calculate Change Over Baseline
#'
#' Calculate change over baseline for each subject using a provided method.
#'
#' @param fc The Fc Array data frame, in longitudinal format.
#' @param method The method used to calculate change over baseline. The default is a simple division.
#'
#' @return The Fc Array data frame, in longitudinal format (sans the first time point, which is baseline), with all feature MFI values replaced by change over baseline values as calculated by method.
#'
#' @export
fcst_covb = function(fc, method=fch_covb_simple){
	times = unique(fc$times)
	subjects = unique(fc$subject)

	baseline = fcff_only_attr(fc, times[1], times, "times")
	output = fcff_exclude_attr(fc, times[1], "times")

	for (s in subjects){
		idx = which(output$subject == s)[1]

		# If subject present after baseline
		if (!is.na(idx)){
			for (t in times[-1]){
				for (i in fccu_first_feat_col(output):ncol(output)){
					# No check, because presumably all subjects present at baseline
					base_measure = baseline[which(baseline$subject == s),i]
					output[idx,i] = method(base_measure, output[idx,i])
				}

				next_s = output$subject[idx + 1]

				# If next row is same subject
				if (!is.na(next_s) && next_s == s){
					idx = idx + 1

				} else {
					break
				}
			}
		}
	}

	return(output)
}

#' Perform Mann-Whitney U Test On Fc Array Features
#'
#' For each feature in an Fc Array data frame, perform the Mann-Whitney U (Wilcoxon Rank-Sum) test to test the null hypothesis that a randomly selected value from one population (group) will be less than or greater than a randomly selected value from a second population (group). Here, this is to test whether the two sets of samples, which should be independent (treatment group A should not have an impact on treatment group B), were selected from populations having the same distribution, and there thus being no difference in the groups.
#'
#' @param fc The Fc Array data frame. Must have two groups present.
#' @param adj_method The method you want to use to adjust p-values. Use "none" as a passthrough.
#' @param alternative What alternative hypothesis we're testing. Supports "two.sided", "one.sided", and "either" (result of "two.sided"/2).
#' @param exact A boolean to control whether or not exact p-values are calculated.
#'
#' @return A vector of p-values of the same length as fccu_first_feat_col(fc):ncol(fc).
#'
#' @export
fcst_mannwhit = function(fc, adj_method, alternative="two.sided", exact=FALSE){
	p_vals = vector("numeric")

	if (!fch_assert_input(fc)){
		warning("Returning NULL!")
		return(NULL)
	}

	for (i in fccu_first_feat_col(fc):ncol(fc)){
		# Remove any NAs from this measurement data
		measure = as.numeric(na.omit(fc[,i]))
		remove = is.na(fc[,i])
		if (length(remove) > 0){
			group = fc$group[!remove]
		} else {
			group = fc$group
		}

		if (alternative == "either"){
			test = wilcox.test(measure ~ group, alternative="two.sided", exact=exact)
			p_val = c(p_vals, (test$p.value)/2)

		} else {
			test = wilcox.test(measure ~ group, alternative=alternative, exact=exact)
			p_val = c(p_vals, (test$p.value))
		}
	}

	return(p_vals)
}

#' Calculate Fold Change For Fc Array Features
#'
#' For each feature in an Fc Array data frame, calculate the fold change between the measurements of one group and another using method. Which group is the numerator in the fold change ratio is determined by the structure of comp_name.
#'
#' @param fc The Fc Array data frame. Must have two groups present.
#' @param comp_name A string representing the name of the group comparison. This should be a contraction of the groups you are comparing, separated by an underscore. The first group in this string will be the numerator in the fold change ratios.
#' @param method What method will be used to calculate fold change. The default is the median of each group's values.
#'
#' @return A vector of fold change values of the same length as fccu_first_feat_col(fc):ncol(fc).
#'
#' @export
fcst_fochs = function(fc, comp_name, method=median){
	fochs = vector("numeric")

	if (!fch_assert_input(fc)){
		warning("Returning NULL!")
		return(NULL)
	}

	for (i in fccu_first_feat_col(fc):ncol(fc)){
		# Remove any NAs from this measurement data
		measure = as.numeric(na.omit(fc[,i]))
		remove = is.na(fc[,i])
		if (length(remove) > 0){
			group = fc$group[!remove]
		} else {
			group = fc$group
		}

		# placebo_vaccine -> [1] "placebo" [2] "vaccine"
		name = strsplit(comp_name,'_')[[1]]
		first_name = name[1]
		second_name = name[2]

		id_first = which(tolower(as.character(group)) == tolower(first_name))
		id_second = which(tolower(as.character(group)) == tolower(second_name))
		measures_first = measure[id_first]
		measures_second = measure[id_second]

		fochs = c(fochs, method(measures_first)/method(measures_second))
	}

	return(fochs)
}
