# -------------------------------------------
# -------------------------------------------
# 2. DATASET CLEANUP
# -------------------------------------------
# -------------------------------------------

# TODO: More advanced imputation methods

#' Missing Value Imputation
#'
#' Impute missing values in an Fc Array data frame on a per-group basis. Thus requires fc$group exists. This function can handle both single time point and longitudinal Fc Array data frames.
#' 
#' @param fc The Fc Array data frame. Can be single time point or longitudinal.
#' @param start_col Where you want the imputing to start iterating along columns from, presumably fccu_first_feat_col(fc).
#' @param type String denoting the type of imputation. Currently, only "det" for deterministic is suppoted.
#' @param method For calculating imputed values. Default is median.
#'
#' @return The Fc Array data frame, with missing values replaced based on method. By default, this will be the median of a group's values for that feature.
#'
#' @export
fcdc_impute = function(fc, start_col=fccu_first_feat_col(fc), type="det", method=median){
	groups = levels(fc$group)
	times = unique(fc$times)

	if (!is.null(times)){
		for (point in times){
			fc_time = fcff_only_attr(fc, point, times, "times")
			fc_time = fch_impute_fc_time(fc, start_col, method)

			# Insert the imputed values into the correct rows
			for (s in fc_time$subject){
				# Find the row in fc that corresponds to the row in fc_time
				indices = which(fc$subject == s)
				idx = indices[which(fc$times[indices] == point)]

				fc[idx,] = fc_time[which(fc_time$subject == s),]
			}
		}

	} else {
		fc = fch_impute_fc_time(fc, start_col, method)
	}

	return(fc)
}

#' Normalize An Fc Array
#'
#' Perform normalization on an Fc Array data frame. Just a glorified wrapper for scale() that calls fcdc_impute if necessary.
#'
#' @param fc The Fc Array data frame.
#' @param center A boolean or vector. See ?scale
#' @param scale A boolean or vector. See ?scale
#'
#' @return The Fc Array data frame, with feature vectors scaled/centered as specified.
#'
#' @export
fcdc_normalize = function(fc, center=TRUE, scale=TRUE){
	if (!all(complete.cases(fc))){
		fc = fcdc_impute(fc, fccu_first_feat_col(fc))
	}

	fc[,fccu_first_feat_col(fc):ncol(fc)] = scale(fc[,fccu_first_feat_col(fc):ncol(fc)], center, scale)

	return(fc)
}

#' Background Signal Subtraction
#'
#' Perform background signal subtraction from an Fc Array data frame.
#'
#' @param fc The Fc Array data frame.
#' @param base_ids The row indices in fc that correspond to the subjects who represent baseline (say, those in the PLACEBO group).
#'
#' @return The Fc Array data frame, with fc[base_ids,] having all feature measurements set to 0, and with non-baseline subjects having had their MFIs subtracted by the median of the baseline group's MFI values for each feature, and then scaled by the sd of those baseline MFI values.
#'
#' @export
fcdc_sub_back = function(fc, base_ids){
	for (i in fccu_first_feat_col(fc):ncol(fc)){
		baseline = fc[base_ids,i]
		measures = fc[!base_ids,i]

		# Anyone who is negative gets set to 0 (because negative MFI makes no sense)
		baseline[baseline < 0] = 0 
		background = median(baseline, na.rm=TRUE)
		back_sd = sd(baseline, na.rm=TRUE) 

		measures = measures - background
		measures[measures < 0] = 0

		# Everyone who is PLACEBO for a feature gets set to 0 
		fc[base_ids,i] = 0
		
		# Scale everyone according to the variation in the measurement itself (i.e. in baseline)
		measures = measures/back_sd 
		fc[!base_ids,i] = measures
	}

	return(fc)
}
