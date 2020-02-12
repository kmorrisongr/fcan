# -------------------------------------------
# -------------------------------------------
# 9. HELPER FUNCTIONS
# -------------------------------------------
# -------------------------------------------

#' Check The Validity Of Fc Array
#'
#' Confirm that an Fc Array data frame has only 2 groups and only one time point.
#'
#' @param fc The Fc Array data frame.
#'
#' @return A boolean determining whether the Fc Array passed inspection.
#'
#' @export
fch_assert_input = function(fc){
	if (is.null(fc$group) || length(levels(fc$group)) != 2){
		warning("Your Fc Array data frame is supposed to have 2 group levels!")
		return(FALSE)

	} else if (length(unique(fc$times)) > 1){
		warning("Your Fc Array data frame is supposed to be single time point!")
		return(FALSE)
	}
	
	return(TRUE)
}

#' Calculate Simple Change Over Baseline
#'
#' @param base_measure The value for this data point at baseline
#' @param measure The value for this data point at some later time point
#'
#' @return A single value: measure/base_measure
#'
#' @export
fch_covb_simple = function(base_measure, measure){
	return(measure/base_measure)
}

#' Generate A Single Fc Array Data Frame Of Just Two Groups
#'
#' @param fc The Fc Array data frame. Must be a single time point.
#' @param group_combs The group combination matrix obtained from combn(groups, 2).
#' @param to_comp The column in group_combs that is to be kept.
#'
#' @return The Fc Array data frame with just the two groups in to_comp.
#'
#' @export
fch_fc_comb = function(fc, group_combs, to_comp){
	# Exclude the other groups
	to_exclude = setdiff(group_combs, to_comp)

	for (j in to_exclude){
		fc = fcff_exclude_attr(fc, j, "group")
	}

	return(fc)
}

#' Check The Validity Of Many Fc Arrays
#'
#' Confirm that each Fc Array in the list has only 2 groups, only one time point, and the list has been given names.
#'
#' @param fc The Fc Array data frame.
#'
#' @return A boolean determining whether the Fc Array passed inspection.
#'
#' @export
fch_fcs_assert_input = function(fcs){
	if (is.null(names(fcs))){
		warning("Your Fc Array data frames list is supposed to have names that are contractions of the groups to be compared!")
		return(FALSE)
	}

	for (fc in fcs){
		if (!fch_assert_input(fc)){
			return(FALSE)
		}
	}

	return(TRUE)
}

#' Caluclate Geometric Mean
#'
#' @param vec A vector of numeric values.
#'
#' @return A single value: the nth root of the product of vec
#'
#' @export
fch_geom_mean = function(vec){
	vec = na.omit(vec)
	return( prod(vec)**(1/length(vec)) )
}

#' Missing Value Imputation For A Single Time Point
#'
#' @param fc The Fc Array data frame, but only at a single time point.
#' @param start_col Where you want the imputing to start iterating along columns from, presumably fccu_first_feat_col(fc).
#' @param method For calculating imputed values. Default is median.
#'
#' @return The Fc Array data frame, with missing values replaced based on method. By default, this will be the median of a group's values for that feature.
#'
#' @export
fch_impute_fc_time = function(fc, start_col, method=median){
	groups = levels(fc$group)

	for (i in start_col:ncol(fc)){
		for (g in groups){
			current = (fc$group == g)
			measures = fc[current,i]

			measures[is.na(measures)] = method(measures, na.rm=TRUE)
			measures[measures < 0] = 0

			fc[current,i] = measures
		}
	}

	return(fc)
}

#' Make A Reagant/Antigen Plot Legend
#'
#' Add a reagant/antigen legend to a plot.
#'
#' @param gopts A list of graphical options. Must have specified the subparameters:
#' \itemize{
#'	\item gopts$reag_legend_inset The inset for the legend relative to the plot.
#'	\item gopts$reag_cats The categories for your detection reagants.
#'	\item gopts$reag_cols The colors for the detection reagant categories.
#'	\item gopts$ant_shapes The shapes for the antigens being detected.
#' }
#'
#' @return Nothing!
#'
#' @export
fch_make_reag_legend = function(gopts){
	# All the rep calls are so we can draw the reagant colors and the antigen shapes together
	legend("topright", inset=gopts$reag_legend_inset, cex=0.85, xpd=TRUE,
	       legend=c(gopts$reag_cats, "other", gopts$feat_ant_legend, "other"),

	       bty="n", text.col=c(gopts$reag_cols, rep( "black",length(gopts$ant_shapes)+1 )),

	       col=c(gopts$reag_cols, rep( "black",length(gopts$ant_shapes) )),
	       pt.bg=c(gopts$reag_cols, rep( "white",length(gopts$ant_shapes) )),

	       # Just make the shapes for the reagant colors circles
	       pch=c(rep(21,length(gopts$reag_cols)), c(gopts$ant_shapes)))
}

#' Calculate Many Min-Max Z-Scores
#'
#' @param vec A vector you wish to min-max z-scale.
#'
#' @return A vector of min-max z-scores.
#'
#' @export
fch_min_max_feat = function(vec){
	vec = na.omit(vec)
	return(as.vector(sapply(vec, fch_z_score, vec=vec)))
}

#' Calculate A Single Min-Max Z-Score
#'
#' @param x A single data point.
#' @param vec The vector from which x derives.
#'
#' @return A single min-max z-score for x.
#'
#' @export
fch_z_score = function(x, vec){
	vec = na.omit(vec)
	return( (x - min(vec))/(max(vec) - min(vec)) )
}
