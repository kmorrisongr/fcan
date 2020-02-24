# -------------------------------------------
# -------------------------------------------
# GRAPHICAL OUTPUT
# -------------------------------------------
# -------------------------------------------

#' Generate Boxplots
#'
#' Generate boxplots for a subset of features in an Fc Array data frame, separated out by group membership. Plots both boxes and whiskers and jitter plot of the points. Supports more than 2 groups. Resultant plots will live in a pdf, and be saved in results_dir.
#'
#' @param fc The Fc Array data frame
#' @param to_box A vector of strings of feature names that you wish to generate boxplots for.
#' @param results_dir Where you want your results to be saved.
#' @param flags A list of analysis booleans. Must have specified the subparameters:
#' \itemize{
#'	\item flags$do_same_range Do you want all boxplots to have the same y-limits?
#'	\item flags$nicer_legend Do you want a nicer legend? Currently applies only to the hardcoded case of RV144, where the binning generates groups that are "yes" and "no". If this flag is TRUE, these get replaced in the plot legends with "Infected" and "Uninfected".
#' }
#' @param gopts A list of graphical options. Must have specified the subparameters, and then the subsequent optional parameters:
#' \itemize{
#'	\item gopts$bpar_mar The margins for the plots. Sensible starting place is c(5.1, 4.1, 4.1, 4).
#'	\item gopts$jplot_group_cols The group colors for the jitter plot points that get drawn on top of the boxes. Should be the same length as unique(fc$group) (i.e. the number of groups you have).
#'	\item gopts$bplot_group_cols The same as jplot_group_cols, but for the boxes themselves. Will need to be slightly different if you want points to show up. If one of jplot_group_cols is black, for example, the corresponding bplot_group_cols should be white. Experiment to get it right.
#'	\item gopts$box_legend_inset The inset for the legend in the boxplots. Sensible starting place is c(-0.40, 0).
#' }
#' \itemize{
#'	\item gopts$pdf_name_suffix If specified, will change the output pdf name from "group1_group2_boxplots.pdf" to "group1_group2_[pdf_name_suffix]_boxplots.pdf".
#'	\item gopts$main If specified, sets the title of all the boxplots to the provided string. Else, the name of the feature is used.
#'	\item gopts$ylab If specified, sets the title of all boxplots to the provided string. Else, "MFI" will be used.
#' }
#'
#' @return Nothing! Just outputs plots.
#'
#' @examples
#' flags = list() # add some things
#' gopts = list() # add some more things
#' # produce a pdf containing boxplots of all IgG3 features, with points colored by fc$group membership
#' fcgo_boxplot(fc, colnames(fc)[grepl("IgG3", colnames(fc), fixed=TRUE)], "/home/me/results/", flags, gopts)
#'
#' @export
fcgo_boxplot = function(fc, to_box, results_dir, flags, gopts, seed=1337){
	# For jittery locationy
	set.seed(seed)

	fc = fcff_filter_features(fc, to_box, "keep")
	fc = fcdc_impute(fc, fccu_first_feat_col(fc))
	feats = sort(to_box)

	groups = levels(fc$group)
	times = unique(fc$times)
	num_groups = length(groups)


	# -------------------------------------------
	# Dataframe list broken down by group, then time point
	# g1a, g1b, g2a, g2b, g3a, g3b, etc...
	# -------------------------------------------
	df_list = fcmu_fcs_combs(fc, "groups", num_groups=1)

	if (!is.null(gopts$pdf_name_suffix)){
		pdf_name = paste(results_dir, paste(tolower(groups), collapse='_'),
				 '_', gopts$pdf_name_suffix, "_boxplots.pdf", sep='')

	} else {
		pdf_name = paste(results_dir, paste(tolower(groups), collapse='_'), "_boxplots.pdf", sep='')
	}
	pdf(pdf_name)
	# Use default if not provided
	if (is.null(gopts$bpar_mar)){
		gopts$bpar_mar = c(5, 4, 4, 2)
	}
	par(mar=gopts$bpar_mar, mfrow=c(2, 2))


	# -------------------------------------------
	# For each feature,  make a boxplot
	# -------------------------------------------

	for (f in feats){
		if (!is.null(flags$do_same_range) && flags$do_same_range){
			ylim = fcmu_fcs_find_ylim(master_foch, f)

		} else {
			ylim = NULL
		}

		# -------------------------------------------
		# Get all the measurements for this feature

		measures = vector("numeric")
		# For proper grouping of measurements for plotting
		box_groups = vector("numeric")
		for (j in 1:length(df_list)){
			local = df_list[[j]]
			# Remove NAs from the measurements
			to_merge = as.vector(na.omit(local[, f]))
			measures = c(measures, to_merge)
			box_groups = c(box_groups, rep(j, times=length(to_merge)))
		}
		
		#measures = log10(measures)


		# -------------------------------------------
		# Grouping and coloring

		# at tells boxplot how to color/group the boxplots
		box_at = vector("numeric")
		for (j in 1:num_groups){
			to_add = seq(j, (length(df_list)+j),  by=(num_groups+1))
			box_at = c(box_at, to_add)
		}
		if (!is.null(times) && (length(times) == 1)){
			new_box_at = vector("numeric")
			for (j in 1:num_groups){
				new_box_at = c(new_box_at, box_at[j + ((j-1)*(num_groups-1))])
			}
			box_at = new_box_at
		}


		# -------------------------------------------
		# Plot!

		# Might have to write your own pretty_title,  depending on what you want
		#main = cust_ed608_pretty_title(i)
		main = ifelse(is.null(gopts$main), f, gopts$main)

		# ?ifelse to understand [1]
		xlab = ifelse(!grepl("dummy", times, fixed=TRUE), "Time Point", '')[1]

		# y axis labels
		if (grepl("score", f, fixed=TRUE)){
			ylab = "Score"

		} else {
			ylab = ifelse(is.null(gopts$ylab), "MFI", gopts$ylab)
		}

		boxplot(measures ~ box_groups, at=box_at, names=NA, col="white", 
			main=main, outline=FALSE, xlab=xlab, ylab=ylab, ylim=ylim)
		if (!is.null(times)){
			each = length(times)
			time_mod = length(times)
			labels = times

		} else {
			each = 1
			time_mod = 1
			labels = ''
		}
		stripchart(measures ~ box_groups, at=box_at, vertical=TRUE, method="jitter", add=TRUE, pch=20, 
			   col=rep(gopts$jplot_group_cols, each=each))

		axis(1, at=seq((0.5*(num_groups+1)), (num_groups+1)*time_mod, by=num_groups+1), labels=labels)

		# print(groups) before reassigning to make sure they line properly - no reversed findings!
		if (!is.null(flags$nicer_legend) && flags$nicer_legend){
			groups = c("Uninfected", "Infected")
		}
		if (is.null(gopts$box_legend_inset)){
			gopts$box_legend_inset = c(-0.40, 0)
		}
		legend("topright", inset=gopts$box_legend_inset, legend=groups, fill=gopts$bplot_group_cols, cex=0.7, xpd=TRUE)
	}

	dev.off()
}

# A function for getting the properties for a feature (reagant or antigen)

#' Get Properties For A Feature
#'
#' Get properties (shape, color) for a feature (case-insensitive) from the graphical presets you've specified.
#'
#' @param feat The feature string you want to get something for. Case-insensitive. Matching is determined by substrings, so if you want something to be separated out (say, IgG subclasses) while calling this multiple times, put the subclasses first.
#' @param prop_cats The categories you have split features into.
#' @param prop_quals The shape, color, etc. you are querying into.
#'
#' @return A singular value from prop_quals corresponding to which entry in prop_cats feat corresponds to. If no match is found, the last entry in prop_quals is returned (useful for an "other" category if you don't want colors, etc. for every single type of feature).
#'
#' @examples
#' gopts$reag_cats = c("IgG", "FcgRIIa", "FcgRIIb", "FcgRIIIa", "FcgRIIIb")
#' # orchid is for anything that's not in the above categories
#' gopts$reag_cols = c("black", "forestgreen", "purple", "gold", "orange", "orchid")
#' feats = c("IgG3.gp120", "FcgRIIa.gp41", "FcgRIIIa.gp70.V1.V2", "C1q.aHuIgG")
#' (colors = as.vector(sapply(feats, fcgo_get_prop, gopts$reag_cats, gopts$reag_cols)))
#' [1] "black" "forestgreen" "gold" "orchid"
#'
#' # Subclasses
#' gopts$reag_cats = c("IgG1", "IgG2", "IgG3", "IgG4", "C1q", "IgG")
#' gopts$reag_cols = c("blue", "pink", "turquoise1", "saddlebrown", "green", "black", "orchid")
#' feats = c("IgG3.gp120", "IgG.gp41.MN", "Ig2.Vif", "FcgRIIa.gp41", "FcgRIIIa.gp70.V1.V2", "C1q.aHuIgG")
#' (colors = as.vector(sapply(feats, fcgo_get_prop, gopts$reag_cats, gopts$reag_cols)))
#' [1] "black" "turquoise1" "black" "pink" "orchid" "orchid" "green"
#'
#' @export
fcgo_get_prop = function(feat, prop_cats, prop_quals){
	feat = tolower(feat); prop_cats = tolower(prop_cats);

	for (i in 1:length(prop_cats)){
		if (any(grepl(prop_cats[i], feat, fixed=TRUE))){
			return(prop_quals[i])
		}

		# if this prop_cats has commas, split it, and then sub-check
		split_feat = strsplit(prop_cats[i], ',')[[1]]
		if (length(split_feat) > 1){
			# For any of these, return the same ant_shape
			for (sub_feat in split_feat){
				if (any(grepl(sub_feat, feat, fixed=TRUE))){
					return(prop_quals[i])
				}
			}
		}
	}

	return(prop_quals[length(prop_quals)])
}

#' Plot -log10 Transformed p-values
#'
#' Plots p-values after -log10 transforming them. This makes a plot where points higher up on the y-axis are lower p-values. Points are colored by detection reagant and shaped by antigen being detected.
#'
#' @param p_vals A vector of p-values.
#' @param pdf_name The name of the output file where the plot will be saved. Should be a full file path.
#' @param main The title of the plot.
#' @param line_method A string indicating the method you want to use to draw the arbitrary cutoff line. "bonferroni" (for p < 0.05), "fdr" (for q < 0.2), and "raw" (just -log10(0.05)) are supported.
#' @param gopts A list of graphical options. Must have specified the subparameters:
#' \itemize{
#'	\item gopts$ppar_mar The margins for this plot. See ?par for more info. Default is c(5.1, 4.1, 4.1, 6).
#'	\item gopts$reag_legend_inset The inset for the legend relative to the plot.
#'	\item gopts$reag_cats The categories for your detection reagants.
#'	\item gopts$reag_cols The colors for the detection reagant categories.
#'	\item gopts$ant_shapes The shapes for the antigens being detected.
#' }
#' @param pt_colors Colors for each point in p_vals. These should be generated from gopts$reag_cols. See ?fcgo_get_prop for an example of how to generate this.
#' @param pt_shapes Shapes for each point in p_vals. These should be generated from gopts$ant_shapes. See ?fcgo_get_prop for an example of how to generate this.
#' @param ylim An optional parameter for setting the y-limits of the plot.
#'
#' @return Nothing!
#'
#' @export
fcgo_l10_pplot = function(p_vals, pdf_name, main, line_method, gopts, pt_colors, pt_shapes, ylim=NULL){
	pdf(pdf_name)
	par(mar=gopts$ppar_mar)

	p_vals = -log10(p_vals)

	plot(p_vals, ylim=ylim, xlab="Feature index", ylab="-log10(p-value)", pch=pt_shapes, col=pt_colors, bg=pt_colors, main=main)

	if (line_method == "bonferroni"){
		abline(h=-log10(0.05/length(p_vals)), col="red", lty=2)

	} else if (line_method == "fdr"){
		abline(h=-log10(0.2), col="red", lty=2)

	} else if (line_method == "raw"){
		abline(h=-log10(0.05), col="red", lty=2)
	}

	fch_make_reag_legend(gopts)

	dev.off()
}

#' Make A Volcano Plot
#'
#' Plots p-values versus fold change values. Points are colored by detection reagant and shaped by antigen being detected.
#'
#' @param foch A vector of fold change values.
#' @param p_vals A vector of p-values.
#' @param line_method A string indicating the method you want to use to draw the arbitrary cutoff line. "bonferroni" (for p < 0.05), "fdr" (for q < 0.2), and "raw" (just -log10(0.05)) are supported.
#' @param main The title of the plot.
#' @param xlab The label for the x-axis.
#' @param ylab The label for the y-axis.
#' @param gopts A list of graphical options. Must have specified the subparameters:
#' \itemize{
#'	\item gopts$volc_xlim The x-limits for the volcano plot.
#'	\item gopts$volc_ylim The y-limits for the volcano plot.
#'	\item gopts$reag_legend_inset The inset for the legend relative to the plot. Reasonable place to start is c(-0.18, 0).
#'	\item gopts$reag_cats The categories for your detection reagants.
#'	\item gopts$reag_cols The colors for the detection reagant categories.
#'	\item gopts$ant_shapes The shapes for the antigens being detected.
#' }
#' @param pt_colors Colors for each point in p_vals. These should be generated from gopts$reag_cols. See ?fcgo_get_prop for an example of how to generate this.
#' @param pt_shapes Shapes for each point in p_vals. These should be generated from gopts$ant_shapes. See ?fcgo_get_prop for an example of how to generate this.
#'
#' @return Nothing!
#'
#' @export
fcgo_volcano_plot = function(foch, p_vals, line_method, main, xlab, ylab, gopts, pt_colors, pt_shapes){
		
	# Put them all together
	merged = data.frame(foch, p_vals)

	plot(merged, xlim=gopts$volc_xlim, ylim=gopts$volc_ylim, xlab=xlab, ylab=ylab, main=main, 
	     pch=pt_shapes, col=pt_colors, bg=pt_colors)

	fch_make_reag_legend(gopts)

	abline(v=0, col="black", lty=2)
	abline(v=-1, col="black", lty=2)
	abline(v=1, col="black", lty=2)

	# p-value cutoff - draw a new line based on the method used
	if (line_method == "bonferroni"){
		abline(h=-log10(0.05/length(p_vals)), col="red", lty=2)

	} else if (line_method == "raw"){
		abline(h=-log10(0.05), col="red", lty=2)
	}
}

#' Plot Reduced Fc Array Components
#'
#' Plot the results of fcdr_dimred.
#'
#' @param fc The Fc Array data frame that resulted from fcdr_dimred.
#' @param dims The number of dimensions used for fcdr_dimred.
#' @param pdf_name The name of the output file where the plot will be saved. Should be a full file path.
#' @param main The title of the plot.
#' @param legend The text that will fill in the legend. Should be a vector of strings.
#' @param group_cols The colors for the legend text. Should be a vector of strings of the same length as legend.
#' @param cols The colors for each of the points in the plot. Should be of the same length as the number of points plotted (i.e. nrow(fc)). Should probably be the same colors as group_cols if you want your plot to be interpretable.
#' @param shapes The shapes for each of the points in the plot. Should be of the same length as the number of points plotted (i.e. nrow(fc)).
#'
#' @return Nothing!
#'
#' @examples
#' fc_tnse = fcdr_dimred(fc, "tsne", dims=2)
#' group_cols = c("blue", "red")
#' cols = sapply(fc$group, function(x){ if (x == "PLACEBO"){ return("blue") } else { return("red") } })
#' shapes = sapply(fc$group, function(x){ if (x == "PLACEBO"){ return(21) } else { return(22) } })
#' fcgo_plot_dimred(fc_tsne, 2, "results/tsne_plot.pdf", "tSNE Components for PLAC/VACC", c("PLAC", "VACC"), group_cols, cols, shapes)
#'
#' @export
fcgo_plot_dimred = function(fc, dims, pdf_name, main, legend, group_cols, cols, shapes){
	fc = fc[,fccu_first_feat_col(fc):(fccu_first_feat_col(fc)+(dims-1))]
	inset = ifelse(dims == 2, -0.3, -0.25)

	if (!(dims == 2 || dims == 3)){
		warning("No support for plotting dimensions that aren't 2 or 3 in fcgo_plot_dimred!")
		return(NULL)
	}

	pdf(pdf_name)
	par(mar=c(8,3,4,3))
	if (dims == 2){
		plot(fc, main=main, pch=shapes, col=cols, bg=cols)

	} else if (dims == 3){
		library(scatterplot3d)
		scatterplot3d(fc, main=main, pch=shapes, type="p", highlight.3d=TRUE, angle=120,
			      bg=cols, mar=c(7,3,4,3))
	}
	legend("bottom", legend=legend, col=group_cols, pch=c(21,22), inset=inset, xpd=TRUE, horiz=TRUE)

	dev.off()
	par(mar=c(5,4,4,2))
}

# TODO: Documentation
# This function generates volcano plots based on the output of fc_foch_pvals()

# gopts defaults: volc_xlim=c(-2.5,2.5), volc_ylim=c(0,8)), volc_legend_inset=c(-0.17,0)
	# fcs must be from fcmu_fcs_combs(fc, "group")
fcgo_fcs_volcano = function(fcs, master_fochs, master_p_vals, results_dir, adj_method, line_method, gopts){
	if (!fch_fcs_assert_input(fcs)){
		warning("Returning NULL!")
		return(NULL)
	}

	old_pair_name = "BANANAS"
	old_pair_time = times[2]

	for (i in 1:length(fcs)){
		fc = fcs[[i]]

		full_name = strsplit(names(fcs)[i], '_')[[1]]
		pair_name = paste(full_name[1:2], collapse='_')
		pair_time = full_name[3]

		if (pair_name != old_pair_name){
			# If this is because we hit the next comparison
			if (any(names(dev.list()) == "pdf")){
				dev.off(dev.list()[which(names(dev.list()) == "pdf")])
			}

			old_pair_name = pair_name
		}

		# Start a new plot if we are back to the first non-baseline time point
		if (pair_time == times[2]){
			volc_name = paste(results_dir, pair_name, "_volcano.pdf", sep='')
			pdf(volc_name)
			par(mar=gopts$covb_vpar_mar)

			old_pair_time = pair_time
		}

		fochs = -log10(master_fochs[[i]]["fochs",])
		p_vals = -log10(master_p_vals[[i]]["p_vals",])

		main = paste("Volcano Plot for ", pair_name,
			   "\nTime Point", pair_time, "\n",
			   "adjust: ", adj_method, " - line: ", line_method, sep='')

		# Decide the colors and the shapes for the features
		pt_colors = as.vector(sapply(colnames(fc)[fccu_first_feat_col(fc):ncol(fc)],
					     fcgo_get_prop, gopts$reag_cats, gopts$reag_cols))
		pt_shapes = as.vector(sapply(colnames(fc)[fccu_first_feat_col(fc):ncol(fc)],
					     fcgo_get_prop, gopts$feat_ant, gopts$ant_shapes))

		xlab = "-log10(Fold Change)"
		ylab = "-log10(p-value)"
		fcgo_volcano_plot(fochs, p_vals, line_method, main, xlab, ylab, gopts, pt_colors, pt_shapes)
	}

	# If we reached the end, close the pdf
	if (any(names(dev.list()) == "pdf")){
		dev.off(dev.list()[which(names(dev.list()) == "pdf")])
	}
}
