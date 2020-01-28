# Written by Kyle Morrison for analysis of Fc Array Data
# This function generates nice box plots 
	# to_box are the features that you want boxplotted - get these from cons_final,  or generate a vector yourself

fc_boxplot = function(fc, to_box, stats_dir, flags, gopts){
	fc = fcu_filter_features(fc, to_box, "keep")
	fc = fcu_impute(fc, fcu_first_feat_col(fc))
	feats = sort(to_box)

	groups = levels(fc$group)
	times = unique(fc$times)
	num_groups = length(groups)


	# -------------------------------------------
	# Dataframe list broken down by time point
	# -------------------------------------------

	if (!grepl("dummy", times)){
		datasets = list()
		for (i in 1:length(times)){
			datasets[[i]] = fcu_only_attr(fc, times[i], times, "times")
			names(datasets)[i] = times[i]
		}

	} else {
		datasets = list(fc)
	}


	# -------------------------------------------
	# Dataframe list further broken down by group
	# -------------------------------------------

	df_list = list()
	k = 1
	for (i in 1:length(groups)){
		for (j in 1:length(times)){
			df_list[[k]] = fcu_only_attr(datasets[[j]], groups[i], groups, "group")
			k = k + 1
		}
	}

	pdf_name = paste(stats_dir, paste(tolower(groups), collapse='_'), "_boxplots.pdf", sep='')
	pdf(pdf_name)
	par(mar=gopts$bpar_mar, mfrow=c(2, 2))


	# -------------------------------------------
	# For each feature,  make a boxplot
	# -------------------------------------------

	for (f in feats){
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
			to_add = seq(j, (length(df_list)+num_groups+j),  by=(num_groups+1))
			box_at = c(box_at, to_add)
		}
		if (length(times) == 1){
			new_box_at = vector("numeric")
			for (j in 1:num_groups){
				new_box_at = c(new_box_at, box_at[j + ((j-1)*(num_groups-1))])
			}
			box_at = new_box_at
		}


		# -------------------------------------------
		# Plot!

		# Might have to write your own pretty_title,  depending on what you want
		#main = ruth_pretty_title(i)
		main = f

		# ?ifelse to understand [1]
		xlab = ifelse(!grepl("dummy",times), "Time Point", '')[1]
		#ylab="log10(MFI)"
		ylab = "MFI"

		boxplot(measures ~ box_groups, at=box_at, names=NA, col="white", 
			main=main, outline=FALSE, xlab=xlab, ylab=ylab)
		stripchart(measures ~ box_groups, at=box_at, vertical=TRUE, method="jitter", add=TRUE, pch=20, 
			   col=rep(gopts$jplot_group_cols, each=length(times)))

		labels = ifelse(grepl("dummy",times), '', times)
		axis(1, at=seq((0.5*(num_groups+1)), (num_groups+1)*length(times), by=num_groups+1), labels=labels)

		# print(groups) before reassigning to make sure they line properly - no reversed findings!
		if (flags$nicer_legend){
			groups = c("Uninfected", "Infected")
		}
		par(xpd=TRUE)
		legend("topright", inset=gopts$box_legend_inset, legend=groups, fill=gopts$bplot_group_cols, cex=0.7)
		par(xpd=FALSE)
	}

	dev.off()
}
