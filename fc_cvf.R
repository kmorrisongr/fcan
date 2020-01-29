# Written by Kyle Morrison for analysis of Fc Array Data
# This function plots feature values relative to challenge outcome

# gopts$cpar_mfrow = c(2,3) ; gopts$cpar_mar = c(5.1,4.1,4.1,3.0)
# gopts$clegend_inset = c(-0.22,0) ; flags$do_same_range = TRUE
# trans_method = I for no transformations
fc_cvf = function(fc, to_cvf, surv, stats_dir, flags, gopts, trans_method=log10){
	fc = fcu_filter_features(fc,to_cvf,"keep")
	fc = fcu_impute(fc,fcu_first_feat_col(fc))
	feats = sort(to_cvf)

	groups = levels(fc[,"group"])
	times = unique(fc[,"times"])
	num_groups = length(groups)

	if (flags$do_same_range){
		file_name = paste(stats_dir,"chals_vs_feats_srange.pdf",sep='')

	} else {
		file_name = paste(stats_dir,"chals_vs_feats.pdf",sep='')
	}

	pdf(file_name)	
	par(mar=gopts$cpar_mar,mfrow=gopts$cpar_mfrow)
	for (f in feats){
		#feature_title = ruth_pretty_title(f)
		feature_title = f

		if (flags$do_same_range){
			# Find the min and max values
			feat_measures = trans_method(fc[,f])

			y_min = min(na.omit(feat_measures))
			y_max = max(na.omit(feat_measures))
			ylim = c(y_min,y_max)

		} else {
			ylim = NULL
		}

		# Keep all time points for a feature on the same page
		layout_counter = 1

		for (point in times){
			# Get our feature's measurements and challenge data
			fc_time = fcu_only_attr(fc, point, times, "times")
			measures = fc_time[,f]

			measures = trans_method(measures)

			# Which challenges correspond to the subjects at this time point?
			local_chals = vector("numeric")
			for (s in fc_time$subject){
				local_chals = c(local_chals,surv$challenges[s == surv$id])
			}

			# Colors for each subject based on group
			time_subj_cols = vector("character")
			for (s in 1:length(fc_time$subject)){
				time_subj_cols = c(time_subj_cols,gopts$plot_group_cols[which(groups == fc_time$group[s])])
			}

			main = paste(feature_title, "\nTime Point ", point, sep='')
			xlab = "Number of Challenges"
			# deparse()?
			ylab = "log10(Measurement)"

			plot(x=local_chals, y=measures, main=main, pch=21, col=gopts$time_subj_cols, 
			     bg=gopts$time_subj_cols, xlab=xlab, ylab=ylab, ylim=ylim)

			layout_counter = layout_counter + 1

			par(xpd=TRUE)
			legend("topright", inset=legend_inset, cex=0.85, legend=groups, 
			       bty="n",  text.col=plot_group_cols,  col=plot_group_cols,  pt.bg=plot_group_cols, 
			       pch=c(rep(21,length(plot_group_cols))))
			par(xpd=FALSE)

			if (layout_counter > length(times)){
				plot.new()
				layout_counter = 1
			}
		}
	}

	dev.off()
}
