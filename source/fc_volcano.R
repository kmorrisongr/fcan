# Written by Kyle Morrison for analysis of Fc Array Data
# This function generates volcano plots based on the output of fc_foch_pvals()

# gopts defaults: volc_xlim=c(-2.5,2.5), volc_ylim=c(0,8)), volc_legend_inset=c(-0.17,0)
fc_volcano = function(fcs_names, stats_dir, times, adj_method, line_method, gopts, bin_groups, do_bin=FALSE){
	# Remove "all" - a three way volcano plot is silly
	if (any(fcs_names == "all")){
		all_index = which(fcs_names == "all")
		fcs_names = fcs_names[-all_index]
	}

	if (do_bin){
		fcs_names = c(fcs_names, paste(bin_groups, collapse='_'))
	}

	# --------------------------
	# For each comparison

	for (fc in fcs_names){
		fc_dir = paste(stats_dir, fc, "/", sep='')
		volcano_dir = paste(fc_dir, "volcano/", sep='')
		dir.create(volcano_dir)
		foch_dir = paste(fc_dir, "fold_changes/", sep='')
		dir.create(foch_dir)

		# Save all the plots together
		volc_name = paste(volcano_dir, fc, "_volcano.pdf", sep='')
		pdf(file=volc_name)
		par(mar=gopts$vpar_mar)

		for (point in times){
			# Open that fold change file
			foch_name = paste(foch_dir, point, "_fold_change.csv", sep='')
			foch_table = read.table(file=foch_name, sep=',', stringsAsFactors=FALSE, header=TRUE)
			# We only ever have one column for non-all comparisons
			foch = -log10(foch_table[,1])

			# Pull those p-values from the merged file
			p_name = paste(fc_dir, adj_method, "_p_values.csv", sep='')
			p_vals_table = read.table(file=p_name, sep=',', stringsAsFactors=FALSE, header=TRUE)
			keep = sapply(colnames(p_vals_table), function(x,point) { grepl(point,x) }, point=point)
			p_vals = p_vals_table[,keep]
			p_vals = -log10(p_vals)

			name = strsplit(fc,"_")[[1]]
			first_name = name[1]
			second_name = name[2]
			main = paste("Volcano Plot for ",toupper(first_name)," vs ",toupper(second_name),
				   "\nTime Point",point,'\n',
				   "adjust: ",adj_method," - line: ",line_method,sep='')

			# Decide the colors and the shapes for the features
			pt_colors = as.vector(sapply(rownames(foch_table), fcu_get_color, gopts$reag_cats, gopts$reag_cols))
			pt_shapes = as.vector(sapply(rownames(foch_table), fcu_get_shape, gopts$feat_ant, gopts$ant_shapes))

			xlab = "-log10(Fold Change)"
			ylab = "-log10(p-value)"
			fcu_volcano_plot(foch, p_vals, line_method, feat_names=rownames(foch_table), main, 
				     xlab, ylab, gopts, pt_colors, pt_shapes)
		}

		dev.off()
	}
}
