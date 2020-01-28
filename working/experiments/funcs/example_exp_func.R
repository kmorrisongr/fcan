exp_20200120 = function(fc, results_dir, adj_method, flags, gopts){
	# ---------------------------
	# Perform univariate analysis
	uni_df = fc_uni(fc, adj_method, flags)
	file_name = paste(results_dir, "univariate_results.csv", sep='')
	write.table(uni_df, file=file_name, sep=', ', quote=FALSE, row.names=FALSE, col.names=TRUE)

	# ---------------------------
	# Plot p-values and q-values
	pt_colors = as.vector(sapply(uni_df[,1], fcu_get_color, gopts$reag_cats, gopts$reag_cols))
	pt_shapes = as.vector(sapply(uni_df[,1], fcu_get_shape, gopts$feat_ant, gopts$ant_shapes))

	pdf_name = paste(results_dir, "univariate_p_vals.pdf", sep='')
	fcu_l10_pplot(pdf_name, uni_df[,2], main="-log10(p_vals) for univariate predictions", line_method="raw", gopts, pt_colors, pt_shapes)

	pdf_name = paste(results_dir, "univariate_q_vals.pdf", sep='')
	fcu_l10_pplot(pdf_name, uni_df[,3], main="-log10(q_vals) for univariate predictions", line_method="fdr", gopts, pt_colors, pt_shapes, ylim=gopts$lpp_ylim)

	# ---------------------------
	# Significant boxplots I/A
	to_plot = uni_df[,1][which(uni_df[,2] < 0.05)]
	if (length(to_plot) > 0){
		fc_boxplot(fc, to_box=to_plot, results_dir, flags, gopts)

	} else {
		print("No features were significant")
	}
}
