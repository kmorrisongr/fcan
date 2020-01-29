# Compilation of utility functions for working with Fc Arrays


# -------------------------------------------
# -------------------------------------------
# SECTIONS
# -------------------------------------------
# -------------------------------------------
# 1. COMMON UTILITES
# 2. DATASET CLEANUP
# 3. FEATURE FILTERING
# 4. DIMENSION REDUCTION
# 5. STATISTICAL TESTING
# 6. STATISTICAL MODELING
# 7. GRAPHICAL OUTPUT
# 8. MULTI-FC UTILITIES
# 9. HELPER FUNCTIONS


# -------------------------------------------
# -------------------------------------------
# 1. COMMON UTILITIES
# -------------------------------------------
# -------------------------------------------

# Extend vector to match max_len by appending to_app
fcu_extend = function(vec, max_len, to_app=NA){
	output = vec

	if (length(output) < max_len){
		output = c(output, rep(to_app,max_len-length(output)))
	}

	return(output)
}

# Which column has the first feature (i.e. not times or group or subject)
fcu_first_feat_col = function(fc){
	return(which(colnames(fc) != "group" & colnames(fc) != "times" & colnames(fc) != "subject")[1])
}

# Remove entries from data and a until a is the same length as b
# By checking if a[i] == b[i] until you find one that isn't
# a must be from data as well (say, rownames())
fcu_match_rows = function(a, b, data){
	local_a = a; local_b = b; local_data = data;
	while (length(local_a) > length(local_b)){
		toss = 0
		for (i in 1:length(local_a)){
			if ( is.na(local_a[i] != local_b[i]) ||
			    (local_a[i] != local_b[i]) ){
				toss = i
				break
			}
		}

		local_a = local_a[-toss]
		local_data = local_data[-toss,]
	}

	return(local_data)
}


# -------------------------------------------
# -------------------------------------------
# 2. DATASET CLEANUP
# -------------------------------------------
# -------------------------------------------

# Impute values in the Fc Array based on medians within group measurements
# Make sure you have group information in the Fc Array
# Provide the column number where features start
fcu_impute = function(fc, start_col, type="det", method=median){
	local = fc
	groups = levels(local[,"group"])

	for (i in start_col:ncol(local)){
		for (g in groups){
			current = (local[,"group"] == g)
			measures = local[current,i]

			measures[is.na(measures)] = method(measures, na.rm=TRUE)

			local[current,i] = measures
		}
	}

	return(local)
}

# Do background subtraction using base_ids = (local[,"group"] == "PLACEBO") to identify the baseline group
fcu_sub_back = function(fc, base_ids){
	local = fc

	for (i in fcu_first_feat_col(local):ncol(local)){
		# Background subtract PLACEBO from VACCINE
		placebo = local[base_ids,i]
		vaccine = local[!base_ids,i]

		# Anyone who is negative gets set to 0 (because negative makes no sense)
		placebo[placebo < 0] = 0 
		background = median(placebo, na.rm=TRUE)
		back_sd = sd(placebo, na.rm=TRUE) 

		vaccine = vaccine - background
		vaccine[vaccine < 0] = 0

		# Everyone who is PLACEBO for a feature gets set to 0 
		local[base_ids,i] = 0
		
		# Scale everyone according to the variation in the measurement itself (i.e. in PLACEBO)
		vaccine = vaccine/back_sd 
		local[!base_ids,i] = vaccine
	}

	return(local)
}


# -------------------------------------------
# -------------------------------------------
# 3. FEATURE FILTERING
# -------------------------------------------
# -------------------------------------------

fcu_exclude_attr = function(fc, attr, clm){
	toss = grepl(attr, fc[,clm])

	return(fc[!toss,])
}

# Do feature filtering using base_ids = (local[,"group"] == "PLACEBO") to keep only features that differ from baseline
#diff_opts = list()
#diff_opts$base_ids = base_ids
#diff_opts$test = t.test
#diff_opts$alternative = "two.sided"
#diff_opts$adj_method = "none"
fcu_feats_diff = function(local, diff_opts){
	p_vals = as.vector(apply(local, 2, function(x, base_ids, test, alternative){
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

# Returns a data frame that has only the features you want
# Pass a vector of the features you want to keep or toss c("gp41","gp140")
# If you want to filter by difference from baseline, pass a diff_opts (see fcu_feats_diff)
fcu_filter_features = function(fc, to_keep, action, behavior="permissive", diff_opts=NULL){
	local = data.frame(fc[,fcu_first_feat_col(fc):ncol(fc)])
	feats = colnames(local) = colnames(fc)[fcu_first_feat_col(fc):ncol(fc)]
	subjects = fc$subject

	keep = vector("logical")

	if ( (behavior == "differs") && (!is.null(diff_opts)) ){
		p_vals = feats_diff(local,diff_opts)
		keep = p_vals < 0.05

	} else {
		# Written this way to accomdate & and |
		for (i in 1:length(to_keep)){
			i_keep = grepl(to_keep[i],feats)

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

	indices = ifelse(rep(action == "toss",length(keep)), !keep, keep)
	output = as.data.frame(local[,indices])
	colnames(output) = colnames(local)[indices]

	if (fcu_first_feat_col(fc) != 1){
		output = cbind(as.data.frame(fc[,1:(fcu_first_feat_col(fc)-1)]), output)
		colnames(output)[1:(fcu_first_feat_col(fc)-1)] = colnames(fc)[1:(fcu_first_feat_col(fc)-1)]
	}

	if (!is.null(output$group)){
		output$group = factor(output$group)
	}

	return(output)
}

fcu_only_attr = function(fc, attr, attrs, clm){
	output = fc
	for (a in attrs){
		if (a != attr){
			output = fcu_exclude_attr(output, a, clm)
		}
	}

	return(output)
}

# Remove correlated features from data based on cor_cutoff
fcu_remove_cor = function(data, cor_cutoff, keep_behavior){
	cor_mat = cor(data)
	high_cor = findCorrelation(cor_mat, cutoff=cor_cutoff)

	# Randomly keep one of the highly correlated features - they should,  by definition,  be interchangeable
	if (keep_behavior == "sample"){
		keep = sample(high_cor, 1)
		keep = which(high_cor == keep)

	} else if (keep_behavior == "first"){
		keep = high_cor[1]
	}
	high_cor = high_cor[-keep]

	data = data[,-high_cor]

	return(data)
}

# wrapper for feature filtering in the master file
fcu_wrap_filter = function(fc, results_dir, keep_filter, discard_filter, k_behavior, d_behavior, diff_opts="NULL"){
	local = fc

	cat("Fc Array has", ncol(local)-fcu_first_feat_col(local)-1, "features", '\n')
	if (k_behavior == "differs" | d_behavior == "differs"){
		local = fcu_filter_features(local, keep_filter, "keep", "differs", diff_opts)
		results_dir = paste(results_dir, "differs", '/', sep='')

	} else if (length(keep_filter) > 0 & length(discard_filter) > 0){
		stamp = toString(Sys.time())
		stamp = paste(strsplit(stamp,' ')[[1]], collapse='_')
		stamp = gsub('-','',stamp)
		stamp = gsub(':','',stamp)

		results_dir = paste(results_dir, "complex_", stamp, '/', sep='')
		dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

		local = fcu_filter_features(local, keep_filter, "keep", k_behavior)
		local = fcu_filter_features(local, discard_filter, "toss", d_behavior)

		filter_file = paste(results_dir, "feat_filter.txt", sep='')
		write(paste("only_", keep_filter, sep=''), filter_file)
		write(paste("no_", discard_filter, sep=''), filter_file, append=TRUE)

	} else if (length(keep_filter) > 0){
		local = fcu_filter_features(local, keep_filter, "keep", k_behavior)
		results_dir = paste(results_dir, "only_", paste(keep_filter, collapse="-"), '/', sep='')

	} else if (length(discard_filter) > 0){
		local = fcu_filter_features(local, discard_filter, "toss", d_behavior)
		results_dir = paste(results_dir, "no_", paste(discard_filter, collapse="-"), '/', sep='')
	}

	dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

	cat("Filtered Fc Array has", ncol(local)-fcu_first_feat_col(local)-1, "features", '\n')
	return(list(fc=local, results_dir=results_dir))
}


# -------------------------------------------
# -------------------------------------------
# 4. DIMENSION REDUCTION
# -------------------------------------------
# -------------------------------------------

# Replace all fc[,to_condense] with a single score representing the geometric mean for each subject across those features
fcu_breadth_score = function(fc, to_condense){
	# No categorical variables
	local = output = fc[,fcu_first_feat_col(fc):ncol(fc)]
	local = local[,to_condense]

	local = data.frame(fc[,fcu_first_feat_col(fc):ncol(fc)])
	local = local[,to_condense]
	output = fcu_filter_features(output,to_condense,"toss")

	scores = as.vector(apply(local,1,geom_mean))

	# Replace to_condense features with combined score
	if (fcu_first_feat_col(fc) != 1){
		output = cbind(as.data.frame(fc[,1:(fcu_first_feat_col(fc)-1)]), output, as.data.frame(scores))
		colnames(output)[1:(fcu_first_feat_col(fc)-1)] = colnames(fc)[1:(fcu_first_feat_col(fc)-1)]

	} else {
		output = cbind(output, as.data.frame(scores))
	}

	# We assume these are all related features
	new_title = strsplit(to_condense[1], '[.]')[[1]][1]
	new_title = paste(new_title, "_score", sep='')

	colnames(output)[ncol(output)] = new_title

	return(output)
}

# Function for doing multiple kinds of dimensionality reduction
# Here, unlike elsewhere, method is a string flag, not a function pointer
fcu_dimred = function(fc, method="pca", dims=NULL, perplexity=30, init_pca=TRUE, max_iters=500, theta=0.5, seed=1337, verbose=FALSE){
	set.seed(seed)

	# No categorical variables
	local = fc[,fcu_first_feat_col(fc):ncol(fc)]

	if (method == "pca"){
		# Don't center because negative values are dumb
		result = prcomp(local, center=FALSE, scale=TRUE, retx=TRUE)
		new_vals = result$x

	} else if (method == "tsne"){
		library(Rtsne)
		if (is.null(dims)){
			dims = length(gopts$reag_cats + gopts$feat_ant)
		}

		if (!(3 * perplexity < nrow(local) - 1)){
			print("perplexity was changed because the value you provided was too large for the number of samples")
			print("see fc_tsne for details")
			perplexity = floor((nrow(local) - 1) / 3)
		}

		result = Rtsne(local, dims=dims, perplexity=perplexity, pca=init_pca, verbose=verbose, max_iter=max_iters, theta=theta)
		new_vals = result$Y

	} else if (method == "umap"){
		library(umap)

		custom_config = umap.defaults
		custom_config$n_components = dims
		custom_config$epochs = max_iters
		custom_config$verbose = verbose
		custom_config$random_state = seed
		custom_config$transform_state = seed
		if (custom_config$n_neighbors > nrow(local)){
			# These seem reasonable, I guess?
			custom_config$n_neighbors = sqrt(nrow(local))
			if (nrow(local) >= 6){
				custom_config$n_neighbors = ceiling(sqrt(nrow(local)))
			}
		}

		result = umap(local, config=custom_config)
		new_vals = result$layout
	}

	# Replace all features with components
	if (fcu_first_feat_col(fc) != 1){
		local = cbind(as.data.frame(fc[,1:(fcu_first_feat_col(fc)-1)]), as.data.frame(new_vals))
		colnames(local)[1:(fcu_first_feat_col(fc)-1)] = colnames(fc)[1:(fcu_first_feat_col(fc)-1)]

	} else {
		local = as.data.frame(new_vals)
	}
	colnames(local)[fcu_first_feat_col(local):ncol(local)] = paste(method, c(1:ncol(new_vals)),sep='.')

	return(local)
}

fcu_plot_dimred = function(fc, dims, pdf_name, main, legend, group_cols, cols, shapes){
	pdf(pdf_name)
	par(mar=c(8,3,4,3))
	if (dims == 2){
		plot(fc[,fcu_first_feat_col(fc):(fcu_first_feat_col(fc)+(dims-1))], main=main, pch=shapes, col=cols, bg=cols)
		legend("bottom", legend=legend, col=group_cols, pch=c(21,22), inset=-0.3, xpd=TRUE, horiz=TRUE)

	} else if (dims == 3){
		library(scatterplot3d)
		scatterplot3d(fc[,fcu_first_feat_col(fc):(fcu_first_feat_col(fc)+(dims-1))], main=main, pch=shapes, type="p",
			      highlight.3d=TRUE, angle=120, bg=cols, mar=c(7,3,4,3))

		legend("bottom", legend=legend, col=group_cols, pch=c(21,22), inset=-0.25, xpd=TRUE, horiz=TRUE)
	}
	dev.off()
	par(mar=c(5,4,4,2))
}


# -------------------------------------------
# -------------------------------------------
# 5. STATISTICAL TESTING
# -------------------------------------------
# -------------------------------------------

# Calculate change over baseline for each subject using method
# Returns a longitudinal dataframe with values replaced by covb
# 	Obviously, the first time point is absent from output
fcu_covb = function(fc, method=fch_covb_simple){
	times = unique(fc$times)
	subjects = unique(fc$subject)

	baseline = fcu_only_attr(fc, times[1], times, "times")
	output = fcu_exclude_attr(fc, times[1], "times")

	for (s in subjects){
		idx = which(output$subject == s)[1]

		# If subject not present after baseline
		if (!is.na(idx)){
			for (t in times[-1]){
				for (i in fcu_first_feat_col(output):ncol(output)){
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

# Generate fold change between two groups of measurements
# Perform Mann-Whitney U test comparing these two groups of measurements
# Pass local_data - columns are measures and groups that correspond
fcu_sing_comp_foch = function(local_data, comp_name, method=median, alternative="two.sided", exact=FALSE){
	# Remove any NAs from this measurement data
	measure = as.numeric(na.omit(local_data$measure))
	remove = is.na(local_data$measure)
	if (length(remove) > 0){
		local_groups = local_data$groups[!remove]
	} else {
		local_groups = local_data$groups
	}

	# m_k -> [1] "m" [2] "k"
	name = strsplit(comp_name,"_")[[1]]
	first_name = name[1]
	second_name = name[2]

	id_first = which(tolower(as.character(local_groups)) == tolower(first_name))
	id_second = which(tolower(as.character(local_groups)) == tolower(second_name))
	local_first = measure[id_first]
	local_second = measure[id_second]

	sing_foch = method(local_first)/method(local_second)

	# Perform Mann-Whitney U test
	if (alternative == "either"){
		test = wilcox.test(measure ~ local_groups, alternative="two.sided", exact=exact)
		p_val = (test$p.value)/2

	} else {
		test = wilcox.test(measure ~ local_groups, alternative=alternative, exact=exact)
		p_val = test$p.value
	}

	output = list(sing_foch=sing_foch,p_val=p_val)
	return(output)
}


# -------------------------------------------
# -------------------------------------------
# 6. STATISTICAL MODELING
# -------------------------------------------
# -------------------------------------------

fcu_test_model = function(model, test, test_group){
	test_results = predict(model, test)
	test_results = cbind(actual=test_group, test_results)
	print(test_results)

	score = sum(as.vector(test_results$pred) == as.vector(test_results$actual))/nrow(test_results)
	cat("Overall accuracy when predicting on the test data was", score*100, "%\n")

	return(score)
}


# -------------------------------------------
# -------------------------------------------
# 7. GRAPHICAL OUTPUT
# -------------------------------------------
# -------------------------------------------

# A function for getting the color for a reagant
fcu_get_color = function(feat, reag_cats, reag_cols){
	feat = tolower(feat); reag_cats = tolower(reag_cats);

	for (i in 1:length(reag_cats)){
		if (any(grepl(reag_cats[i], feat))){
			return(reag_cols[i])
		}

		# if this reag_cats has commas, split it, and then sub-check
		split_reag = strsplit(reag_cats[i], ',')[[1]]
		if (length(split_reag) > 1){
			# For any of these, return the same reag_cols
			for (sub_reag in split_reag){
				if (any(grepl(sub_reag, feat))){
					return(reag_cols[i])
				}
			}
		}
	}

	return(reag_cols[length(reag_cols)])
}

# A function for getting the shape for an antigen
fcu_get_shape = function(feat, feat_ant, ant_shapes){
	feat = tolower(feat); feat_ant = tolower(feat_ant);

	for (i in 1:length(feat_ant)){
		if (any(grepl(feat_ant[i], feat))){
			return(ant_shapes[i])
		}

		# if this feat_ant has commas, split it, and then sub-check
		split_feat = strsplit(feat_ant[i], ',')[[1]]
		if (length(split_feat) > 1){
			# For any of these, return the same ant_shape
			for (sub_feat in split_feat){
				if (any(grepl(sub_feat, feat))){
					return(ant_shapes[i])
				}
			}
		}
	}

	return(ant_shapes[length(ant_shapes)])
}

# Plot -log10(p_vals)
# mar=c(5.1,4.1,4.1,6)
fcu_l10_pplot = function(pdf_name, p_vals, main, line_method, gopts, pt_colors, pt_shapes, ylim=NULL){
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

	par(xpd=TRUE)
	fcu_make_reag_legend(gopts)
	par(xpd=FALSE)

	dev.off()
}

fcu_make_reag_legend = function(gopts){
	# All the rep calls are so we can draw the reagant colors and the antigen shapes together
	legend("topright", inset=gopts$reag_legend_inset, cex=0.85,
	       legend=c(gopts$reag_cats, "other", gopts$feat_ant_legend, "other"),

	       bty="n", text.col=c(gopts$reag_cols, rep( "black",length(gopts$ant_shapes)+1 )),

	       col=c(gopts$reag_cols, rep( "black",length(gopts$ant_shapes) )),
	       pt.bg=c(gopts$reag_cols, rep( "white",length(gopts$ant_shapes) )),

	       # Just make the shapes for the reagant colors circles
	       pch=c(rep(21,length(gopts$reag_cols)), c(gopts$ant_shapes)))
}

# A function that does some serious volcano plotting, mate
# Defaults: volc_xlim=c(-2.5,2.5), volc_ylim=c(0,8), inset=c(-0.18,0)
fcu_volcano_plot = function(foch, p_vals, line_method, feat_names, main, xlab, ylab, gopts, pt_colors, pt_shapes, outside_legend=TRUE){
		
	# Put them all together
	merged = data.frame(foch, p_vals)

	plot(merged, xlim=gopts$volc_xlim, ylim=gopts$volc_ylim, xlab=xlab, ylab=ylab, main=main, 
	     pch=pt_shapes, col=pt_colors, bg=pt_colors)

	if (outside_legend){
		par(xpd=TRUE)
		fcu_make_reag_legend(gopts)
		par(xpd=FALSE)

	} else {
		gopts$reag_legend_inset = 0
		fcu_make_reag_legend(gopts)
	}

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


# -------------------------------------------
# -------------------------------------------
# 8. MULTI-FC UTILITES
# -------------------------------------------
# -------------------------------------------

# Returns a list of datasets representing pairwise combinations of the groups
fcu_fcs_combs = function(fc, groups){
	fcs = list()

	# For each combination of groups - pairwise because fold change stuff
	# Won't change anything if only two groups
	group_combs = combn(groups,2)
	for (i in 1:ncol(group_combs)){
		# Exclude the other ones
		to_exclude = setdiff(groups, group_combs[,i])

		local = fc
		for (j in to_exclude){
			local = fcu_exclude_attr(local, j, "group")
		}

		fcs[[i]] = local
		names(fcs)[i] = tolower(paste(group_combs[,i], collapse="_"))
	}

	return(fcs)
}

# Iterate over a list of data.frames, find the min and max values for the specified feature column
fcu_fcs_find_ylim = function(data_list, data_clm){
	y_min = NULL
	y_max = NULL

	for (d in data_list){
		measures = na.omit(d[,data_clm])

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
		write("No limits found - check your data, fc_find_ylim(), and potentially fc_covb.R",file='')

	# If y_max but no y_min
	} else if (is.null(y_min)){
		ylim = c(0, y_max)
	}

	return(ylim)
}


# -------------------------------------------
# -------------------------------------------
# 9. HELPER FUNCTIONS
# -------------------------------------------
# -------------------------------------------

# Helper for fc_covb
fch_covb_simple = function(base_measure, measure){
	return(measure/base_measure)
}

fch_geom_mean = function(vec){
	return( prod(vec)**(1/length(vec)) )
}

fch_min_max_feat = function(vec){
	return(as.vector(sapply(vec, z_score, vec=vec)))
}

fch_z_score = function(x, vec){
	return( (x - min(vec))/(max(vec) - min(vec)) )
}
