# -------------------------------------------
# -------------------------------------------
# DIMENSION REDUCTION
# -------------------------------------------
# -------------------------------------------

#' Calculate A Breadth Score
#'
#' Replace a subset of features with a single score representing the geometric mean for each subject across those features. Do NOT pass standardized data.
#'
#' @param fc The Fc Array data frame.
#' @param to_condense A vector of strings of features that you want to collapse into a single score.
#' @param name An optional name parameter to specify for the new score feature. Otherwise, the members of to_condense will be presumed to be of a similar class, and the first substring of the first feature (when split by [.]) will be prepended to "_score".
#'
#' @return The Fc Array data frame, with the features provided replaced with a single score.
#'
#' @examples
#' # pass only the colnames that contain IgG3
#' fc_new = fcdr_breath_score(fc, colnames(fc)[grepl("IgG3", colnames(fc), fixed=TRUE)])
#' # contains the new score feature
#' fc_new[,"IgG3_score"]
#'
#' @export
fcdr_breadth_score = function(fc, to_condense, name=NULL){
	# No categorical variables
	local = output = data.frame(fc[,fccu_first_feat_col(fc):ncol(fc)])
	local = local[,to_condense] # keep
	output = fcff_filter_features(output, to_condense, "toss")

	if (ncol(local) < 1){
		scores = (rep(NA, nrow(local)))

	} else {
		scores = as.vector(apply(local, 1, fch_geom_mean))
	}

	# Replace to_condense features with combined score
	if (fccu_first_feat_col(fc) != 1){
		output = cbind(as.data.frame(fc[,1:(fccu_first_feat_col(fc)-1)]), output, as.data.frame(scores))
		colnames(output)[1:(fccu_first_feat_col(fc)-1)] = colnames(fc)[1:(fccu_first_feat_col(fc)-1)]

	} else {
		output = cbind(output, as.data.frame(scores))
	}

	if (is.null(name)){
		# We assume these are all related features
		new_title = strsplit(to_condense[1], '[.]')[[1]][1]
		new_title = paste(new_title, length(to_condense), "score", sep='_')

	} else {
		new_title = name
	}

	colnames(output)[ncol(output)] = new_title

	return(output)
}

# TODO: Handle subclasses

#' Calculate Many Breadth Scores
#'
#' Replace each subset of features with a single score representing the geometric mean for each subject across those features. Do NOT pass standardized data. Does not work with ordering subclasses (IgG1 before IgG2, etc. - see how fcgo_get_prop works).
#'
#' @param fc The Fc Array data frame.
#' @param to_condense A vector of strings of features that you want to collapse into individual scores.
#'
#' @return The Fc Array data frame, with each class of features provided replaced with individual scores. None of the other feature data remains.
#'
#' @examples
#' # pass only the colnames that contain IgG3
#' gopts$feat_ant = c("gp120","gp140","V1.V2.","gag","p17,p24,p51,p55,p66")
#' # contains gp120_score, gp140_score, V1.V2._score, etc.
#' fc_new = fcdr_breath_score(fc, gopts$feat_ant)
#'
#' @export
fcdr_fc_breadth = function(fc, to_condense){
	fc_breadth = fc

	for (f in to_condense){
		# if has commas, split it, and then sub-eval
		split_feat = strsplit(f, ',')[[1]]
		if (length(split_feat) > 1){
			for (sub_feat in split_feat){
				to_condense = colnames(fc_breadth)[grepl(sub_feat, colnames(fc_breadth), fixed=TRUE)]
				fc_breadth = fcdr_breadth_score(fc_breadth, to_condense, name=paste(sub_feat, length(to_condense), "score", sep='_'))
			}

		} else {
			to_condense = colnames(fc_breadth)[grepl(f, colnames(fc_breadth), fixed=TRUE)]
			fc_breadth = fcdr_breadth_score(fc_breadth, to_condense, name=paste(f, length(to_condense), "score", sep='_'))
		}
	}

	if (fccu_first_feat_col(fc) != 1){
		fc_breadth = cbind(as.data.frame(fc_breadth[,1:(fccu_first_feat_col(fc_breadth)-1)]),
				       as.data.frame(fc_breadth[,grepl("score", colnames(fc_breadth), fixed=TRUE)]))

		colnames(fc_breadth)[1:(fccu_first_feat_col(fc_breadth)-1)] =
			colnames(fc)[1:(fccu_first_feat_col(fc)-1)]

	} else {
		fc_breadth = fc_breadth[,grepl("score", colnames(fc_breadth), fixed=TRUE)]
	}

	toss = vector("numeric")
	for (i in fccu_first_feat_col(fc_breadth):ncol(fc_breadth)){
		if (all(is.na(fc_breadth[,i]))){
			toss = c(toss, i)
		}

		# Inf cannot be handled, but a large number can be (Inf essentially means outlier city)
		measures = fc_breadth[,i][!is.infinite(fc_breadth[,i])]
		fc_breadth[,i][is.infinite(fc_breadth[,i])] = max(measures) + 6*sd(measures)
	}

	if (length(toss) > 0){
		fc_breadth = fc_breadth[,-toss]
	}


	return(fc_breadth)
}

#' Perform Dimension Reduction
#'
#' Function for doing multiple kinds of dimension reduction.
#'
#' @param fc The Fc Array data frame
#' @param method A string "pca", "tsne", or "umap" the denotes what dimension reduction method will be used.
#' @param dims The number of dimensions to reduce the dataset down to. Pass NULL if you want all components back (PCA), or if you feel like letting tSNE and UMAP yell at you that they picked 2 for you.
#' @param perplexity A parameter for tSNE. It will be adjusted for you if you don't provide it/change it and it needs providing/to be changed.
#' @param init_pca A boolean that determines whether or not tSNE performs an initial pca to perform a first-pass on reducing dimensions.
#' @param max_iters The maximum number of iterations you want tSNE or UMAP to perform before terminating.
#' @param theta Speed/accuracy tradeoff for tSNE. 0.0 is exact tSNE, 1.0 is "just get it done ASAP idfc if it's wildly inaccurate".
#' @param seed The random seed you wish to use to facilitate reproducibility.
#' @param verbose A boolean that determines how much you get yelled at by tSNE and UMAP about how they're doing.
#'
#' @return The Fc Array data frame with all the features replaced by all the components generated by the method you used.
#'
#' @examples
#' colnames(fc)
#' # [1] "subject" "group" "IgG3.gp41" "IgG.gp70.V1.V2" # etc.
#' fc_tsne = fcdr_dimred(fc, "tsne", dims=2)
#' colnames(fc_tsne)
#' # [1] "subject" "group" "tsne.1" "tsne.2"
#'
#' @export
fcdr_dimred = function(fc, method, dims=NULL, perplexity=30, init_pca=TRUE, max_iters=500, theta=0.5, seed=1337, verbose=FALSE){
	set.seed(seed)

	# No categorical variables
	output = fc[,fccu_first_feat_col(fc):ncol(fc)]

	if (method == "pca"){
		# Don't center because negative values are dumb
		result = prcomp(output, center=FALSE, scale=TRUE, retx=TRUE)
		new_vals = result$x

		if (!is.null(dims)){
			new_vals = new_vals[,1:dims]
		}

	} else if (method == "tsne"){
		library(Rtsne)
		if (is.null(dims)){
			dims = 2
			print("dims was changed to 2 because you didn't provide a dims")
		}

		if ( is.null(perplexity) || !(3 * perplexity < nrow(output) - 1)){
			print("perplexity was changed because the value you provided was too large for the number of samples")
			print("see fc_tsne for details")
			perplexity = floor((nrow(output) - 1) / 3)
		}

		result = Rtsne(output, dims=dims, perplexity=perplexity, pca=init_pca, verbose=verbose, max_iter=max_iters, theta=theta)
		new_vals = result$Y

	} else if (method == "umap"){
		library(umap)

		if (is.null(dims)){
			dims = 2
			print("dims was changed to 2 because you didn't provide a dims")
		}

		if (umap.defaults$n_neighbors > nrow(output)){
			# These seem reasonable, I guess?
			n_neighbors = sqrt(nrow(output))
			if (nrow(output) >= 6){
				n_neighbors = ceiling(sqrt(nrow(output)))
			}

		} else {
			n_neighbors = umap.defaults$n_neighbors
		}

		result = umap(output, n_components=dims, epochs=max_iters, verbose=verbose,
			      random_state=seed, transform_state=seed, n_neighbors=n_neighbors)
		new_vals = result$layout
	}

	# Replace all features with components
	if (fccu_first_feat_col(fc) != 1){
		output = cbind(as.data.frame(fc[,1:(fccu_first_feat_col(fc)-1)]), as.data.frame(new_vals))
		colnames(output)[1:(fccu_first_feat_col(fc)-1)] = colnames(fc)[1:(fccu_first_feat_col(fc)-1)]

	} else {
		output = as.data.frame(new_vals)
	}
	colnames(output)[fccu_first_feat_col(output):ncol(output)] = paste(method, c(1:ncol(new_vals)),sep='.')

	return(output)
}
