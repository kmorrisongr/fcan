# Written by Kyle Morrison for analysis of Fc Array Data
# This function performs univariate analysis for each feature

# flags$do_standard = TRUE - You should standardize if you want to compare apples to apples
fc_uni = function(fc, adj_method, flags){
	fc = fcu_impute(fc, fcu_first_feat_col(fc))

	if (flags$do_standard){
		# Normalize each feature
		for (i in fcu_first_feat_col(fc):ncol(fc)){
			fc[,i] = fch_min_max_feat(fc[,i])
		}
	}

	results = as.data.frame(matrix(nrow=ncol(fc)-fcu_first_feat_col(fc)-1,ncol=3))
	colnames(results) = c("feature","p_val","q_val")

	for (i in fcu_first_feat_col(fc):ncol(fc)){
		idx = i - (fcu_first_feat_col(fc) - 1)
		model = glm(fc$group ~ fc[,i],family=binomial(link='logit'))
		output = c(colnames(fc)[i], summary(model)$coefficients[,"Pr(>|z|)"][2], NA)
		results[idx,] = output
	}

	results[,"p_val"] = as.numeric(results[,"p_val"])
	results[,"q_val"] = p.adjust(results[,"p_val"],method=adj_method)

	return(results)
}
