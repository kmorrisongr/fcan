# -------------------------------------------
# -------------------------------------------
# STATISTICAL MODELING
# -------------------------------------------
# -------------------------------------------

#' Test A Statistical Model
#'
#' A wrapper for testing statistical models.
#'
#' @param model The model you generated, be it from lm, glm, etc.
#' @param test The test data frame you'll be using. Should not include any group membership information.
#' @param test_group The actual group membership for your test data.
#'
#' @return The accuracy of your model when predicting on the test data.
#'
#' @export
fcsm_test_model = function(model, test, test_group){
	test_results = predict(model, test)
	test_results = cbind(actual=test_group, test_results)
	print(test_results)

	score = sum(as.vector(test_results$pred) == as.vector(test_results$actual))/nrow(test_results)
	cat("Overall accuracy when predicting on the test data was", score*100, "%\n")

	return(score)
}

# TODO: Documentation
# Written by Kyle Morrison for analysis of Fc Array Data
# This function performs univariate analysis for each feature

#' Perform Univariate Analysis
#'
#' Perform univariate analysis on each feature in the Fc Array data frame. Return both p-values and q-values (based on adj_method).
#'
#' @param fc The Fc Array data frame. Will be standardized for you so you are comparing apples to apples.
#' @param adj_method The method you want to use to adjust p-values. Use "none" as a passthrough.
#'
#' @export
fcsm_uni = function(fc, adj_method){
	fc = fcdc_normalize(fc)

	results = as.data.frame(matrix(nrow=ncol(fc)-fccu_first_feat_col(fc)-1, ncol=3))
	colnames(results) = c("feature", "p_val", "q_val")

	for (i in fccu_first_feat_col(fc):ncol(fc)){
		idx = i - (fccu_first_feat_col(fc) - 1)
		model = glm(fc$group ~ fc[,i], family=binomial(link='logit'))
		output = c(colnames(fc)[i], summary(model)$coefficients[,"Pr(>|z|)"][2], NA)
		results[idx,] = output
	}

	results[,"p_val"] = as.numeric(results[,"p_val"])
	results[,"q_val"] = p.adjust(results[,"p_val"], method=adj_method)

	return(results)
}
