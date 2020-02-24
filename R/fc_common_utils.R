# -------------------------------------------
# -------------------------------------------
# COMMON UTILITIES
# -------------------------------------------
# -------------------------------------------

#' Extend A Vector
#'
#' Extend vector to match max_len by appending to_app.
#'
#' @param vec A vector.
#' @param max_len How long you want vec to be.
#' @param to_app What you want to be appended to vec.
#'
#' @return A vector of length max_len.
#'
#' @export
fccu_extend = function(vec, max_len, to_app=NA){
	if (length(vec) < max_len){
		vec = c(vec, rep(to_app, max_len-length(vec)))
	}

	return(vec)
}

# TODO: Less hard-coded (Chris' METACOLS)

#' Find the First Column That Is A Feature
#'
#' Find the first column that is a feature.
#'
#' @return The first column that is not "subject", "group", or "times".
#'
#' @export
fccu_first_feat_col = function(fc){
	return(which(colnames(fc) != "group" & colnames(fc) != "times" & colnames(fc) != "subject")[1])
}

#' Match The Length Of A Data Frame And Vector
#'
#' Make data and vector same length/nrow as other vector. Remove entries from data and a until a is the same length as b. By checking if a[i] == b[i] until you find one that isn't.
#'
#' @param a A vector to be trimmed, must be from data as well (say, rownames()).
#' @param b The vector of shorter length that you want length(a) and nrow(data) to be.
#' @param data Some data frame you want trimmed.
#'
#' @return A trimmed data frame.
#'
#' @export
fccu_match_rows = function(a, b, data){
	while (length(a) > length(b)){
		toss = 0
		for (i in 1:length(a)){
			if ( is.na(a[i] != b[i]) || (a[i] != b[i]) ){
				toss = i
				break
			}
		}

		a = a[-toss]
		data = data[-toss,]
	}

	return(data)
}
