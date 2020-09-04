# -------------------------------------------
# -------------------------------------------
# COMMON UTILITIES
# -------------------------------------------
# -------------------------------------------

#' Black Magic Add/Subtract for Lists
#'
#' This is probably dangerous. Combines the objects in my_list according to my_op_str. See example.
#'
#' @param my_op_str A string representation of the arithmetic operator you want to use, "+" or "-".
#' @param my_list The list you want to perform this sorcery on.
#' @param max_idx The highest index you want to go to. Defaults to length(my_list).
#'
#' @return The result of the mungo addition or subtraction.
#'
#' @examples
#' > data1 = data.frame(a=c(1,2,3), b=c(4,5,6))
#' > data2 = data.frame(a=c(9,8,7), b=c(6,5,4))
#' > data3 = data.frame(a=c(1,1,1), b=c(2,2,2))
#'
#' > fccu_as_evil("+", list(data1, data2, data3))
#'    a  b
#' 1 11 12
#' 2 11 12
#' 3 11 12
#'
#' # evaluates to the equivalent of writing out:
#' > data1 + data2 + data3
#'
#' @export
fccu_as_evil = function(my_op_str, my_list, max_idx=length(my_list)){
	do_not_do_this = paste(paste("my_list[[", 1:max_idx, "]]", sep=''),
			       collapse=my_op_str)

	return(eval(parse(text=do_not_do_this)))
}

#' Extend a Vector
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

#' Find the First Column that is a Feature
#'
#' Find the first column that is a feature.
#'
#' @return The first column that is not "subject", "group", or "times".
#'
#' @export
fccu_first_feat_col = function(fc){
	return(which(colnames(fc) != "group" & colnames(fc) != "times" & colnames(fc) != "subject")[1])
}

#' Match the Length of a Data Frame and Vector
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

#' Black Magic mapply for Nested Lists
#'
#' This is probably dangerous. It generates and evaluates an mapply call for a nested list, which is really useful if you don't know how long your list will be before runtime.
#'
#' @param func_str The name of the function you want to evaluate in the mapply call.
#' @param my_list The list you want to perform this sorcery on.
#' @param max_idx The highest index you want to go to. Defaults to length(my_list).
#' @param other_args_str A string of the other arguments you want sent to mapply, formatted like "arg1, arg2, arg3". You need the commas. See examples. Defaults to NULL.
#'
#' @return The result of the mapply call.
#'
#' @examples
#' list1 = list(a=data.frame(a1=rep(1.1, 10), a2=rep(1.2, 10)),
#'		b=data.frame(b1=rep(1.1, 10), b2=rep(1.2, 10)))
#' list2 = list(a=data.frame(a1=rep(2.1, 10), a2=rep(2.2, 10)),
#'		b=data.frame(b1=rep(2.1, 10), b2=rep(2.2, 10)))
#' list3 = list(a=data.frame(a1=rep(3.1, 10), a2=rep(3.2, 10)),
#'		b=data.frame(b1=rep(3.1, 10), b2=rep(3.2, 10)))
#' 
#' big_list = list(list1, list2, list3)
#' 
#' # USE.NAMES makes no difference here, included to demonstrate
#' #	other_args_str syntax
#' combined = fccu_mapply_evil("rbind", big_list,
#'		other_args_str="SIMPLIFY=FALSE, USE.NAMES=FALSE")
#' 
#' # combined produces the same result as having written:
#' mapply(rbind, big_list[[1]], big_list[[2]],
#'	big_list[[3]], SIMPLIFY=FALSE, USE.NAMES=FALSE)
#' 
#' # or having written this:
#' mapply(rbind, list1, list2, list3, SIMPLIFY=FALSE,
#'	USE.NAMES=FALSE)
#' 
#' # the result of each is a list of 2 data.frames:
#' #	[[1]] is all the a's rbind-ed together
#' #	[[2]] is all the b's rbind-ed together
#'
#' @export
fccu_mapply_evil = function(func_str, my_list, max_idx=length(my_list), other_args_str=NULL){
	if (!is.null(other_args_str)){
		other_args_str = paste(', ', other_args_str, sep='')
	}

	do_not_do_this = paste("mapply(", func_str, ", ",
			       paste(paste("my_list[[", 1:max_idx, "]]", sep=''),
				     collapse=', '),
			       other_args_str, ')', sep='')

	return(eval(parse(text=do_not_do_this)))
}

#' Get Levels for a Non-Factor
#'
#' This is a shortcut for levels(as.factor(vec)).
#'
#' @param vec The vector you want the "levels" for.
#'
#' @return The result of levels(as.factor(vec)).
#'
#' @export
fccu_plevels = function(vec){
	return(levels(as.factor(vec)))
}

#' Atomic Vector and Matrix Addition with NA Tolerance
#'
#' On a list of matrices/data.frames or on a vector, perform addition while allowing for some NAs to be ignored. Useful if you want to address questions of prevalence where NA is a meaningful result that should not be propagated, but substituting it with 0, 1, etc. doesn't make sense.
#'
#' As a warning, this currently uses R for loops, so it's not fast.
#'
#' @param r_obj An atomic vector acceptable by sum, or a list of matrices or data.frames that have entries acceptable by sum. The latter case also requires that all list elements have the same dimensions.
#' @param ratio How many NAs to tolerate as a ratio of vector length. In the case of matrix/data.frame addition, this is the vector in the third direction, i.e. the Z when the matrices/data.frames are stacked on top of one another.
#'
#' @return For an atomic vector, NA if NAs make up more than the ratio of the vector (i.e. >50 percent of the vector), otherwise the sum of the non-NA elements. For a list of matrices/data.frames, a matrix/data.frame equivalent to manually adding them together with, subject to the behavior for atomic vectors.
#'
#' @examples
#' > m1 = matrix(1, 3, 3)
#' > m2 = matrix(2, 3, 3)
#' > m0 = matrix(NA, 3, 3)
#'
#' > fccu_add(list(m1,m2,m0), 0.5)
#'      [,1] [,2] [,3]
#' [1,]    3    3    3
#' [2,]    3    3    3
#' [3,]    3    3    3
#'
#' > fccu_add(list(m1,m0,m0), 0.5)
#'      [,1] [,2] [,3]
#' [1,]   NA   NA   NA
#' [2,]   NA   NA   NA
#' [3,]   NA   NA   NA
#'
#' > m0[1,] = 3
#' > m0[,3] = 5
#' > fccu_add(list(m1,m0,m0), 0.5)
#'      [,1] [,2] [,3]
#' [1,]    7    7   11
#' [2,]   NA   NA   11
#' [3,]   NA   NA   11
#'
#' @export
fccu_add = function(r_obj, ratio=0.5){
	.fccu_internal_sum = function(vec, ratio){
		if ((sum(is.na(vec)) / length(vec)) > ratio){
			return(NA)
		} 

		return(sum(vec, na.rm=TRUE))
	}

	# If list of matrices, can add in place
	if ("list" %in% class(r_obj)){
		classes = sapply(r_obj, function(x){ return(class(x)[1]) })

		if (all(classes == "matrix") || all(classes == "data.frame")){
			# check dimensions
			dims = lapply(r_obj, dim)

			if (length(unique(dims)) > 1){
				warning("Not all list elements have the same dimensions. Returning NULL,", call.=TRUE)
				return(NULL)
			}

			# add all items on top together
			# for each position, add across all list elements
			nrows = nrow(r_obj[[1]])
			ncols = ncol(r_obj[[1]])

			# TODO: no for loops?
			result = matrix(0, nrows, ncols)
			for (i in 1:nrows){
				for (j in 1:ncols){
					# Grab all the i,j for the list of matrices
					i_j = sapply(r_obj, function(x, i, j){
						return(x[i,j])
					}, i, j)
					
					result[i,j] = .fccu_internal_sum(i_j, ratio)
				}
			}

			# TODO:
			colnames(result) = colnames(r_obj[[1]])
			return(result)

		} else {
			warning("Not all elements of list are matrices or data.frames. Returning NULL.", call.=TRUE)
			return(NULL)
		}

	# If is vector - will throw argument if wrong type so no need to handle
	#	it myself!
	} else if (any(c("logical", "integer", "numeric", "complex", "double", "raw", "character") %in% class(r_obj))){
		return(.fccu_internal_sum(r_obj, ratio))
	}

	warning("Please pass an atomic vector or a matrix. Returning NULL.", call.=TRUE)

	return(NULL)
}
