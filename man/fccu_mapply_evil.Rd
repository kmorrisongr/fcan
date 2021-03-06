% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fc_common_utils.R
\name{fccu_mapply_evil}
\alias{fccu_mapply_evil}
\title{Black Magic mapply for Nested Lists}
\usage{
fccu_mapply_evil(
  func_str,
  my_list,
  max_idx = length(my_list),
  other_args_str = NULL
)
}
\arguments{
\item{func_str}{The name of the function you want to evaluate in the mapply call.}

\item{my_list}{The list you want to perform this sorcery on.}

\item{max_idx}{The highest index you want to go to. Defaults to length(my_list).}

\item{other_args_str}{A string of the other arguments you want sent to mapply, formatted like "arg1, arg2, arg3". You need the commas. See examples. Defaults to NULL.}
}
\value{
The result of the mapply call.
}
\description{
This is probably dangerous. It generates and evaluates an mapply call for a nested list, which is really useful if you don't know how long your list will be before runtime.
}
\examples{
list1 = list(a=data.frame(a1=rep(1.1, 10), a2=rep(1.2, 10)),
	b=data.frame(b1=rep(1.1, 10), b2=rep(1.2, 10)))
list2 = list(a=data.frame(a1=rep(2.1, 10), a2=rep(2.2, 10)),
	b=data.frame(b1=rep(2.1, 10), b2=rep(2.2, 10)))
list3 = list(a=data.frame(a1=rep(3.1, 10), a2=rep(3.2, 10)),
	b=data.frame(b1=rep(3.1, 10), b2=rep(3.2, 10)))

big_list = list(list1, list2, list3)

# USE.NAMES makes no difference here, included to demonstrate
#	other_args_str syntax
combined = fccu_mapply_evil("rbind", big_list,
	other_args_str="SIMPLIFY=FALSE, USE.NAMES=FALSE")

# combined produces the same result as having written:
mapply(rbind, big_list[[1]], big_list[[2]],
big_list[[3]], SIMPLIFY=FALSE, USE.NAMES=FALSE)

# or having written this:
mapply(rbind, list1, list2, list3, SIMPLIFY=FALSE,
USE.NAMES=FALSE)

# the result of each is a list of 2 data.frames:
#	[[1]] is all the a's rbind-ed together
#	[[2]] is all the b's rbind-ed together

}
