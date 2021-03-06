% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fc_common_utils.R
\name{fccu_as_evil}
\alias{fccu_as_evil}
\title{Black Magic Add/Subtract for Lists}
\usage{
fccu_as_evil(my_op_str, my_list, max_idx = length(my_list))
}
\arguments{
\item{my_op_str}{A string representation of the arithmetic operator you want to use, "+" or "-".}

\item{my_list}{The list you want to perform this sorcery on.}

\item{max_idx}{The highest index you want to go to. Defaults to length(my_list).}
}
\value{
The result of the mungo addition or subtraction.
}
\description{
This is probably dangerous. Combines the objects in my_list according to my_op_str. See example.
}
\examples{
> data1 = data.frame(a=c(1,2,3), b=c(4,5,6))
> data2 = data.frame(a=c(9,8,7), b=c(6,5,4))
> data3 = data.frame(a=c(1,1,1), b=c(2,2,2))

> fccu_as_evil("+", list(data1, data2, data3))
   a  b
1 11 12
2 11 12
3 11 12

# evaluates to the equivalent of writing out:
> data1 + data2 + data3

}
