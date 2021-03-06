% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fc_stat_test.R
\name{fcst_fochs}
\alias{fcst_fochs}
\title{Calculate Fold Change For Fc Array Features}
\usage{
fcst_fochs(fc, comp_name, method = median)
}
\arguments{
\item{fc}{The Fc Array data frame. Must have two groups present.}

\item{comp_name}{A string representing the name of the group comparison. This should be a contraction of the groups you are comparing, separated by an underscore. The first group in this string will be the numerator in the fold change ratios.}

\item{method}{What method will be used to calculate fold change. The default is the median of each group's values.}
}
\value{
A vector of fold change values of the same length as fccu_first_feat_col(fc):ncol(fc).
}
\description{
For each feature in an Fc Array data frame, calculate the fold change between the measurements of one group and another using method. Which group is the numerator in the fold change ratio is determined by the structure of comp_name.
}
