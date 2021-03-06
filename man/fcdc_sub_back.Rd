% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fc_data_clean.R
\name{fcdc_sub_back}
\alias{fcdc_sub_back}
\title{Background Signal Subtraction}
\usage{
fcdc_sub_back(fc, base_ids)
}
\arguments{
\item{fc}{The Fc Array data frame.}

\item{base_ids}{The row indices in fc that correspond to the subjects who represent baseline (say, those in the PLACEBO group).}
}
\value{
The Fc Array data frame, with fc[base_ids,] having all feature measurements set to 0, and with non-baseline subjects having had their MFIs subtracted by the median of the baseline group's MFI values for each feature, and then scaled by the sd of those baseline MFI values.
}
\description{
Perform background signal subtraction from an Fc Array data frame.
}
