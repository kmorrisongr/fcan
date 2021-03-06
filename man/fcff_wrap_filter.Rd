% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fc_feat_filt.R
\name{fcff_wrap_filter}
\alias{fcff_wrap_filter}
\title{Higher-Level Interface For Filtering Features}
\usage{
fcff_wrap_filter(
  fc,
  results_dir,
  keep_filter,
  discard_filter,
  k_behavior,
  d_behavior,
  diff_opts = "NULL"
)
}
\arguments{
\item{fc}{The Fc Array data frame}

\item{results_dir}{A string representing the directory where you want to store the results of your analyses.}

\item{keep_filter}{A vector of strings of features that you want to retain in fc. If no keeping is to be done, pass an empty vector c().}

\item{discard_filter}{A vector of strings of features that you want to remove from fc. If no discarding is to be done, pass an empty vector c().}

\item{k_behavior}{A string "permissive" or "strict" controlling keeping behavior. Can also be set to "differs" to perform fcff_feats_diff. See ?fcff_filter_features for more.}

\item{d_behavior}{A string "permissive" or "strict" controlling discarding behavior. Can also be set to "differs" to perform fcff_feats_diff. See ?fcff_filter_features for more.}

\item{diff_opts}{A list of subparameters for fcff_feats_diff. Pass "NULL" if you don't want to do any fcff_feats_diff things. See ?fcff_feats_diff for more.}
}
\value{
A list containing the filtered Fc Array data frame (accessed by $fc), and the new results_dir string based on the filtering performed (accessed by $results_dir). If both keeping and discarding are performed, results_dir will be of the form complex_YYYYMMDD_HHMMSS, and the directory will contain a txt file with the filtering that was performed.
}
\description{
Filter features in such a way as to set up experiments. Will filter feautres/perform fcff_feats_diff as necessary, and will create files/directories (as necessary) to facilitate result reproduction and keep different versions of fc analysis sequestered.
}
\examples{
results_dir = "/home/me/science/experiments/"
filtered = fcff_wrap_filter(fc, results_dir, c("gp41", "gp120"), c(), "permissive", "NULL")
fc = filtered$fc
# is now "/home/me/science/experiments/only_gp41_gp120"
results_dir = filtered$results_dir

}
