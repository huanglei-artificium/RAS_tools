\name{get_mean_and_std}
\alias{get_mean_and_std}
\title{Compute a weighted mean and jackknife standard error}
\usage{
get_mean_and_std(Nr_sites, f)
}
\arguments{
\item{Nr_sites}{Numeric vector of block weights.}
\item{f}{Numeric vector of block-level statistics.}
}
\description{
Compute a weighted mean and jackknife standard error from block-level statistics.
}
\value{
A list with number of blocks, total sites, weighted statistic, and jackknife
standard error.
}
