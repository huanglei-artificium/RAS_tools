\name{get_f_from_fblock}
\alias{get_f_from_fblock}
\title{Summarize block-level statistics}
\usage{
get_f_from_fblock(
  block_object = NA,
  block_files = NA,
  by_columns_extra = NULL,
  all_sites_block = FALSE,
  out_file = NA,
  return_result = FALSE
)
}
\arguments{
\item{block_object}{A block table returned by \code{get_blockf()}.}
\item{block_files}{Optional vector of block-table files to read and combine.}
\item{by_columns_extra}{Optional additional grouping columns.}
\item{all_sites_block}{Use \code{Nr_sites_all} instead of \code{Nr_sites_nonmissing}.}
\item{out_file}{Optional output path for the summary table.}
\item{return_result}{Whether to return the summary table.}
}
\description{
Summarize block-level statistics into weighted means and jackknife standard errors.
}
\value{
Returns the summary table when \code{return_result} is \code{TRUE}; otherwise
returns it invisibly.
}
