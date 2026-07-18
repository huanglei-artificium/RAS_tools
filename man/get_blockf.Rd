\name{get_blockf}
\alias{get_blockf}
\title{Compute block-based rare allele sharing statistics}
\usage{
get_blockf(
  frequency_files_list,
  pop_left1,
  pop_right1,
  asc_outgroup,
  base_dir = NULL,
  position_file = NA,
  pop_def_file = NA,
  pop_left2 = NA,
  pop_left_fillna = NA,
  pop_right2 = NA,
  pop_right_fillna = NA,
  jackknife = NA,
  is_ascertained = FALSE,
  asc_refpop = NA,
  asc_freq_table = NA,
  all_sites_block = TRUE,
  refpop_max_miss = NA,
  refpop_remove_homo = FALSE,
  outgroup_max_hetero = 0,
  outgroup_round_freq = FALSE
)
}
\description{
Reads population frequency/genotype inputs, applies optional site filters and
ascertainment conditions, then returns block-level rare allele sharing statistics.
}
\arguments{
\item{frequency_files_list}{Path to a tab-delimited frequency-file manifest.}
\item{base_dir}{Optional base directory for relative file paths. When set,
relative paths passed to this function and relative paths in the manifest are
resolved against \code{base_dir}. When omitted, manifest entries are resolved
against the manifest file's directory.}
\item{pop_left1}{Comma-separated left population names, or a path to a file with one population per line.}
\item{pop_right1}{Comma-separated right population names, or a path to a file with one population per line.}
\item{asc_outgroup}{Outgroup population name.}
\item{position_file}{Optional two-column \code{CHR}, \code{POS} site subset file.}
\item{pop_def_file}{Optional custom population definition file.}
\item{pop_left2}{Optional secondary left population set for RASD contrasts.}
\item{pop_left_fillna}{Optional missing-data replacement for left populations.}
\item{pop_right2}{Optional secondary right population set for RASD contrasts.}
\item{pop_right_fillna}{Optional missing-data replacement for right populations.}
\item{jackknife}{Either \code{"CHR"} or an integer genomic window length.}
\item{is_ascertained}{Whether to compute ascertained frequency-bin results.}
\item{asc_refpop}{Reference population for derived-frequency ascertainment.}
\item{asc_freq_table}{Frequency-bin table with \code{min_freq} and \code{max_freq}.}
\item{all_sites_block}{Retained for compatibility; summarization uses this option in \code{get_f_from_fblock()}.}
\item{refpop_max_miss}{Maximum tolerated missingness in \code{asc_refpop}.}
\item{refpop_remove_homo}{Whether to exclude homozygous reference sites.}
\item{outgroup_max_hetero}{Maximum tolerated outgroup heterozygosity.}
\item{outgroup_round_freq}{Whether to round outgroup frequency to 0 or 1.}
}
\value{
A data frame/data table of block-level statistics.
}
