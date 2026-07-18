#' Compute a weighted mean and jackknife standard error
#'
#' @param Nr_sites Numeric vector of block weights.
#' @param f Numeric vector of block-level statistics.
#'
#' @return A list with number of blocks, total sites, weighted statistic, and
#'   jackknife standard error.
#' @export
get_mean_and_std <- function(Nr_sites, f) {
    if ((length1 <- length(Nr_sites)) != (length2 <- length(f))) {
        stop("Error: unequal length between vectors on number of sites and f.", call. = FALSE)
    }

    jackknife_weight_sum <- sum(Nr_sites)
    f_value <- sum(f * Nr_sites) / jackknife_weight_sum
    list(
        Nr_blocks = length1,
        Nr_sites = jackknife_weight_sum,
        f_value = f_value,
        std = sqrt(mean(Nr_sites / (jackknife_weight_sum - Nr_sites) * (f_value - f)^2))
    )
}

#' Summarize block-level statistics
#'
#' @param block_object A block table returned by [get_blockf()].
#' @param block_files Optional vector of block-table files to read and combine.
#' @param by_columns_extra Optional additional grouping columns.
#' @param all_sites_block Use `Nr_sites_all` instead of `Nr_sites_nonmissing`.
#' @param out_file Optional output path for the summary table.
#' @param return_result Whether to return the summary table.
#'
#' @return Invisibly writes a table when `out_file` is supplied. Returns the
#'   summary table when `return_result` is `TRUE`.
#' @export
get_f_from_fblock <- function(
    block_object = NA,
    block_files = NA,
    by_columns_extra = NULL,
    all_sites_block = FALSE,
    out_file = NA,
    return_result = FALSE
) {
    if (!is.na(block_object)[1]) {
        fblock <- block_object
    } else {
        fblock <- purrr::map(block_files, data.table::fread, header = TRUE, sep = "\t") %>%
            dplyr::bind_rows()
    }

    by_columns <- base::colnames(fblock) %>%
        .[. %in% c("min_freq", "max_freq", "asc_outgroup", "pop_left1", "pop_left2", "pop_right1", "pop_right2")] %>%
        c(., by_columns_extra)
    weight_column <- if (all_sites_block) "Nr_sites_all" else "Nr_sites_nonmissing"
    if (length(by_columns)) {
        group_keys <- lapply(fblock[by_columns], function(x) {
            x <- as.character(x)
            x[is.na(x)] <- "<NA>"
            x
        })
        groups <- split(fblock, group_keys, drop = TRUE)
    } else {
        groups <- list(all = fblock)
    }
    f3 <- base::do.call(rbind, lapply(groups, function(x) {
        stats <- get_mean_and_std(x[[weight_column]], x[["f"]])
        cbind(
            x[1, by_columns, drop = FALSE],
            data.frame(
                Nr_blocks = stats$Nr_blocks,
                Nr_sites = stats$Nr_sites,
                f_value = stats$f_value,
                std = stats$std
            )
        )
    }))
    row.names(f3) <- NULL

    if (!is.na(out_file)) {
        utils::write.table(f3, out_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }

    if (return_result) {
        return(f3)
    }
    invisible(f3)
}
