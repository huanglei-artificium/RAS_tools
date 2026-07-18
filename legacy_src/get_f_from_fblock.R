get_mean_and_std <- function(Nr_sites, f)
{
    if ( (length1 <- length(Nr_sites)) != (length2 <- length(f)) )
        stop("Error: unequal length between vectors on number of sites and f.")

    jackknife_weight_sum <- sum(Nr_sites)
    f_value <- sum(f * Nr_sites) / jackknife_weight_sum
    return(list(Nr_blocks = length1, Nr_sites = jackknife_weight_sum, f_value = f_value, std = sqrt(mean(Nr_sites / (jackknife_weight_sum - Nr_sites) * (f_value - f) ^ 2))))
}


get_f_from_fblock <- function(block_object=NA, block_files=NA, by_columns_extra=NULL, all_sites_block=FALSE, out_file=NA, return_result=FALSE)
{
    # by_columns_extra : define extra columns other than standard columns.

    if (!is.na(block_object))
    {
        fblock <- block_object
    } else {
        fblock <- purrr::map(block_files, data.table::fread, header = TRUE, sep = "\t") %>% dplyr::bind_rows()
    }

    by_columns <- base::colnames(fblock) %>% .[. %in% c("min_freq", "max_freq", "asc_outgroup", "pop_left1", "pop_left2", "pop_right1", "pop_right2")] %>% c(., by_columns_extra)
    f3 <- fblock %>% .[, get_mean_and_std(if (all_sites_block) Nr_sites_all else Nr_sites_nonmissing, f), by = by_columns]

    if (!is.na(out_file))
        write.table(f3, out_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

    if (return_result)
        return(f3)
}