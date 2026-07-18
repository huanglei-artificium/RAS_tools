get_blockRAS_asc <- function(asc_outgroup, pop_left, pop_right, full_info) {
    asc_outgroup_freq <- full_info$freq[[asc_outgroup]]
    pop_left_freq <- full_info$freq[[pop_left]] %>%
        {if (!is.na(full_info$pop_left_fillna)) fill_na_with_ref(., ref = if (full_info$pop_left_fillna == "Outgroup") full_info$outgroup_freq else 0) else .}
    pop_right_freq <- full_info$freq[[pop_right]] %>%
        {if (!is.na(full_info$pop_right_fillna)) fill_na_with_ref(., ref = if (full_info$pop_right_fillna == "Outgroup") full_info$outgroup_freq else 0) else .}

    f_value <- (asc_outgroup_freq - pop_left_freq) * (asc_outgroup_freq - pop_right_freq)

    site_values <- data.frame(
        block_id = full_info$block_id,
        derived_freq_level = full_info$derived_freq_level,
        f_value = f_value
    )
    fsum_block_asc <- base::do.call(rbind, lapply(
        split(site_values, list(site_values$block_id, site_values$derived_freq_level), drop = TRUE),
        function(x) {
            data.frame(
                block_id = x$block_id[1],
                derived_freq_level = x$derived_freq_level[1],
                Nr_sites_all = nrow(x),
                Nr_sites_nonmissing = sum(!is.na(x$f_value)),
                fsum = sum(x$f_value, na.rm = TRUE)
            )
        }
    ))
    fsum_block_asc <- fsum_block_asc[base::order(fsum_block_asc$block_id, fsum_block_asc$derived_freq_level), ]
    fsum_block_asc_cutoff <- base::do.call(rbind, lapply(
        split(fsum_block_asc, fsum_block_asc$block_id),
        function(x) {
            x <- x[base::order(x$derived_freq_level), ]
            data.frame(
                block_id = x$block_id,
                derived_freq_level_max = x$derived_freq_level,
                Nr_sites_all = cumsum(x$Nr_sites_all),
                Nr_sites_nonmissing = cumsum(x$Nr_sites_nonmissing),
                fsum = cumsum(x$fsum)
            )
        }
    ))

    f_value_block_table <- data.table::data.table()

    for (i in seq_len(nrow(full_info$asc_conds))) {
        min_freq_level <- full_info$asc_conds$min_freq_level[i]
        max_freq_level <- full_info$asc_conds$max_freq_level[i]

        block_remain <- full_info$block_id_unique
        max_freq_level_curr <- max_freq_level
        fsum_block_max <- data.table::data.table()
        repeat {
            fsum_block_max <- fsum_block_max %>%
                rbind(fsum_block_asc_cutoff %>% filter(block_id %in% block_remain & derived_freq_level_max == max_freq_level_curr))
            block_remain <- block_remain %>% base::setdiff(fsum_block_max$block_id)
            if (length(block_remain) == 0) {
                break
            }
            max_freq_level_curr <- max_freq_level_curr - 1
        }
        fsum_block_max <- fsum_block_max[base::order(fsum_block_max$block_id), c("Nr_sites_all", "Nr_sites_nonmissing", "fsum")]

        fsum_block_min <- if (min_freq_level > 1) {
            fsum_block_min <- fsum_block_asc_cutoff[
                fsum_block_asc_cutoff$derived_freq_level_max == min_freq_level - 1,
                c("block_id", "Nr_sites_all", "Nr_sites_nonmissing", "fsum")
            ]
            fsum_block_min <- fsum_block_min[match(full_info$block_id_unique, fsum_block_min$block_id), c("Nr_sites_all", "Nr_sites_nonmissing", "fsum")]
            fsum_block_min[is.na(fsum_block_min)] <- 0
            fsum_block_min
        } else {
            0
        }

        fsum_block <- fsum_block_max - fsum_block_min

        f_value_block_table <- f_value_block_table %>%
            rbind(data.table::data.table(
                min_freq = full_info$asc_conds$min_freq[i],
                max_freq = full_info$asc_conds$max_freq[i],
                block_id = full_info$block_id_unique,
                Nr_sites_all = fsum_block$Nr_sites_all,
                Nr_sites_nonmissing = fsum_block$Nr_sites_nonmissing,
                f = fsum_block$fsum / fsum_block$Nr_sites_nonmissing
            ))
    }

    f_value_block_table
}

get_RASD_from_RAS <- function(asc_outgroup, pop_left1, pop_left2 = NULL, pop_right1, pop_right2 = NULL, full_info) {
    id_L1_R1 <- full_info$all_stats_RAS %>%
        .[.$asc_outgroup == asc_outgroup & .$pop_left == pop_left1 & .$pop_right == pop_right1, ] %>%
        row.names() %>%
        as.numeric()
    f_value_block_table <- full_info$RAS_table_all[[id_L1_R1]][, c("min_freq", "max_freq", "block_id", "Nr_sites_all", "Nr_sites_nonmissing", "f")]

    if (!is.null(pop_left2)) {
        id_L2_R1 <- full_info$all_stats_RAS %>%
            .[.$asc_outgroup == asc_outgroup & .$pop_left == pop_left2 & .$pop_right == pop_right1, ] %>%
            row.names() %>%
            as.numeric()
        f_value_block_table$f <- f_value_block_table$f - full_info$RAS_table_all[[id_L2_R1]]$f
    }

    if (!is.null(pop_right2)) {
        id_L1_R2 <- full_info$all_stats_RAS %>%
            .[.$asc_outgroup == asc_outgroup & .$pop_left == pop_left1 & .$pop_right == pop_right2, ] %>%
            row.names() %>%
            as.numeric()
        f_value_block_table$f <- f_value_block_table$f - full_info$RAS_table_all[[id_L1_R2]]$f
    }

    if (!is.null(pop_left2) && !is.null(pop_right2)) {
        id_L2_R2 <- full_info$all_stats_RAS %>%
            .[.$asc_outgroup == asc_outgroup & .$pop_left == pop_left2 & .$pop_right == pop_right2, ] %>%
            row.names() %>%
            as.numeric()
        f_value_block_table$f <- f_value_block_table$f + full_info$RAS_table_all[[id_L2_R2]]$f
    }
    f_value_block_table
}

get_blockf_nonasc <- function(asc_outgroup, pop_left1, pop_left2 = NULL, pop_right1, pop_right2 = NULL, full_info) {
    asc_outgroup_freq <- full_info$freq[[asc_outgroup]]
    pop_left1_freq <- full_info$freq[[pop_left1]] %>%
        {if (!is.na(full_info$pop_left_fillna)) fill_na_with_ref(., ref = if (full_info$pop_left_fillna == "Outgroup") full_info$outgroup_freq else 0) else .}
    if (!is.null(pop_left2)) {
        pop_left2_freq <- full_info$freq[[pop_left2]] %>%
            {if (!is.na(full_info$pop_left_fillna)) fill_na_with_ref(., ref = if (full_info$pop_left_fillna == "Outgroup") full_info$outgroup_freq else 0) else .}
        pop_left_freq <- pop_left2_freq - pop_left1_freq
    } else {
        pop_left_freq <- asc_outgroup_freq - pop_left1_freq
    }

    pop_right1_freq <- full_info$freq[[pop_right1]] %>%
        {if (!is.na(full_info$pop_right_fillna)) fill_na_with_ref(., ref = if (full_info$pop_right_fillna == "Outgroup") full_info$outgroup_freq else 0) else .}
    if (!is.null(pop_right2)) {
        pop_right2_freq <- full_info$freq[[pop_right2]] %>%
            {if (!is.na(full_info$pop_right_fillna)) fill_na_with_ref(., ref = if (full_info$pop_right_fillna == "Outgroup") full_info$outgroup_freq else 0) else .}
        pop_right_freq <- pop_right2_freq - pop_right1_freq
    } else {
        pop_right_freq <- asc_outgroup_freq - pop_right1_freq
    }
    f_value <- pop_left_freq * pop_right_freq

    site_values <- data.frame(block_id = full_info$block_id, f_value = f_value)
    fsum_block <- base::do.call(rbind, lapply(split(site_values, site_values$block_id), function(x) {
        data.frame(
            block_id = x$block_id[1],
            Nr_sites_all = nrow(x),
            Nr_sites_nonmissing = sum(!is.na(x$f_value)),
            fsum = sum(x$f_value, na.rm = TRUE)
        )
    }))
    fsum_block <- fsum_block[match(full_info$block_id_unique, fsum_block$block_id), ]

    data.table::data.table(
        block_id = full_info$block_id_unique,
        Nr_sites_all = fsum_block$Nr_sites_all,
        Nr_sites_nonmissing = fsum_block$Nr_sites_nonmissing,
        f = fsum_block$fsum / fsum_block$Nr_sites_nonmissing
    )
}
