
# private
# single frequency file
get_frequency_file_info <- function(frequency_file, header_info = NA)
{
    is_geno <- tools::file_ext(frequency_file) == "geno"
    if (base::is.na(header_info))
        header_info <- base::paste0(tools::file_path_sans_ext(frequency_file), if (is_geno) ".ind" else ".header")
    
    header <- ( if (base::file.exists(header_info)) utils::read.table(header_info, header = FALSE, sep = "\t")[[1]] else base::strsplit(header_info, ",")[[1]] ) %>% base::strsplit(split="[()]")
    return(data.frame(frequency_file = frequency_file, field = 1:base::length(header), column_name = base::sapply(header, `[`, 1), denominator = if (is_geno) 2 else base::as.numeric(base::sapply(header, `[`, 2))))
}


# private
is_two_column_pop <- function(pop, tbl)
{
    frac2_name <- paste0(pop, c("~num", "~deno"))
    if (! all(frac2_name %in% rownames(tbl)))
        return(FALSE)
    
    deno <- tbl[frac2_name, "denominator"] %>% base::unique()
    if (base::length(deno) > 1 || is.na(deno[[1]]))
    {
        message("Population ", pop, " doesn't have valid definition of total count in its numerator (~num) and denominator (~deno).")
        return(FALSE)
    }
    return(TRUE)
}

# private
# multiple frequency files
get_frequency_files_info <- function(frequency_files_list, full_info)
{
    if (!file.exists(frequency_files_list))
    {
        message("Frequency files list", frequency_files_list, "does not exist.")
        stop()
    }
    frequency_files <- utils::read.table(frequency_files_list, header = FALSE, sep = "\t", col.names = c("frequency_file", "header_info"), fill = TRUE, na.strings = c("NA", ""))
    # frequency_file; field; column_name; denominator
    frequency_files_info <- purrr::pmap(frequency_files, get_frequency_file_info) %>% do.call(base::rbind, .)

    column_names_duplicated <- frequency_files_info$column_name %>% .[base::duplicated(.)] %>% base::unique()
    if (base::length(column_names_duplicated))
    {
        message("The following column(s) in frequency file(s) is/are duplicated: ", base::paste(column_names_duplicated, collapse=", "), ".")
        stop()
    }

    frequency_files_info <- frequency_files_info %>% dplyr::mutate(pop_name = .$column_name %>% base::strsplit(split="~") %>% base::sapply(`[`, 1) )
    # rename the row names to "column names" of frequency files.
    row.names(frequency_files_info) <- frequency_files_info$column_name
    frequency_files_info$column_name <- NULL

    full_info$columns <- frequency_files_info %>% base::split(base::ifelse(.$pop_name %in% c("CHR", "POS", "ALT", "REF"), "posi", "pop") )
    
    # pick out the position columns specifically.
    full_info$columns$posi <- full_info$columns$posi %>% .[base::row.names(.) %in% c("CHR", "POS", "ALT", "REF"), ]

    # full_info$columns$pop : frequency_files_info
    pop_name_count <- base::table(full_info$columns$pop$pop_name)

    full_info$pop_two_column <- pop_name_count %>% .[. == 2] %>% base::names() %>% .[base::sapply(., is_two_column_pop, tbl=full_info$columns$pop)]
    full_info$pop_one_column <- pop_name_count %>% .[. == 1] %>% base::names() %>% .[. %in% base::row.names(full_info$columns$pop)]

    full_info$pop_original <- c(full_info$pop_two_column, full_info$pop_one_column)

    # filter the populations. c(pop_two_column, pop_one_column)
    full_info$columns$pop <- full_info$columns$pop %>% .[.$pop_name %in% full_info$pop_original, ]
}


# could be optimized on recording. for example, in the case using defined population to define defined population
# could add the judgment: whether freq data (denominator 1) is included. if so, delete this item. update: do not need to message error and delete. message warning now.
get_pop_def <- function(pop_def_file, pop_input)
{
    if (is.na(pop_def_file))
    {
        return(list())
    }

    if (!file.exists(pop_def_file))
    {
        message("The population definition file ", pop_def_file , " does not exist.")
        return(list())
    }

    message("Process information on customized population(s) in ", pop_def_file, " .")

    pop_def_table_raw <- utils::read.table(pop_def_file, header = FALSE, sep = "\t")[[1]] %>% base::strsplit(split=":")
    # can use %>% .[! lapply(., `[`, 1) %in% names(pop_info)] to filter, but not applied anymore

    pop_def <- list()
    for (pop_def_item in pop_def_table_raw)
    {
        message("Customize population ", pop_def_item[1], ":")

        pop_exist <- c(pop_input, names(pop_def))
        if (pop_def_item[1] %in% pop_exist)
        {
            message("Population ", pop_def_item[1], " already exists! Skip.")
            next
        }
        pop_comp <- pop_def_item[2] %>% base::strsplit(split=",") %>% .[[1]] %>% base::unique()
        pop_comp_noexist <- pop_comp %>% .[! . %in% pop_exist]
        if (length(pop_comp_noexist) > 0)
        {
            message("Population(s) not existing: ", base::paste(pop_comp_noexist, collapse=", "), ". Not included.")
            pop_comp <- pop_comp %>% .[. %in% pop_exist]
        }
        message(pop_def_item[1], ": ", base::paste(pop_comp, collapse=", "), ".")

        pop_def[[pop_def_item[1]]] <- pop_comp
    }
    return(pop_def)
}

# ...
pop_to_freq <- function(pop_all, column_info, mask = NULL)
{
    # optimize: use read_lines() ?
    pop_all_file_info <- column_info %>% .[.$pop_name %in% pop_all, ]

    file_all <- pop_all_file_info$frequency_file %>% base::unique()

    geno_all <- list()
    for (file in file_all)
    {
        pop_file_info <- pop_all_file_info %>% filter(frequency_file == file)
        field <- pop_file_info$field
        pop <- base::row.names(pop_file_info)

        is_geno <- tools::file_ext(file) == "geno"
        geno <- if (is_geno) readr::read_fwf(file, col_positions=readr::fwf_positions(field, field, col_names=pop), col_types=readr::cols(.default="i"), na="9") else data.table::fread(file=file, header=FALSE, select=field, na.strings=c("-1", "*"), col.names=pop)
        if (!is.null(mask))
            geno <- geno[mask, ]
        geno_all <- geno_all %>% append(geno)
    }

    return(geno_all)
}

##
get_freq <- function(pops, lth, mask, full_info)
{
    # names(full_info$freq) == full_info$pops[!is.na(lth)] always holds

    pop_exist <- pops %>% .[. %in% names(full_info$freq)]
    # (1) already in full_info$freq (2) length is updated
    pop_to_mask <- pop_exist %>% .[full_info$pops[., "lth"] < lth]

    # could use lapply to optimize:
    # my_list[c("key2", "key4")] <- lapply(my_list[c("key2", "key4")], function(x) x * 2)
    for (pop in pop_to_mask)
    {
        full_info$freq[[pop]] <- full_info$freq[[pop]][mask]
        if (full_info$pops[pop, "keep_frac"])
        {
            full_info$frac[[pop]] <- full_info$frac[[pop]][mask, ]
        }
    }

    pop_to_read <- pops %>% .[! . %in% pop_exist]
    pop_defined_to_read <- pop_to_read %>% .[. %in% full_info$pop_newdefined]
    # having already removed existing pops
    # process pop_original + pop_defined -> pop_original
    pop_original_to_read <- c(pop_to_read %>% .[. %in% full_info$pop_original], full_info$defs[pop_defined_to_read] %>% base::unlist() %>% .[! . %in% pop_exist]) %>% base::unique()

    # read data in  a list
    columns <- pop_to_freq(pop_original_to_read, full_info$columns$pop, mask)

    for (pop in pop_original_to_read)
    {
        # faster? mutate or cbind?
        if (full_info$pops[pop, "Nr_columns"] == 2) # frac2
        {
            full_info$freq[[pop]] <- ( columns[[paste0(pop, "~num")]] / columns[[paste0(pop, "~deno")]] ) %>% replace(is.nan(.), NA)
            if (full_info$pops[pop, "keep_frac"])
            {
                full_info$frac[[pop]] <- data.frame(num = columns[[paste0(pop, "~num")]], deno = columns[[paste0(pop, "~deno")]])
            }
        } else # frac1 || freq, do not distinguish total>1 or total == 1 now. Problem: what if haploid calling with total 1? it's not frequency and can be added!
        {
            full_info$freq[[pop]] <- columns[[pop]] / full_info$pops[pop, "total"]
            if (full_info$pops[pop, "keep_frac"])
            {
                full_info$frac[[pop]] <- data.frame(num = dplyr::coalesce(columns[[pop]], 0), deno = ifelse(is.na(columns[[pop]]), 0, full_info$pops[pop, "total"]))
            }
            # already have NA in columns[[pop]]

        }
        
    }
    rm(columns)
    
    for (pop in pop_defined_to_read)
    {
        full_info$frac[[pop]] <- base::Reduce(`+`, full_info$frac[full_info$defs[[pop]]])
        full_info$freq[[pop]] <- (full_info$frac[[pop]]$num / full_info$frac[[pop]]$deno) %>% replace(is.nan(.), NA)
    }

    full_info$pops[c(pop_to_mask, pop_original_to_read, pop_defined_to_read), "lth"] <- lth

}


###
get_blockfsum <- function(asc_outgroup, pop_left, pop_right, full_info)
{
    asc_outgroup_freq <- full_info$freq[[asc_outgroup]]
    pop_left_freq <- full_info$freq[[pop_left]] %>% {if (!is.na(full_info$pop_left_fillna)) fill_na_with_ref(., ref=if (full_info$pop_left_fillna=="Outgroup") full_info$outgroup_freq else 0) else .}
    pop_right_freq <- full_info$freq[[pop_right]] %>% {if (!is.na(full_info$pop_right_fillna)) fill_na_with_ref(., ref=if (full_info$pop_right_fillna=="Outgroup") full_info$outgroup_freq else 0) else .}

    f_value <- (pop_left_freq - asc_outgroup_freq) * (pop_right_freq - asc_outgroup_freq)

    if (!is.na(asc_freq_table))
    {
        fsum_block_asc_cutoff <- data.table::data.table(block_id=full_info$block_id, derived_freq_level=full_info$derived_freq_level, f_value=f_value)[, .(Nr_sites_all=.N, Nr_sites_nonmissing=sum(!is.na(f_value)), fsum=sum(f_value, na.rm=TRUE)), .(block_id, derived_freq_level)][base::order(block_id, derived_freq_level)][, .(derived_freq_level_max=derived_freq_level, Nr_sites_all=cumsum(Nr_sites_all), Nr_sites_nonmissing=cumsum(Nr_sites_nonmissing), fsum=cumsum(fsum)), .(block_id)]

        f_value_block_table <- data.table::data.table()

        for (i in 1:nrow(full_info$asc_conds))
        {
            min_freq_level <- full_info$asc_conds$min_freq_level[i]
            max_freq_level <- full_info$asc_conds$max_freq_level[i]

            block_remain <- full_info$block_id_unique
            max_freq_level_curr <- max_freq_level
            fsum_block_max <- data.table::data.table()
            repeat
            {
                fsum_block_max <- fsum_block_max %>% rbind(fsum_block_asc_cutoff %>% filter(block_id %in% block_remain & derived_freq_level_max == max_freq_level_curr))
                block_remain <- block_remain %>% base::setdiff(fsum_block_max$block_id)
                if (length(block_remain) == 0)
                    break
                max_freq_level_curr <- max_freq_level_curr - 1
            }
            fsum_block_max <- fsum_block_max[base::order(block_id)][, c("Nr_sites_all", "Nr_sites_nonmissing", "fsum")]
            
            fsum_block_min <- if (min_freq_level > 1) fsum_block_asc_cutoff[fsum_block_asc_cutoff$derived_freq_level_max == min_freq_level - 1, c("Nr_sites_all", "Nr_sites_nonmissing", "fsum")] else 0

            fsum_block <- fsum_block_max - fsum_block_min
            
            f_value_block_table <- f_value_block_table %>% rbind(data.table::data.table(min_freq = full_info$asc_conds$min_freq[i], max_freq = full_info$asc_conds$max_freq[i], block_id = full_info$block_id_unique, Nr_sites_all = fsum_block$Nr_sites_all, Nr_sites_nonmissing = fsum_block$Nr_sites_nonmissing, f = fsum_block$fsum / fsum_block$Nr_sites_nonmissing))


        }
    } else
    {
        fsum_block <- data.table::data.table(block_id=full_info$block_id, f_value=f_value)[, .(Nr_sites_all=.N, Nr_sites_nonmissing=sum(!is.na(f_value)), fsum=sum(f_value, na.rm=TRUE)), .(block_id)][base::order(block_id)]
        
        f_value_block_table <- data.table::data.table(block_id = full_info$block_id_unique, Nr_sites_all = fsum_block$Nr_sites_all, Nr_sites_nonmissing = fsum_block$Nr_sites_nonmissing, f = fsum_block$fsum / fsum_block$Nr_sites_nonmissing)

    }

    # return(f_value_block_table[, c("asc_outgroup", "pop_left", "pop_right") := data.frame(asc_outgroup, pop_left, pop_right)])
    return(f_value_block_table)
}


##
get_RASD_from_RAS <- function(asc_outgroup, pop_left1, pop_left2=NULL, pop_right1, pop_right2=NULL, full_info)
{
    id_L1_R1 <- full_info$all_stats_RAS %>% .[.$asc_outgroup==asc_outgroup & .$pop_left==pop_left1 & .$pop_right==pop_right1, ] %>% row.names() %>% as.numeric()
    f_value_block_table <- full_info$RAS_table_all[[id_L1_R1]][, c("min_freq", "max_freq", "block_id", "Nr_sites_all")]
    f_L1_R1 <- full_info$RAS_table_all[[id_L1_R1]]$f
    if (is.null(pop_left2))
    {
        id_L1_R2 <- full_info$all_stats_RAS %>% .[.$asc_outgroup==asc_outgroup & .$pop_left==pop_left1 & .$pop_right==pop_right2, ] %>% row.names() %>% as.numeric()
        f_L1_R2 <- full_info$RAS_table_all[[id_L1_R2]]$f
        f_value_block_table$f <- f_L1_R1 - f_L1_R2
    } else if (is.null(pop_right2))
    {
        id_L2_R1 <- full_info$all_stats_RAS %>% .[.$asc_outgroup==asc_outgroup & .$pop_left==pop_left2 & .$pop_right==pop_right1, ] %>% row.names() %>% as.numeric()
        f_L2_R1 <- full_info$RAS_table_all[[id_L2_R1]]$f
        f_value_block_table$f <- f_L1_R1 - f_L2_R1
    } else
    {
        id_L1_R2 <- full_info$all_stats_RAS %>% .[.$asc_outgroup==asc_outgroup & .$pop_left==pop_left1 & .$pop_right==pop_right2, ] %>% row.names() %>% as.numeric()
        id_L2_R1 <- full_info$all_stats_RAS %>% .[.$asc_outgroup==asc_outgroup & .$pop_left==pop_left2 & .$pop_right==pop_right1, ] %>% row.names() %>% as.numeric()
        id_L2_R2 <- full_info$all_stats_RAS %>% .[.$asc_outgroup==asc_outgroup & .$pop_left==pop_left2 & .$pop_right==pop_right2, ] %>% row.names() %>% as.numeric()
        f_L1_R2 <- full_info$RAS_table_all[[id_L1_R2]]$f
        f_L2_R1 <- full_info$RAS_table_all[[id_L2_R1]]$f
        f_L2_R2 <- full_info$RAS_table_all[[id_L2_R2]]$f
        f_value_block_table$f <- f_L1_R1 - f_L1_R2 - f_L2_R1 + f_L2_R2
    }
    return(f_value_block_table)
    # return(f_value_block_table[, c("asc_outgroup", "pop_left1", "pop_left2", "pop_right1", "pop_right2") := data.frame(asc_outgroup, pop_left1, pop_left2, pop_right1, pop_right2)])
}