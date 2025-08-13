get_blockf <- function(
    position_file = NA,
    frequency_files_list,
    pop_def_file = NA,
    pop_left1,
    pop_left1_fillna = NA,
    pop_left2 = NA,
    pop_left2_fillna = NA,
    pop_right1,
    pop_right1_fillna = NA,
    pop_right2 = NA,
    pop_right2_fillna = NA,
    jackknife = NA,
    is_ascertained = FALSE,
    asc_refpop = NA,
    asc_freq_table = NA,
    all_sites_block = TRUE,
    refpop_max_miss = NA,
    refpop_remove_homo = NA,
    asc_outgroup,
    outgroup_max_hetero = 0,
    outgroup_round_freq = FALSE
)
{
    
frequency_files_info <- get_frequency_files_info(frequency_files_list)

pop_names_all <- frequency_files_info$column_name %>% .[! . %in% c("CHR", "POS", "ALT", "REF")] %>% base::strsplit(split="~") %>% base::sapply(`[`, 1) %>% base::unique()


pop_info <- list()
columns_selected <- c()


for (pop in pop_names_all)
{
    frac2_name <- paste0(pop, c("~num", "~deno"))
    if (all(frac2_name %in% frequency_files_info$column_name))
    {
        deno <- frequency_files_info %>% dplyr::filter(column_name %in% frac2_name) %>% .$denominator %>% base::unique()
        if (base::length(deno) > 1)
        {
            message("Population ", pop, " doesn't have the same total count in its numerator and denominator.")
            next
        }
        if (is.na(deno))
        {
            message("Population ", pop, " doesn't have the total count in its numerator or denominator.")
            next
        }
        pop_info[[pop]] <- list(data_type = "frac2", total = deno, lth = NA)
        columns_selected <- columns_selected %>% append(frac2_name)
        next
    }
    
    if (pop %in% frequency_files_info$column_name)
    {
        #
        deno <- frequency_files_info %>% dplyr::filter(column_name == pop) %>% .$denominator
        pop_info[[pop]] <- if (is.na(deno) || deno == 1) list(data_type = "freq", total = 1, lth = NA) else list(data_type = "frac1", total = deno, lth = NA)
        columns_selected <- columns_selected %>% append(pop)
        next
    }

    message("Population ", pop, " doesn't apply either fraction or frequency format.")
}

frequency_files_info_filtered <- frequency_files_info %>% dplyr::filter(column_name %in% c("CHR", "POS", "ALT", "REF", columns_selected))


pop_def <- get_pop_def(pop_def_file, names(pop_info))
for (pop in names(pop_def))
    pop_info[[pop]] <- list(data_type = "frac2", total = NA, lth = NA)


get_freq <- function(pop, lth, columns, keep_frac = FALSE, pop_info)
{
    #
    if (!is.na(pop_info[[pop]]$lth))
    {
        if (pop_info[[pop]]$lth == lth)
            return()
        if (pop_info[[pop]]$lth == "full" && lth == "asc")
        {
            if (keep_frac)
            {
                pop_info[[pop]]$num <<- pop_info[[pop]]$num[filter_all]
                pop_info[[pop]]$deno <<- pop_info[[pop]]$deno[filter_all]
            }
            pop_info[[pop]]$freq <<- pop_info[[pop]]$freq[filter_all]
            pop_info[[pop]]$lth <<- lth
            return()
        }
    }
    
    if (! pop %in% names(pop_def))
    {
        data_type <- pop_info[[pop]]$data_type
        if (data_type == "frac2")
        {
            pop_num <- paste0(pop, "~num")
            pop_deno <- paste0(pop, "~deno")
            
            if (keep_frac)
            {
                pop_info[[pop]]$num <<- columns[[pop_num]]
                pop_info[[pop]]$deno <<- columns[[pop_deno]]
            }
            pop_info[[pop]]$freq <<- ( columns[[pop_num]] / columns[[pop_deno]] ) %>% replace(is.nan(.), NA)
            pop_info[[pop]]$lth <<- lth
        } else if (data_type == "frac1")
        {
            geno <- columns[[pop]]
            geno_na <- is.na(geno)
            
            if (keep_frac)
            {
                pop_info[[pop]]$num <<- replace(geno, geno_na, 0)
                pop_info[[pop]]$deno <<- ifelse(geno_na, 0, pop_info[[pop]]$total)
            }
            pop_info[[pop]]$freq <<- geno / pop_info[[pop]]$total
            pop_info[[pop]]$lth <<- lth
        } else if (data_type == "freq")
        {
            # message("Can not be added!")
            pop_info[[pop]]$freq <<- columns[[pop]]
            pop_info[[pop]]$lth <<- lth
        }
    } else    # pop %in% names(pop_def)
    {
        pop_info[[pop]]$num <<- 0
        pop_info[[pop]]$deno <<- 0
        pop_info[[pop]]$total <<- 0
        for (pop0 in pop_def[[pop]])
        {
            get_freq(pop0, lth, columns, keep_frac = TRUE, pop_info)
            pop_info[[pop]]$num <<- pop_info[[pop]]$num + pop_info[[pop0]]$num
            pop_info[[pop]]$deno <<- pop_info[[pop]]$deno + pop_info[[pop0]]$deno
            pop_info[[pop]]$total <<- pop_info[[pop]]$total + pop_info[[pop0]]$total
        }
        pop_info[[pop]]$freq <<- (pop_info[[pop]]$num / pop_info[[pop]]$deno) %>% replace(is.nan(.), NA)
        pop_info[[pop]]$lth <<- lth
    }
    
}

column_name_outgroup <- asc_outgroup %>% lapply(pop_to_column, pop_def) %>% unlist() %>% base::unique()
column_outgroup <- pop_to_freq(column_name_outgroup, frequency_files_info_filtered)
get_freq(asc_outgroup, lth = "full", column_outgroup, keep_frac = FALSE, pop_info)
rm(column_outgroup)
filter_all <- pop_info[[asc_outgroup]]$freq <= outgroup_max_hetero | pop_info[[asc_outgroup]]$freq >= 1 - outgroup_max_hetero # filter_outgroup_max_hetero only


if (is_ascertained)
{
    column_name_refpop <- asc_refpop %>% lapply(pop_to_column, pop_def) %>% unlist() %>% base::unique()
    column_refpop <- pop_to_freq(column_name_refpop, frequency_files_info_filtered)
    get_freq(asc_refpop, lth = "full", column_refpop, keep_frac = TRUE, pop_info)
    rm(column_refpop)
    filter_refpop_max_miss <- pop_info[[asc_refpop]]$deno >= pop_info[[asc_refpop]]$total * (1 - refpop_max_miss)

    filter_all <- (filter_all & filter_refpop_max_miss) %>% {if (refpop_remove_homo) . & (!(pop_info[[asc_refpop]]$freq == 0 | pop_info[[asc_refpop]]$freq == 1)) else .} %>% replace(is.na(.), FALSE)

}

# the position file is to define a subset
if (!is.na(position_file) && !file.exists(position_file))
{
    position_file <- NA
    message("The position file", position_file, "does not exist.")
}

if (is.na(position_file))
{
    position <- pop_to_freq(c("CHR", "POS"), frequency_files_info_filtered, mask = filter_all)
} else
{
    position_all <- pop_to_freq(c("CHR", "POS"), frequency_files_info_filtered)
    position_sub <- utils::read.table(position_file, header = FALSE, sep = "\t", col.names = c("CHR", "POS"))
    filter_all <- filter_all & ( (position_all[["POS"]] * 100 + position_all[["CHR"]]) %in% (position_sub[["POS"]] * 100 + position_sub[["CHR"]]) )
    position <- list("CHR" = position_all[["CHR"]][filter_all], "POS" = position_all[["POS"]][filter_all])
}

pop_all <- ( if (is_ascertained) c(asc_refpop, asc_outgroup, pop_left1, pop_right1) else c(asc_outgroup, pop_left1, pop_right1) ) %>% base::unique()

column_name_test_related <- pop_all %>% lapply(pop_to_column, pop_def) %>% unlist() %>% base::unique()

column_test_related <- pop_to_freq(column_name_test_related, frequency_files_info_filtered, mask = filter_all)

lapply(pop_all, get_freq, lth = "asc", columns = column_test_related, keep_frac = FALSE, pop_info=pop_info)

rm(column_test_related)


fill_na_with_ref <- function(tgt, ref)
{
    pos <- is.na(tgt)
    tgt[pos] <- if (length(ref)==1) ref else ref[pos]
    return(tgt)
}

outgroup_freq <- pop_info[[asc_outgroup]]$freq %>% {if (outgroup_round_freq) round(.) else .}
if (is_ascertained)
{

refpop_freq <- pop_info[[asc_refpop]]$freq

derived_freq <- ifelse(outgroup_freq < 0.5, refpop_freq, 1 - refpop_freq)

if (file.exists(asc_freq_table))
    asc_conds <- utils::read.table(asc_freq_table, header = FALSE, sep = "\t", col.names = c("min_freq", "max_freq"))

asc_freq_unique <- c(asc_conds[["min_freq"]], asc_conds[["max_freq"]]) %>% base::unique() %>% base::sort()
asc_freq_unique_with_nega <- c(-1, asc_freq_unique, 2)

derived_freq_level <- cut(derived_freq, breaks = asc_freq_unique_with_nega, labels = 0:(length(asc_freq_unique_with_nega)-2)) %>% as.numeric()

asc_conds$min_freq_level <- cut(asc_conds$min_freq, breaks = asc_freq_unique_with_nega, labels = 0:(length(asc_freq_unique_with_nega)-2)) %>% as.numeric()
asc_conds$max_freq_level <- cut(asc_conds$max_freq, breaks = asc_freq_unique_with_nega, labels = 0:(length(asc_freq_unique_with_nega)-2)) %>% as.numeric()
}

all_stats <- tidyr::crossing(asc_outgroup, pop_left1, pop_right1)

block_id <- if (jackknife=="CHR") position[["CHR"]] else position[["CHR"]] * 100 + position[["POS"]] %/% jackknife

block_id_unique <- block_id %>% base::unique()


get_blockfsum <- function(asc_outgroup, pop_left1, pop_right1)
{
    asc_outgroup_freq <- pop_info[[asc_outgroup]]$freq
    pop_left1_freq <- pop_info[[pop_left1]]$freq %>% {if (!is.na(pop_left1_fillna)) fill_na_with_ref(., ref=if (pop_left1_fillna=="Outgroup") outgroup_freq else 0) else .}
    pop_right1_freq <- pop_info[[pop_right1]]$freq %>% {if (!is.na(pop_right1_fillna)) fill_na_with_ref(., ref=if (pop_right1_fillna=="Outgroup") outgroup_freq else 0) else .}

    f_value <- (pop_left1_freq - asc_outgroup_freq) * (pop_right1_freq - asc_outgroup_freq)

    if (!is.na(asc_freq_table))
    {
        fsum_block_asc_cutoff <- data.table::data.table(block_id=block_id, derived_freq_level=derived_freq_level, f_value=f_value)[, .(Nr_sites_all=.N, Nr_sites_nonmissing=sum(!is.na(f_value)), fsum=sum(f_value, na.rm=TRUE)), .(block_id, derived_freq_level)][base::order(block_id, derived_freq_level)][, .(derived_freq_level_max=derived_freq_level, Nr_sites_all=cumsum(Nr_sites_all), Nr_sites_nonmissing=cumsum(Nr_sites_nonmissing), fsum=cumsum(fsum)), .(block_id)]

        f_table <- data.table::data.table()
        f_value_block_table <- data.table::data.table()

        for (i in 1:nrow(asc_conds))
        {
            min_freq_level <- asc_conds$min_freq_level[i]
            max_freq_level <- asc_conds$max_freq_level[i]

            block_remain <- block_id_unique
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
            
            jackknife_weight_block <- if (all_sites_block) fsum_block$Nr_sites_all else fsum_block$Nr_sites_nonmissing
            f_value_block <- fsum_block$fsum / fsum_block$Nr_sites_nonmissing
            f_value_block_table <- f_value_block_table %>% rbind(data.table::data.table(min_freq = asc_conds$min_freq[i], max_freq = asc_conds$max_freq[i], Nr_sites = jackknife_weight_block, f = f_value_block, block_id = block_id_unique))


        }
    } else
    {
        fsum_block <- data.table::data.table(block_id=block_id, f_value=f_value)[, .(Nr_sites_all=.N, Nr_sites_nonmissing=sum(!is.na(f_value)), fsum=sum(f_value, na.rm=TRUE)), .(block_id)][base::order(block_id)]
        
        jackknife_weight_block <- if (all_sites_block) fsum_block$Nr_sites_all else fsum_block$Nr_sites_nonmissing
        f_value_block <- fsum_block$fsum / fsum_block$Nr_sites_nonmissing
        f_value_block_table <- data.table::data.table(Nr_sites = jackknife_weight_block, f = f_value_block, block_id = block_id_unique)

    }

    return(f_value_block_table[, c("asc_outgroup", "pop_left1", "pop_right1") := data.frame(asc_outgroup, pop_left1, pop_right1)])
}

f_block_table_all <- purrr::pmap(all_stats, get_blockfsum) %>% do.call(base::rbind, .)
return(f_block_table_all)
}
