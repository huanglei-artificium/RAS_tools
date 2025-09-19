get_blockf <- function(
    position_file = NA,
    frequency_files_list,
    pop_def_file = NA,
    pop_left1,
    pop_left2 = NA,
    pop_left_fillna = NA,
    pop_right1,
    pop_right2 = NA,
    pop_right_fillna = NA,
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
    
full_info <- base::new.env()

get_frequency_files_info(frequency_files_list, full_info)

# for self-defined populations, process below.
full_info$defs <- get_pop_def(pop_def_file, full_info$pop_original)

# names(pop_def) == pop_newdefined
full_info$pop_newdefined <- names(full_info$defs)
full_info$pop_input <- c(full_info$pop_original, full_info$pop_newdefined)


pop_str_convert <- function(pop_str, universe)
{
    pop <- ( if (file.exists(pop_str)) utils::read.table(pop_str, header=FALSE, sep = "\t")[[1]] else base::strsplit(pop_str, ",")[[1]] ) %>% .[. %in% universe] %>% base::unique()
    message("Input ", length(pop), " valid population(s): ", base::paste(pop, collapse=", "), ".")
    return(pop)
}

message("pop_left1:")
pop_left1 <- pop_str_convert(pop_left1_str, full_info$pop_input)
message("pop_left2:")
pop_left2 <- if (!is.na(pop_left2_str)) pop_str_convert(pop_left2_str, full_info$pop_input) else NULL
pop_left <- c(pop_left1, pop_left2)

message("pop_right1:")
pop_right1 <- pop_str_convert(pop_right1_str, full_info$pop_input)
message("pop_right2:")
pop_right2 <- if (!is.na(pop_right2_str)) pop_str_convert(pop_right2_str, full_info$pop_input) else NULL
pop_right <- c(pop_right1, pop_right2)


full_info$pops <- base::data.frame(Nr_columns = c(rep(2, length(full_info$pop_two_column)), rep(1, length(full_info$pop_one_column))), total = full_info$columns$pop[c(paste0(full_info$pop_two_column, "~deno"), full_info$pop_one_column), "denominator"], keep_frac = full_info$pop_original %in% base::unique(c(base::unlist(full_info$defs), asc_refpop) ),  row.names = full_info$pop_original)

# also include self-defined populations in full_info$pops
# row.names(full_info$pops) == pop_input  so far!!
full_info$pops <- full_info$pops %>% rbind(base::data.frame(Nr_columns = 2, total = base::sapply(full_info$pop_newdefined, function(pops){sum(full_info$pops[full_info$defs[[pops]], "total"])}), keep_frac = TRUE, row.names = full_info$pop_newdefined) )

# initially no population is read.
# NA; 1; 2
full_info$pops$lth <- NA


full_info$freq <- list()
full_info$frac <- list()

get_freq(asc_outgroup, 1, NULL, full_info)

filter_all <- (!is.na(full_info$freq[[asc_outgroup]]) ) %>% {if (!is.na(outgroup_max_hetero)) . & (full_info$freq[[asc_outgroup]] %>% {. <= outgroup_max_hetero | . >= 1 - outgroup_max_hetero} ) else . }

if (is_ascertained)
{
    get_freq(asc_refpop, 1, NULL, full_info)
    
    filter_refpop_max_miss <- full_info$frac[[asc_refpop]]$deno >= full_info$pops[asc_refpop, "total"] * (1 - refpop_max_miss)

    # %>% is in higher priority than & |
    filter_all <- (filter_all & filter_refpop_max_miss) %>% {if (refpop_remove_homo) . & full_info$freq[[asc_refpop]] %>% {. > 0 & . < 1}  else .}

}


if (!is.na(position_file) && !file.exists(position_file))
{
    position_file <- NA
    message("The position file", position_file, "does not exist.")
}

if (is.na(position_file))
{
    position <- pop_to_freq(c("CHR", "POS"), full_info$columns$posi, mask = filter_all)
} else
{
    position_all <- pop_to_freq(c("CHR", "POS"), full_info$columns$posi)
    position_sub <- utils::read.table(position_file, header = FALSE, sep = "\t", col.names = c("CHR", "POS"))
    filter_all <- filter_all & ( (position_all[["POS"]] * 100 + position_all[["CHR"]]) %in% (position_sub[["POS"]] * 100 + position_sub[["CHR"]]) )
    position <- list("CHR" = position_all[["CHR"]][filter_all], "POS" = position_all[["POS"]][filter_all])
}

pop_all <- c(asc_outgroup, pop_left, pop_right) %>% {if (is_ascertained) c(., asc_refpop) else .} %>% base::unique()

get_freq(pop_all, 2, filter_all, full_info)


fill_na_with_ref <- function(tgt, ref)
{
    pos <- is.na(tgt)
    tgt[pos] <- if (length(ref)==1) ref else ref[pos]
    return(tgt)
}

full_info$outgroup_freq <- full_info$freq[[asc_outgroup]] %>% {if (outgroup_round_freq) base::round(.) else .}
if (is_ascertained)
{

    refpop_freq <- full_info$freq[[asc_refpop]]

    derived_freq <- ifelse(full_info$outgroup_freq < 0.5, refpop_freq, 1 - refpop_freq)

    if (base::file.exists(asc_freq_table))
        full_info$asc_conds <- utils::read.table(asc_freq_table, header = FALSE, sep = "\t", col.names = c("min_freq", "max_freq"))

    asc_freq_unique <- c(full_info$asc_conds[["min_freq"]], full_info$asc_conds[["max_freq"]]) %>% base::unique() %>% base::sort()
    asc_freq_unique_with_nega <- c(-1, asc_freq_unique, 2)

    full_info$derived_freq_level <- cut(derived_freq, breaks = asc_freq_unique_with_nega, labels = 0:(length(asc_freq_unique_with_nega)-2)) %>% as.numeric()

    full_info$asc_conds$min_freq_level <- cut(full_info$asc_conds$min_freq, breaks = asc_freq_unique_with_nega, labels = 0:(length(asc_freq_unique_with_nega)-2)) %>% as.numeric()
    full_info$asc_conds$max_freq_level <- cut(full_info$asc_conds$max_freq, breaks = asc_freq_unique_with_nega, labels = 0:(length(asc_freq_unique_with_nega)-2)) %>% as.numeric()
}

full_info$all_stats_RAS <- tidyr::crossing(asc_outgroup, pop_left, pop_right) %>% as.data.frame()
row.names(full_info$all_stats_RAS) <- 1:nrow(full_info$all_stats_RAS)

full_info$block_id <- if (jackknife=="CHR") position[["CHR"]] else position[["CHR"]] * 100 + position[["POS"]] %/% jackknife

full_info$block_id_unique <- full_info$block_id %>% base::unique()


full_info$pop_left_fillna <- pop_left_fillna
full_info$pop_right_fillna <- pop_right_fillna


full_info$RAS_table_all <- purrr::pmap(full_info$all_stats_RAS, get_blockfsum, full_info = full_info)

full_info$all_stats_RASD <- tidyr::crossing(asc_outgroup, pop_left1, pop_left2, pop_right1, pop_right2) %>% as.data.frame()
row.names(full_info$all_stats_RASD) <- 1:nrow(full_info$all_stats_RASD)

full_info$RASD_table_all <- purrr::pmap(full_info$all_stats_RASD, get_RASD_from_RAS, full_info = full_info)


f_block_table_all <- full_info$RASD_table_all %>% do.call(base::rbind, .)
return(f_block_table_all)
}
