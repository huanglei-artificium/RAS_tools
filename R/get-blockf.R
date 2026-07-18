#' Compute block-based rare allele sharing statistics
#'
#' Reads population frequency/genotype inputs, applies optional site filters and
#' ascertainment conditions, then returns block-level rare allele sharing
#' statistics.
#'
#' @param frequency_files_list Path to a tab-delimited frequency-file manifest.
#' @param base_dir Optional base directory for relative file paths. When set,
#'   relative paths passed to this function and relative paths in the manifest
#'   are resolved against `base_dir`. When omitted, manifest entries are resolved
#'   against the manifest file's directory.
#' @param pop_left1,pop_right1 Comma-separated population names, or paths to
#'   files with one population per line.
#' @param position_file Optional two-column `CHR`, `POS` site subset file.
#' @param pop_def_file Optional custom population definition file.
#' @param pop_left2,pop_right2 Optional secondary population sets for RASD
#'   contrasts.
#' @param pop_left_fillna,pop_right_fillna Optional missing-data replacement.
#'   Use `"Outgroup"` to replace missing values with outgroup frequencies, or
#'   another non-`NA` value to replace with zero.
#' @param jackknife Either `"CHR"` or an integer genomic window length.
#' @param is_ascertained Whether to compute ascertained frequency-bin results.
#' @param asc_refpop Reference population for derived-frequency ascertainment.
#' @param asc_freq_table Frequency-bin table with `min_freq` and `max_freq`.
#' @param all_sites_block Retained for compatibility; summarization uses this
#'   option in [get_f_from_fblock()].
#' @param refpop_max_miss Maximum tolerated missingness in `asc_refpop`.
#' @param refpop_remove_homo Whether to exclude homozygous reference sites.
#' @param asc_outgroup Outgroup population name.
#' @param outgroup_max_hetero Maximum tolerated outgroup heterozygosity.
#' @param outgroup_round_freq Whether to round outgroup frequency to 0 or 1.
#'
#' @return A data frame/data table of block-level statistics.
#' @export
get_blockf <- function(
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
) {
    full_info <- base::new.env(parent = emptyenv())

    if (!is.null(base_dir) && !is.na(base_dir)) {
        base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
    }

    frequency_files_list <- resolve_input_path(frequency_files_list, base_dir = base_dir)
    position_file <- resolve_input_path(position_file, base_dir = base_dir)
    pop_def_file <- resolve_input_path(pop_def_file, base_dir = base_dir)
    asc_freq_table <- resolve_input_path(asc_freq_table, base_dir = base_dir)

    get_frequency_files_info(frequency_files_list, full_info, base_dir = base_dir)

    full_info$defs <- get_pop_def(pop_def_file, full_info$pop_original)
    full_info$pop_newdefined <- names(full_info$defs)
    full_info$pop_input <- c(full_info$pop_original, full_info$pop_newdefined)

    pop_str_convert <- function(pop_str, universe) {
        if (length(pop_str) != 1 || is.na(pop_str)) {
            return(NULL)
        }
        pop <- if (file.exists(pop_str)) {
            utils::read.table(pop_str, header = FALSE, sep = "\t")[[1]]
        } else {
            base::strsplit(pop_str, ",")[[1]]
        }
        pop <- pop %>% .[. %in% universe] %>% base::unique()
        message("Input ", length(pop), " valid population(s): ", base::paste(pop, collapse = ", "), ".")
        pop
    }

    message("pop_left1:")
    pop_left1 <- pop_str_convert(pop_left1, full_info$pop_input)
    message("pop_left2:")
    pop_left2 <- pop_str_convert(pop_left2, full_info$pop_input)
    pop_left <- c(pop_left1, pop_left2)

    message("pop_right1:")
    pop_right1 <- pop_str_convert(pop_right1, full_info$pop_input)
    message("pop_right2:")
    pop_right2 <- pop_str_convert(pop_right2, full_info$pop_input)
    pop_right <- c(pop_right1, pop_right2)

    if (!length(pop_left1) || !length(pop_right1)) {
        stop("Both pop_left1 and pop_right1 must contain at least one valid population.", call. = FALSE)
    }
    if (!asc_outgroup %in% full_info$pop_input) {
        stop("asc_outgroup is not a valid population: ", asc_outgroup, call. = FALSE)
    }
    if (is_ascertained && !asc_refpop %in% full_info$pop_input) {
        stop("asc_refpop is not a valid population: ", asc_refpop, call. = FALSE)
    }

    full_info$pops <- base::data.frame(
        Nr_columns = c(rep(2, length(full_info$pop_two_column)), rep(1, length(full_info$pop_one_column))),
        total = full_info$columns$pop[c(paste0(full_info$pop_two_column, "~deno"), full_info$pop_one_column), "denominator"],
        keep_frac = full_info$pop_original %in% base::unique(c(base::unlist(full_info$defs), asc_refpop)),
        row.names = full_info$pop_original
    )

    if (length(full_info$pop_newdefined)) {
        full_info$pops <- full_info$pops %>%
            rbind(base::data.frame(
                Nr_columns = 2,
                total = base::sapply(full_info$pop_newdefined, function(pop) {
                    sum(full_info$pops[full_info$defs[[pop]], "total"])
                }),
                keep_frac = TRUE,
                row.names = full_info$pop_newdefined
            ))
    }

    full_info$pops$lth <- NA
    full_info$freq <- list()
    full_info$frac <- list()

    get_freq(asc_outgroup, 1, NULL, full_info)

    filter_all <- (!is.na(full_info$freq[[asc_outgroup]])) %>%
        {if (!is.na(outgroup_max_hetero)) . & (full_info$freq[[asc_outgroup]] %>% {. <= outgroup_max_hetero | . >= 1 - outgroup_max_hetero}) else .}

    if (is_ascertained) {
        get_freq(asc_refpop, 1, NULL, full_info)
        filter_refpop_max_miss <- full_info$frac[[asc_refpop]]$deno >=
            full_info$pops[asc_refpop, "total"] * (1 - refpop_max_miss)
        filter_all <- (filter_all & filter_refpop_max_miss) %>%
            {if (refpop_remove_homo) . & full_info$freq[[asc_refpop]] %>% {. > 0 & . < 1} else .}
    }

    if (!is.na(position_file) && !file.exists(position_file)) {
        warning("Position file does not exist; ignoring: ", position_file, call. = FALSE)
        position_file <- NA
    }

    if (is.na(position_file)) {
        position <- data.table::fread(
            file = full_info$columns$posi$frequency_file[1],
            header = FALSE,
            select = c(1, 2),
            na.strings = c("-1", "*"),
            col.names = c("CHR", "POS")
        )
        position <- position[filter_all, ]
    } else {
        position_all <- data.table::fread(
            file = full_info$columns$posi$frequency_file[1],
            header = FALSE,
            select = c(1, 2),
            na.strings = c("-1", "*"),
            col.names = c("CHR", "POS")
        )
        position_sub <- utils::read.table(position_file, header = FALSE, sep = "\t", col.names = c("CHR", "POS"))
        filter_all <- filter_all &
            ((position_all[["POS"]] * 100 + position_all[["CHR"]]) %in%
                (position_sub[["POS"]] * 100 + position_sub[["CHR"]]))
        position <- list("CHR" = position_all[["CHR"]][filter_all], "POS" = position_all[["POS"]][filter_all])
    }

    pop_all <- c(asc_outgroup, pop_left, pop_right) %>%
        {if (is_ascertained) c(., asc_refpop) else .} %>%
        base::unique()
    get_freq(pop_all, 2, filter_all, full_info)

    full_info$outgroup_freq <- full_info$freq[[asc_outgroup]] %>%
        {if (outgroup_round_freq) base::round(.) else .}

    if (is_ascertained) {
        refpop_freq <- full_info$freq[[asc_refpop]]
        derived_freq <- ifelse(full_info$outgroup_freq < 0.5, refpop_freq, 1 - refpop_freq)

        if (!base::file.exists(asc_freq_table)) {
            stop("asc_freq_table does not exist: ", asc_freq_table, call. = FALSE)
        }
        full_info$asc_conds <- utils::read.table(
            asc_freq_table,
            header = FALSE,
            sep = "\t",
            col.names = c("min_freq", "max_freq")
        )

        asc_freq_unique <- c(full_info$asc_conds[["min_freq"]], full_info$asc_conds[["max_freq"]]) %>%
            base::unique() %>%
            base::sort()
        asc_freq_unique_with_nega <- c(-1, asc_freq_unique, 2)

        full_info$derived_freq_level <- cut(
            derived_freq,
            breaks = asc_freq_unique_with_nega,
            labels = 0:(length(asc_freq_unique_with_nega) - 2)
        ) %>%
            as.numeric()
        full_info$asc_conds$min_freq_level <- cut(
            full_info$asc_conds$min_freq,
            breaks = asc_freq_unique_with_nega,
            labels = 0:(length(asc_freq_unique_with_nega) - 2)
        ) %>%
            as.numeric()
        full_info$asc_conds$max_freq_level <- cut(
            full_info$asc_conds$max_freq,
            breaks = asc_freq_unique_with_nega,
            labels = 0:(length(asc_freq_unique_with_nega) - 2)
        ) %>%
            as.numeric()

        full_info$all_stats_RAS <- tidyr::crossing(asc_outgroup, pop_left, pop_right) %>% as.data.frame()
        row.names(full_info$all_stats_RAS) <- seq_len(nrow(full_info$all_stats_RAS))
    }

    full_info$block_id <- if (jackknife == "CHR") {
        position[["CHR"]]
    } else {
        position[["CHR"]] * 100 + position[["POS"]] %/% jackknife
    }
    full_info$block_id_unique <- full_info$block_id %>% base::unique()

    pop_left2_out <- if (is.null(pop_left2)) NA_character_ else pop_left2
    pop_right2_out <- if (is.null(pop_right2)) NA_character_ else pop_right2
    full_info$all_stats_RASD <- tidyr::crossing(
        asc_outgroup,
        pop_left1,
        pop_left2 = pop_left2_out,
        pop_right1,
        pop_right2 = pop_right2_out
    ) %>%
        as.data.frame()
    row.names(full_info$all_stats_RASD) <- seq_len(nrow(full_info$all_stats_RASD))

    full_info$pop_left_fillna <- pop_left_fillna
    full_info$pop_right_fillna <- pop_right_fillna

    if (is_ascertained) {
        full_info$RAS_table_all <- purrr::pmap(full_info$all_stats_RAS, get_blockRAS_asc, full_info = full_info)
        full_info$RASD_table_all <- purrr::pmap(
            full_info$all_stats_RASD,
            function(asc_outgroup, pop_left1, pop_left2, pop_right1, pop_right2) {
                get_RASD_from_RAS(
                    asc_outgroup = asc_outgroup,
                    pop_left1 = pop_left1,
                    pop_left2 = if (is.na(pop_left2)) NULL else pop_left2,
                    pop_right1 = pop_right1,
                    pop_right2 = if (is.na(pop_right2)) NULL else pop_right2,
                    full_info = full_info
                )
            }
        )
    } else {
        full_info$RASD_table_all <- purrr::pmap(
            full_info$all_stats_RASD,
            function(asc_outgroup, pop_left1, pop_left2, pop_right1, pop_right2) {
                get_blockf_nonasc(
                    asc_outgroup = asc_outgroup,
                    pop_left1 = pop_left1,
                    pop_left2 = if (is.na(pop_left2)) NULL else pop_left2,
                    pop_right1 = pop_right1,
                    pop_right2 = if (is.na(pop_right2)) NULL else pop_right2,
                    full_info = full_info
                )
            }
        )
    }

    base::mapply(
        base::cbind,
        full_info$RASD_table_all,
        full_info$all_stats_RASD %>% base::split(seq_len(nrow(.))),
        SIMPLIFY = FALSE
    ) %>%
        dplyr::bind_rows()
}
