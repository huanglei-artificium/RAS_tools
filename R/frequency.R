get_frequency_file_info <- function(frequency_file, header_info = NA) {
    is_geno <- tools::file_ext(frequency_file) == "geno"
    if (base::is.na(header_info)) {
        header_info <- base::paste0(
            tools::file_path_sans_ext(frequency_file),
            if (is_geno) ".ind" else ".header"
        )
    }

    header <- (
        if (base::file.exists(header_info)) {
            utils::read.table(header_info, header = FALSE, sep = "\t")[[1]]
        } else {
            base::strsplit(header_info, ",")[[1]]
        }
    ) %>%
        base::strsplit(split = "[()]")

    data.frame(
        frequency_file = frequency_file,
        field = seq_along(header),
        column_name = base::sapply(header, `[`, 1),
        denominator = if (is_geno) 2 else base::as.numeric(base::sapply(header, `[`, 2))
    )
}

is_two_column_pop <- function(pop, tbl) {
    frac2_name <- paste0(pop, c("~num", "~deno"))
    if (!all(frac2_name %in% rownames(tbl))) {
        return(FALSE)
    }

    deno <- tbl[frac2_name, "denominator"] %>% base::unique()
    if (base::length(deno) > 1 || is.na(deno[[1]])) {
        message("Population ", pop, " does not have a valid numerator/denominator total count.")
        return(FALSE)
    }
    TRUE
}

get_frequency_files_info <- function(frequency_files_list, full_info, base_dir = NULL) {
    if (!file.exists(frequency_files_list)) {
        stop("Frequency files list does not exist: ", frequency_files_list, call. = FALSE)
    }

    manifest_dir <- dirname(normalizePath(frequency_files_list, winslash = "/", mustWork = TRUE))
    path_base_dir <- if (is.null(base_dir) || is.na(base_dir)) {
        manifest_dir
    } else {
        normalizePath(base_dir, winslash = "/", mustWork = TRUE)
    }
    frequency_files <- utils::read.table(
        frequency_files_list,
        header = FALSE,
        sep = "\t",
        col.names = c("frequency_file", "header_info"),
        fill = TRUE,
        na.strings = c("NA", "")
    )

    frequency_files$frequency_file <- vapply(
        frequency_files$frequency_file,
        resolve_input_path,
        character(1),
        base_dir = path_base_dir
    )
    frequency_files$header_info <- vapply(
        frequency_files$header_info,
        resolve_input_path,
        character(1),
        base_dir = path_base_dir
    )

    frequency_files_info <- purrr::pmap(frequency_files, get_frequency_file_info) %>%
        do.call(base::rbind, .)

    column_names_duplicated <- frequency_files_info$column_name %>%
        .[base::duplicated(.)] %>%
        base::unique()
    if (base::length(column_names_duplicated)) {
        stop(
            "Duplicated columns in frequency file(s): ",
            base::paste(column_names_duplicated, collapse = ", "),
            call. = FALSE
        )
    }

    frequency_files_info <- frequency_files_info %>%
        dplyr::mutate(pop_name = .$column_name %>% base::strsplit(split = "~") %>% base::sapply(`[`, 1))
    row.names(frequency_files_info) <- frequency_files_info$column_name
    frequency_files_info$column_name <- NULL

    full_info$columns <- frequency_files_info %>%
        base::split(base::ifelse(.$pop_name %in% c("CHR", "POS", "ALT", "REF"), "posi", "pop"))
    full_info$columns$posi <- full_info$columns$posi %>%
        .[base::row.names(.) %in% c("CHR", "POS", "ALT", "REF"), ]

    pop_name_count <- base::table(full_info$columns$pop$pop_name)
    full_info$pop_two_column <- pop_name_count %>%
        .[. == 2] %>%
        base::names() %>%
        .[base::sapply(., is_two_column_pop, tbl = full_info$columns$pop)]
    full_info$pop_one_column <- pop_name_count %>%
        .[. == 1] %>%
        base::names() %>%
        .[. %in% base::row.names(full_info$columns$pop)]
    full_info$pop_original <- c(full_info$pop_two_column, full_info$pop_one_column)
    full_info$columns$pop <- full_info$columns$pop %>%
        .[.$pop_name %in% full_info$pop_original, ]
}

get_pop_def <- function(pop_def_file, pop_input) {
    if (is.na(pop_def_file)) {
        return(list())
    }
    if (!file.exists(pop_def_file)) {
        warning("Population definition file does not exist: ", pop_def_file, call. = FALSE)
        return(list())
    }

    pop_def_table_raw <- utils::read.table(pop_def_file, header = FALSE, sep = "\t")[[1]] %>%
        base::strsplit(split = ":")

    pop_def <- list()
    for (pop_def_item in pop_def_table_raw) {
        pop_exist <- c(pop_input, names(pop_def))
        if (pop_def_item[1] %in% pop_exist) {
            message("Population ", pop_def_item[1], " already exists; skipping custom definition.")
            next
        }
        pop_comp <- pop_def_item[2] %>%
            base::strsplit(split = ",") %>%
            .[[1]] %>%
            base::unique()
        pop_comp_noexist <- pop_comp %>% .[!. %in% pop_exist]
        if (length(pop_comp_noexist) > 0) {
            message("Population(s) not existing: ", base::paste(pop_comp_noexist, collapse = ", "), ". Not included.")
            pop_comp <- pop_comp %>% .[. %in% pop_exist]
        }
        pop_def[[pop_def_item[1]]] <- pop_comp
    }
    pop_def
}

get_freq <- function(pops, lth, mask, full_info) {
    pop_exist <- pops %>% .[. %in% names(full_info$freq)]
    pop_to_mask <- pop_exist %>% .[full_info$pops[., "lth"] < lth]

    for (pop in pop_to_mask) {
        full_info$freq[[pop]] <- full_info$freq[[pop]][mask]
        if (full_info$pops[pop, "keep_frac"]) {
            full_info$frac[[pop]] <- full_info$frac[[pop]][mask, ]
        }
    }

    pop_to_read <- pops %>% .[!. %in% pop_exist]
    pop_defined_to_read <- pop_to_read %>% .[. %in% full_info$pop_newdefined]
    pop_original_to_read <- c(
        pop_to_read %>% .[. %in% full_info$pop_original],
        full_info$defs[pop_defined_to_read] %>% base::unlist() %>% .[!. %in% pop_exist]
    ) %>%
        base::unique()

    pop_all_file_info <- full_info$columns$pop %>% .[.$pop_name %in% pop_original_to_read, ]
    file_all <- pop_all_file_info$frequency_file %>% base::unique()

    for (file in file_all) {
        pop_file_info <- pop_all_file_info %>% .[.$frequency_file == file, ]
        field <- pop_file_info$field
        column_name <- base::row.names(pop_file_info)

        is_geno <- tools::file_ext(file) == "geno"
        columns <- if (is_geno) {
            readr::read_fwf(
                file,
                col_positions = readr::fwf_positions(field, field, col_names = column_name),
                col_types = readr::cols(.default = "i"),
                na = "9"
            )
        } else {
            data.table::fread(
                file = file,
                header = FALSE,
                select = field,
                na.strings = c("-1", "*"),
                col.names = column_name
            )
        }
        if (!is.null(mask)) {
            columns <- columns[mask, ]
        }

        pop_name_uniq <- pop_file_info$pop_name %>% base::unique()
        for (pop in pop_name_uniq) {
            if (full_info$pops[pop, "Nr_columns"] == 2) {
                full_info$freq[[pop]] <- (columns[[paste0(pop, "~num")]] / columns[[paste0(pop, "~deno")]]) %>%
                    replace(is.nan(.), NA)
                if (full_info$pops[pop, "keep_frac"]) {
                    full_info$frac[[pop]] <- data.frame(
                        num = columns[[paste0(pop, "~num")]],
                        deno = columns[[paste0(pop, "~deno")]]
                    )
                }
            } else {
                full_info$freq[[pop]] <- columns[[pop]] / full_info$pops[pop, "total"]
                if (full_info$pops[pop, "keep_frac"]) {
                    full_info$frac[[pop]] <- data.frame(
                        num = dplyr::coalesce(columns[[pop]], 0),
                        deno = ifelse(is.na(columns[[pop]]), 0, full_info$pops[pop, "total"])
                    )
                }
            }
        }
    }

    for (pop in pop_defined_to_read) {
        full_info$frac[[pop]] <- base::Reduce(`+`, full_info$frac[full_info$defs[[pop]]])
        full_info$freq[[pop]] <- (full_info$frac[[pop]]$num / full_info$frac[[pop]]$deno) %>%
            replace(is.nan(.), NA)
    }

    full_info$pops[c(pop_to_mask, pop_original_to_read, pop_defined_to_read), "lth"] <- lth
}
