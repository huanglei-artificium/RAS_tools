
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
# multiple frequency files
get_frequency_files_info <- function(frequency_files_list)
{
    if (!file.exists(frequency_files_list))
    {
        message("Frequency files list", frequency_files_list, "does not exist.")
        stop()
    }
    frequency_files <- utils::read.table(frequency_files_list, header = FALSE, sep = "\t", col.names = c("frequency_file", "header_info"), fill = TRUE, na.strings = c("NA", ""))
    frequency_files_info <- purrr::pmap(frequency_files, get_frequency_file_info) %>% do.call(base::rbind, .)

    column_names_duplicated <- frequency_files_info$column_name %>% .[base::duplicated(.)] %>% base::unique()
    if (base::length(column_names_duplicated))
    {
        message("The following column(s) in frequency file(s) is/are duplicated: ", base::paste(column_names_duplicated, collapse=", "), ".")
        stop()
    }
    return(frequency_files_info)
}

# get_pop_def currently not here
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

pop_to_column <- function(pop, pop_def)
{
    if (! pop %in% names(pop_def))
    {
        return(if (pop_info[[pop]]$data_type == "frac2") paste0(pop, c("~num", "~deno")) else pop)
    } else
    {
        return(pop_def[[pop]] %>% lapply(pop_to_column, pop_def) %>% unlist() %>% base::unique())
    }
}


pop_to_freq <- function(pop_all, frequency_files_info_filtered, mask = NULL)
{
    pop_all_file_info <- frequency_files_info_filtered %>% dplyr::filter(column_name %in% pop_all)

    file_all <- pop_all_file_info$frequency_file %>% base::unique()

    geno_all <- list()
    for (file in file_all)
    {
        pop_file_info <- pop_all_file_info %>% filter(frequency_file == file)
        field <- pop_file_info$field
        pop <- pop_file_info$column_name

        is_geno <- tools::file_ext(file) == "geno"
        geno <- if (is_geno) readr::read_fwf(file, col_positions=readr::fwf_positions(field, field, col_names=pop), col_types=readr::cols(.default="i"), na="9") else data.table::fread(file=file, header=FALSE, select=field, na.strings=c("-1", "*"), col.names=pop)
        if (!is.null(mask))
            geno <- geno[mask, ]
        geno_all <- geno_all %>% append(geno)
    }

    return(geno_all)
}