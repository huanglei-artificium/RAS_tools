library(dplyr)
library(optparse)

option_list <- list(
    make_option(c("--positionFile"), type = "character", default = NA),
    make_option(c("--frequencyFilesList"), type = "character"),
    make_option(c("-n", "--npops"), type = "integer"),
    make_option(c("--popA"), type = "character"),
    make_option(c("--popAfillna"), type = "character", default = NA),
    make_option(c("--popB"), type = "character"),
    make_option(c("--popBfillna"), type = "character", default = NA),
    make_option(c("--popC"), type = "character"),
    make_option(c("--popCfillna"), type = "character", default = NA),
    make_option(c("--popD"), type = "character"),
    make_option(c("--popDfillna"), type = "character", default = NA),
    make_option(c("--jackknife"), type = "character"),
    make_option(c("--ascFreqFile"), type = "character", default = NA),
    make_option(c("--outgroup"), type = "character", default = NA),
    make_option(c("-o", "--outFile"), type = "character", default = NA)
    )

opt <- parse_args(OptionParser(option_list=option_list))

########### Process the input arguments ##########

position_file <- opt$positionFile

frequency_files_list <- opt$frequencyFilesList

str2num <- function(num_str_each)
{
    num_str <- num_str_each %>% base::strsplit(split="*", fixed=TRUE) %>% unlist()
    NUM <- num_str[1] %>% base::strsplit(split=":") %>% unlist() %>% as.numeric()
    num_vec <- switch(length(NUM), "1"=NUM, "2"=seq(NUM[1], NUM[2]), "3"=seq(NUM[1], NUM[3], by=NUM[2]))
    return(if (length(num_str)==1) num_vec else rep(num_vec, each=as.numeric(num_str[2])))
}

n_pops <- opt$npops %>% as.numeric()
if (!(n_pops %in% c(2,3,4)))
    stop("The number of populations (--npops) should be 2, 3 or 4.")

popA_str <- opt$popA

popA_fillna <- opt$popAfillna

popB_str <- opt$popB

popB_fillna <- opt$popBfillna

popC_str <- opt$popC

popC_fillna <- opt$popCfillna

popD_str <- opt$popD

popD_fillna <- opt$popDfillna

jackknife <- opt$jackknife

asc_freq_file <- opt$ascFreqFile

outgroup <- opt$outgroup

out_file <- opt$outFile

##################

get_frequency_file_info <- function(frequency_file, header_file)
{
    is_geno <- tools::file_ext(frequency_file) == "geno"
    if (is.na(header_file))
        header_file <- paste0(tools::file_path_sans_ext(frequency_file), if (is_geno) ".ind" else ".header")
    
    header <- utils::read.table(header_file, header = FALSE, sep = "\t", fill = TRUE)
    return(data.frame(frequency_file = frequency_file, field = 1:nrow(header), pop_name = header[[1]], denominator = if (is_geno) 2 else header[[2]]))
}
frequency_files <- utils::read.table(frequency_files_list, header = FALSE, sep = "\t", col.names = c("frequency_file", "header_file"), fill = TRUE)
frequency_files_info <- purrr::pmap(frequency_files, get_frequency_file_info) %>% do.call(base::rbind, .)

pop_names_all <- frequency_files_info$pop_name %>% .[! . %in% c("CHR", "POS", "ALT", "REF")]

pop_names_duplicated <- pop_names_all %>% .[base::duplicated(.)] %>% base::unique()
if (base::length(pop_names_duplicated))
{
    message("The following population name(s) is/are duplicated: ", base::paste(pop_names_duplicated, collapse=", "), ".")
    stop()
}

# should modify outgroup case

if (!is.na(asc_freq_file))
{
    if (file.exists(asc_freq_file))
    {
        asc_freq <- utils::read.table(asc_freq_file, header = FALSE, sep = "\t")[[1]]
    } else
    {
        asc_freq_file <- NA
        message("The ascertainment file ", asc_freq_file , " does not exist. No ascertainment is defined.")
    }
}

if (is.na(asc_freq_file) && !is.na(outgroup))
{
    message("Outgroup only applies when ascertainment is defined.")
    outgroup <- NA
}

if (!is.na(outgroup) && ! outgroup %in% pop_names_all)
{
    message("Outgroup ", outgroup, " is not in the population name list.")
    outgroup <- NA
}

pop_str_convert <- function(pop_str, universe)
{
    if (pop_str == "Outgroup")
    {
        pop <- if (!is.na(outgroup)) outgroup else NULL
    } else
        pop <- if (file.exists(pop_str)) utils::read.table(pop_str, header=FALSE, sep = "\t", fill = TRUE)[[1]] else base::strsplit(pop_str, ",")[[1]] %>% .[. %in% universe] %>% base::unique()
    message("Input ", length(pop), " valid population(s): ", base::paste(pop, collapse=", "), ".")
    return(pop)
}

message("PopA:")
popA_all <- pop_str_convert(popA_str, pop_names_all)
message("PopB:")
popB_all <- pop_str_convert(popB_str, pop_names_all)
if (n_pops >= 3)
{
    message("PopC:")
    popC_all <- pop_str_convert(popC_str, pop_names_all)
} else
    popC_all <- NA
if (n_pops == 4)
{
    message("PopD:")
    popD_all <- pop_str_convert(popD_str, pop_names_all)
} else
    popD_all <- NA


##################

fill_na_with_ref <- function(tgt, ref)
{
    pos <- is.na(tgt)
    tgt[pos] <- if (length(ref)==1) ref else ref[pos]
    return(tgt)
}

pop_to_freq <- function(pop_all, frequency_files_info)
{
    pop_all_file_info <- frequency_files_info %>% dplyr::filter(pop_name %in% pop_all)

    file_all <- pop_all_file_info$frequency_file %>% base::unique()

    geno_all <- list()
    for (file in file_all)
    {
        pop_file_info <- pop_all_file_info %>% filter(frequency_file == file)
        field <- pop_file_info$field
        pop <- pop_file_info$pop_name
        denominator <- pop_file_info$denominator %>% replace(is.na(.), 1)

        is_geno <- tools::file_ext(file) == "geno"
        geno_all <- geno_all %>% append( (if (is_geno) readr::read_fwf(file, col_positions=readr::fwf_positions(field, field, col_names=pop), col_types=readr::cols(.default="i"), na="9") else data.table::fread(file=file, header=FALSE, select=field, na.strings=c("-1", "*"), col.names=pop) ) %>% purrr::map2(denominator, .f= ~ .x / .y) )
    }

    return(geno_all)
}

if (!is.na(position_file) && !file.exists(position_file))
{
    position_file <- NA
    message("The position file", position_file, "does not exist.")
}
    
position <- if (!is.na(position_file)) utils::read.table(position_file, header = FALSE, sep = "\t", col.names = c("CHR", "POS")) else pop_to_freq(c("CHR", "POS"), frequency_files_info)

outgroup_freq <- pop_to_freq(outgroup, frequency_files_info)[[outgroup]]

popA_freq <- pop_to_freq(popA_all, frequency_files_info) %>% {if (!is.na(popA_fillna)) base::lapply(., fill_na_with_ref, ref=if (popA_fillna=="Outgroup") outgroup_freq else 0) else .}

popB_freq <- pop_to_freq(popB_all, frequency_files_info) %>% {if (!is.na(popB_fillna)) base::lapply(., fill_na_with_ref, ref=if (popB_fillna=="Outgroup") outgroup_freq else 0) else .}

if (n_pops >= 3)
{
    popC_freq <- pop_to_freq(popC_all, frequency_files_info) %>% {if (!is.na(popC_fillna)) base::lapply(., fill_na_with_ref, ref=if (popC_fillna=="Outgroup") outgroup_freq else 0) else .}
} else
    popC_freq <- NA

if (n_pops == 4)
{
    popD_freq <- pop_to_freq(popD_all, frequency_files_info) %>% {if (!is.na(popD_fillna)) base::lapply(., fill_na_with_ref, ref=if (popD_fillna=="Outgroup") outgroup_freq else 0) else .}
} else
    popD_freq <- NA


all_stats <- tidyr::crossing(popA = popA_all, popB = popB_all, popC = popC_all, popD = popD_all)

block_id <- if (jackknife=="CHR") position[["CHR"]] else position[["CHR"]] * 100 + position[["POS"]] %/% jackknife

get_blockfsum <- function(popA, popB, popC, popD)
{
    f_value <- switch(as.character(n_pops), "2"=(popB_freq[[popB]] - popA_freq[[popA]]) ^ 2, "3"=(popB_freq[[popB]] - popA_freq[[popA]]) * (popC_freq[[popC]] - popA_freq[[popA]]), "4"=(popB_freq[[popB]] - popA_freq[[popA]]) * (popD_freq[[popD]] - popC_freq[[popC]]))

    #return(f_value)

    if (!is.na(asc_freq_file))
    {
        fsum_block_asc_cutoff <- data.table::data.table(block_id=block_id, asc_freq=asc_freq, f_value=f_value)[, .(Nr_sites_all=.N, Nr_sites_nonmissing=sum(!is.na(f_value)), fsum=sum(f_value, na.rm=TRUE)), .(block_id, asc_freq)][base::order(block_id, asc_freq)][, .(asc_freq_max=asc_freq, Nr_sites_all=cumsum(Nr_sites_all), Nr_sites_nonmissing=cumsum(Nr_sites_nonmissing), fsum=cumsum(fsum)), .(block_id)]
    } else
        fsum_block_asc_cutoff <- data.table::data.table(block_id=block_id, f_value=f_value)[, .(Nr_sites_all=.N, Nr_sites_nonmissing=sum(!is.na(f_value)), fsum=sum(f_value, na.rm=TRUE)), .(block_id)][base::order(block_id)]
    
    return(fsum_block_asc_cutoff[, c("popA", "popB", "popC", "popD") := data.frame(popA, popB, popC, popD)])
}

fsum_block_asc_cutoff_all <- purrr::pmap(all_stats, get_blockfsum) %>% do.call(base::rbind, .)


if (!is.na(out_file))
{
    write.table(fsum_block_asc_cutoff_all, out_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else
    print(fsum_block_asc_cutoff_all)
