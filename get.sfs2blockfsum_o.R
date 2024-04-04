library(dplyr)
library(optparse)

option_list <- list( 
    make_option(c("-f", "--sfsFile"), type = "character", help = "Mandatory. The directory of the freqSum/SFS file to process."),
    make_option(c("-b", "--blockCol"), type = "integer"),
    make_option(c("--blockID"), type = "character", default = NA),
    make_option(c("-n", "--npops"), type = "integer"),
    make_option(c("-p", "--popCols"), type = "character"),
    make_option(c("--fillHomoCols"), type = "character", default = NA, help = "Indicate the indexes of columns that should fill the missing sites with frequency 0. Default: NA."), # new
    make_option(c("-s", "--popSizes"), type = "character"),
    make_option(c("--countCol"), type = "integer", default = NA, help = "Indicate the index of column used for counting. If not indicated, count each row once."), # -c
    make_option(c("-a", "--ascCol"), type = "integer", default = NA, help = "Indicate the index of column used for ascertainment. If not indicated, compute without ascertainment."),
    make_option(c("-o", "--outFile"), type = "character", default = NA)
    )

opt <- parse_args(OptionParser(option_list=option_list))

########### Process the input arguments ##########

sfs_file <- opt$sfsFile

block_col <- opt$blockCol %>% as.numeric()

str2num <- function(num_str_each)
{
    num_str <- num_str_each %>% base::strsplit(split="*", fixed=TRUE) %>% unlist()
    NUM <- num_str[1] %>% base::strsplit(split=":") %>% unlist() %>% as.numeric()
    num_vec <- switch(length(NUM), "1"=NUM, "2"=seq(NUM[1], NUM[2]), "3"=seq(NUM[1], NUM[3], by=NUM[2]))
    return(if (length(num_str)==1) num_vec else rep(num_vec, each=as.numeric(num_str[2])))
}

block_id_all <- opt$blockID %>% strsplit(split=",") %>% unlist() %>% purrr::map(str2num) %>% unlist()

n_pops <- opt$npops %>% as.numeric()
if (!(n_pops %in% c(2,3,4)))
    stop("The number of populations (--npops) should be 2, 3 or 4.")

pop_cols <- opt$popCols %>% strsplit(split=",") %>% unlist() %>% as.numeric()
if (length(pop_cols) != n_pops)
    stop(paste0("The length of --popCols does not equal to --npops(", n_pops, ")."))

fill_homo_cols <- if (is.na(opt$fillHomoCols)) NA else opt$fillHomoCols %>% strsplit(split=",") %>% unlist() %>% as.numeric()

pop_sizes <- opt$popSizes %>% strsplit(split=",") %>% unlist() %>% as.numeric()
if (length(pop_sizes) != n_pops)
    stop(paste0("The length of --popSizes does not equal to --npops(", n_pops, ")."))

count_col <- opt$countCol %>% as.numeric() # if opt$countCol is NA, it will also be converted to NA.

asc_col <- opt$ascCol %>% as.numeric()
#if (is.na(asc_col))
#    stop("The argument --ascCol is mandatory.")

out_file <- opt$outFile

##################


sfs <- utils::read.table(sfs_file, header = FALSE, sep = "\t", na.strings = c("-1", "*"))

if (!is.na(fill_homo_cols))
    sfs[is.na(sfs[, fill_homo_cols]), fill_homo_cols] <- 0

allele_freq <- list()
for (i in 1:n_pops)
   allele_freq[[i]] <- if (pop_cols[i]) sfs[, pop_cols[i]] / pop_sizes[i] else 0

# Compute f value (f2/f3/f4). For sites that are not available in all pops, the value would be NA.
sfs$f <- switch(as.character(n_pops), "2"=(allele_freq[[2]] - allele_freq[[1]]) ^ 2, "3"=(allele_freq[[2]] - allele_freq[[1]]) * (allele_freq[[3]] - allele_freq[[1]]), "4"=(allele_freq[[2]] - allele_freq[[1]]) * (allele_freq[[4]] - allele_freq[[3]]))

fsum_block_asc_cutoff <- data.frame()

if (!is.na(asc_col))
{
    asc_freq_all <- base::unique(sfs[, asc_col]) %>% base::sort()

    #cumulative sum for cutoff min_freq_min:asc_freq_max, with asc_freq_max ranging from min_freq_min to max_freq_max
    for (block_id in block_id_all)
    {
        sfs_block <- sfs[sfs[, block_col]==block_id, ]

        Nr_sites_all <- 0
        Nr_sites_nonmissing <- 0
        fsum <- 0

        for (asc_freq_max in asc_freq_all)
        {
            sfs_block_asc <- sfs_block[sfs_block[, asc_col]==asc_freq_max, ]

            if (!is.na(count_col))
            {
                Nr_sites_all <- Nr_sites_all + sum(sfs_block_asc[, count_col])
                Nr_sites_nonmissing <- Nr_sites_nonmissing + sum(sfs_block_asc[!is.na(sfs_block_asc$f), count_col])
                fsum <- fsum + sum(sfs_block_asc[, count_col] * sfs_block_asc$f, na.rm = TRUE)
            } else
            {
                Nr_sites_all <- Nr_sites_all + nrow(sfs_block_asc)
                Nr_sites_nonmissing <- Nr_sites_nonmissing + sum(!is.na(sfs_block_asc$f))
                fsum <- fsum + sum(sfs_block_asc$f, na.rm = TRUE)
            }

            fsum_block_asc_cutoff <- base::rbind(fsum_block_asc_cutoff, list(block_id=block_id, asc_freq_max=asc_freq_max, Nr_sites_all=Nr_sites_all, Nr_sites_nonmissing=Nr_sites_nonmissing, fsum=fsum))
        }
    }
} else
{
    # is.na(asc_col)
    for (block_id in block_id_all)
    {
        sfs_block <- sfs[sfs[, block_col]==block_id, ]

        if (!is.na(count_col))
        {
            Nr_sites_all <- sum(sfs_block[, count_col])
            Nr_sites_nonmissing <- sum(sfs_block[!is.na(sfs_block_asc$f), count_col])
            fsum <- sum(sfs_block[, count_col] * sfs_block$f, na.rm = TRUE)
        } else
        {
            Nr_sites_all <- nrow(sfs_block)
            Nr_sites_nonmissing <- sum(!is.na(sfs_block$f))
            fsum <- sum(sfs_block$f, na.rm = TRUE)
        }

        fsum_block_asc_cutoff <- base::rbind(fsum_block_asc_cutoff, list(block_id=block_id, Nr_sites_all=Nr_sites_all, Nr_sites_nonmissing=Nr_sites_nonmissing, fsum=fsum))
    }
}


if (!is.na(out_file))
{
    write.table(fsum_block_asc_cutoff, out_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else
    print(fsum_block_asc_cutoff)
