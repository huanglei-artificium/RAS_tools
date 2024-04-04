library(dplyr)
library(optparse)

option_list <- list( 
    make_option(c("-f", "--blockFSumFile"), type = "character", help = "Mandatory. The directory of the blockFSum file to process."),
    make_option(c("--blockID"), type = "character", default = NA),
    make_option(c("--minFreq"), type = "character", default = NA),
    make_option(c("--maxFreq"), type = "character", default = NA),
    make_option(c("--allSitesBlock"), action = "store_true", default = FALSE), # -l
    make_option(c("--minSitesBlock"), type = "integer", default = NA),
    make_option(c("-o", "--outFile"), type = "character", default = NA)
    )

opt <- parse_args(OptionParser(option_list=option_list))

########### Process the input arguments ##########

blockfsum_file <- opt$blockFSumFile

str2num <- function(num_str_each)
{
    num_str <- num_str_each %>% base::strsplit(split="*", fixed=TRUE) %>% unlist()
    NUM <- num_str[1] %>% base::strsplit(split=":") %>% unlist() %>% as.numeric()
    num_vec <- switch(length(NUM), "1"=NUM, "2"=seq(NUM[1], NUM[2]), "3"=seq(NUM[1], NUM[3], by=NUM[2]))
    return(if (length(num_str)==1) num_vec else rep(num_vec, each=as.numeric(num_str[2])))
}

block_id_all <- opt$blockID %>% strsplit(split=",") %>% unlist() %>% purrr::map(str2num) %>% unlist()

if (is.na(opt$minFreq) & is.na(opt$maxFreq))
{
    is_ascertained <- FALSE
} else
{
    is_ascertained <- TRUE

    min_freq_all <- opt$minFreq %>% strsplit(split=",") %>% unlist() %>% purrr::map(str2num) %>% unlist()
    max_freq_all <- opt$maxFreq %>% strsplit(split=",") %>% unlist() %>% purrr::map(str2num) %>% unlist()
    len_min <- length(min_freq_all)
    len_max <- length(max_freq_all)
    if (len_min!=1 & len_max!=1 & len_min!=len_max)
        stop("The length of --minFreq and --maxFreq should be equal, or one of them should be 1.")

    min_freq_min <- min(min_freq_all)
    max_freq_max <- max(max_freq_all)
    asc_conds <- data.frame(min_freq=min_freq_all, max_freq=max_freq_all)
}

all_sites_block <- opt$allSitesBlock

min_sites_block <- opt$minSitesBlock %>% as.numeric()

out_file <- opt$outFile

##################

fsum_block_asc_cutoff <- utils::read.table(blockfsum_file, header = TRUE, sep = "\t")

if (is_ascertained)
{
    f_block_all <- data.frame()

    if (all_sites_block & !is.na(min_sites_block))
    {
        # only for all_sites_block=TRUE and !is.na(min_sites_block)
        # divide the cutoff so that in each cutoff there is similar proportion of non-missing SNPs among different frequency.

        cutoff_range <- list()
        asc_freq_all <- base::unique(fsum_block_asc_cutoff[, "asc_freq_max"]) %>% base::sort()

        for (asc_max in base::rev(asc_freq_all[asc_freq_all <= max_freq_max]))
        {
            #print(asc_max)
            cutoff_range[[as.character(asc_max)]] <- list()

            asc_min <- asc_max
            while (asc_min > min_freq_min)
            {
                Nr_sites_all <- fsum_block_asc_cutoff[fsum_block_asc_cutoff$asc_freq_max==asc_max, "Nr_sites_all"] - fsum_block_asc_cutoff[fsum_block_asc_cutoff$asc_freq_max==asc_min-1, "Nr_sites_all"]
                if (all(Nr_sites_all >= min_sites_block))
                    break
                asc_min <- asc_min - 1
            }
            if (asc_min > min_freq_min)
            {
                f_block <- fsum_block_asc_cutoff[fsum_block_asc_cutoff$asc_freq_max==asc_max, c("Nr_sites_all", "Nr_sites_nonmissing", "fsum")] - fsum_block_asc_cutoff[fsum_block_asc_cutoff$asc_freq_max==asc_min-1, c("Nr_sites_all", "Nr_sites_nonmissing", "fsum")]
            } else 
                f_block <- fsum_block_asc_cutoff[fsum_block_asc_cutoff$asc_freq_max==asc_max, c("Nr_sites_all", "Nr_sites_nonmissing", "fsum")]
            f_block <- f_block %>% mutate(f=fsum/Nr_sites_nonmissing)

            cutoff_range[[as.character(asc_max)]]$asc_min <- asc_min
            cutoff_range[[as.character(asc_max)]]$f_block <- f_block
        }


        for (i in 1:nrow(asc_conds))
        {
            min_freq <- asc_conds$min_freq[i]
            max_freq <- asc_conds$max_freq[i]
            
            asc_max <- max_freq
            fsum <- 0
            Nr_sites_all <- 0
            repeat
            {
                f_block_cutoff <- cutoff_range[[as.character(asc_max)]]$f_block
                fsum <- fsum + f_block_cutoff$Nr_sites_all * f_block_cutoff$f
                Nr_sites_all <- Nr_sites_all + f_block_cutoff$Nr_sites_all
                
                asc_max <- cutoff_range[[as.character(asc_max)]]$asc_min - 1
                if (asc_max < min_freq)
                    break
            }
            f_block_all <- f_block_all %>% rbind(data.frame(block_id=block_id_all, jackknife_weight=Nr_sites_all, f=fsum/Nr_sites_all, Asc_min=min_freq, Asc_max=max_freq))
        }
    } else
    {
        # all_sites_block & is.na(min_sites_block) | !all_sites_block

        for (i in 1:nrow(asc_conds))
        {
            min_freq <- asc_conds$min_freq[i]
            max_freq <- asc_conds$max_freq[i]
            
            if (min_freq > min_freq_min)
            {
                f_block <- fsum_block_asc_cutoff[fsum_block_asc_cutoff$asc_freq_max==max_freq, c("Nr_sites_all", "Nr_sites_nonmissing", "fsum")] - fsum_block_asc_cutoff[fsum_block_asc_cutoff$asc_freq_max==min_freq-1, c("Nr_sites_all", "Nr_sites_nonmissing", "fsum")]
            } else 
                f_block <- fsum_block_asc_cutoff[fsum_block_asc_cutoff$asc_freq_max==max_freq, c("Nr_sites_all", "Nr_sites_nonmissing", "fsum")]

            f_block_all <- f_block_all %>% rbind(data.frame(block_id=block_id_all, jackknife_weight=if (all_sites_block) f_block$Nr_sites_all else f_block$Nr_sites_nonmissing, f=f_block$fsum/f_block$Nr_sites_nonmissing, Asc_min=min_freq, Asc_max=max_freq))
        }
    }
} else
    f_block_all <- fsum_block_asc_cutoff %>% mutate(jackknife_weight=if (all_sites_block) Nr_sites_all else Nr_sites_nonmissing, f=fsum/Nr_sites_nonmissing) %>% select(block_id, jackknife_weight, f)


if (!is.na(out_file))
{
    write.table(f_block_all, out_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else
    print(f_block_all)
