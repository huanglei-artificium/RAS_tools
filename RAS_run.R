library(RAStools)


# RAS_empirical; ascertained

position_file <- NA

frequency_files_list <- "./test_example/RAS_empirical/frequency_files_list.txt"
# allowed non-frequency columns: CHR; POS; ALT; REF

pop_def_file <- "./test_example/RAS_empirical/pop_def.txt"

pop_left1_str <- "HGDP01381,HGDP01382"

pop_left2_str <- NA

pop_left_fillna <- NA

pop_right1_str <- "CEU,FIN,GBR"

pop_right2_str <- "IBS,TSI"

pop_right_fillna <- NA

jackknife <- "CHR"

is_ascertained <- TRUE

asc_refpop <- "EUR5"

asc_freq_table <- "./test_example/RAS_empirical/asc_freq_table.txt"

all_sites_block <- FALSE # not applied anymore

refpop_max_miss <- 0.75

refpop_remove_homo <- TRUE

asc_outgroup <- "AFR_all"

outgroup_max_hetero <- 0

outgroup_round_freq <- FALSE
# use the rounded frequency (0/1) when outgroup is involved in F-stats


f_block_table_all <- get_blockf(
    position_file = position_file,
    frequency_files_list = frequency_files_list,
    pop_def_file = pop_def_file,
    pop_left1 = pop_left1_str,
    pop_left2 = pop_left2_str,
    pop_left_fillna = pop_left_fillna,
    pop_right1 = pop_right1_str,
    pop_right2 = pop_right2_str,
    pop_right_fillna = pop_right_fillna,
    jackknife = jackknife,
    is_ascertained = is_ascertained,
    asc_refpop = asc_refpop,
    asc_freq_table = asc_freq_table,
    all_sites_block = all_sites_block,
    refpop_max_miss = refpop_max_miss,
    refpop_remove_homo = refpop_remove_homo,
    asc_outgroup = asc_outgroup,
    outgroup_max_hetero = outgroup_max_hetero,
    outgroup_round_freq = outgroup_round_freq
)


f_block_table_all
