position_file: (optional) The directory of the file consisting of two columns on chromosome and position. If indicated, this will be priority to the position information in “frequency_files_list”.
frequency_files_list: (mandatory) The directory of the file indicating all files on allele frequency.
pop_def_file: (optional) The directory of the file indicating the newly defined populations.
pop_left1: (mandatory) The string separated by “,” indicating all groups served as Left Population 1.
pop_left1_fillna: (optional) String. “HomoRef” / “Outgroup” / NA. Default NA.
pop_left2: (optional) The string separated by “,” indicating all groups served as Left Population 2.
pop_left2_fillna: (optional) String. “HomoRef” / “Outgroup” / NA. Default NA.
pop_right1: (mandatory) The string separated by “,” indicating all groups served as Right Population 1.
pop_right1_fillna: (optional) String. “HomoRef” / “Outgroup” / NA. Default NA.
pop_right2: (optional) The string separated by “,” indicating all groups served as Right Population 2.
pop_right2_fillna: (optional) String. “HomoRef” / “Outgroup” / NA. Default NA.
jackknife: (optional) “CHR” / an integer about the length. Default “CHR”.
asc_refpop: (optional) The name of the population that we would like to ascertain alleles with specific frequency cutoffs.
asc_freq_table: (optional) The directory of the file, with each row indicating the minimum and maximum derived allele frequency cutoffs.
all_sites_block: (optional) Boolean. Whether we use the number of all alleles rather than the number of non-missing alleles in jackknifing. Default False.
refpop_max_miss: (optional) Double between 0 and 1. The maximum allowed missing data for “asc_refpop” at each position. Default 0.
refpop_remove_homo: (optional) Whether we remove homozygote sites in Reference Population. Default False.
asc_outgroup: (optional) The name of the population that we would like to serve as the outgroup.
outgroup_max_hetero: (optional) Double between 0 and 1. The maximum allowed heterozygosity for “asc_outgroup” at each position. Default 0.
outgroup_round_freq: (optional) Boolean. Whether we round the frequency of the outgroup to 0 or 1. Default False.
