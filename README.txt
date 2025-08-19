# Overview of Rare Allele Sharing (RAS) tool

This repository contains scripts used to generate the results of the preprint titled "High resolution analysis of population structure using rare variants" (url: https://www.biorxiv.org/content/10.1101/2025.07.18.665597v1). The RAS tool is implemented in R.

## Brief workflow

The tool accepts `.geno` file as input. Besides, we also define a self-designed allele count file called “freqSum”, consisting of two columns: alternative and non-missing allele count at each position, which is a concise summary for allele frequency.
Overall, the tool accepts the following input formats:
Type 1: genotype file in EIGENSTRAT (.geno)
Type 2: two-column freqSum format for alternative and non-missing allele count at each position
Type 3: one-column format for alternative allele count at each position, with fixed non-missing allele count which is indicated elsewhere
Type 4: one-column format for allele frequency at each position, with non-missing allele count 1 which is indicated elsewhere

We need a file to collect the frequency information for all individuals/populations involved in the analysis (including left/right populations, reference group and outgroup). This file contains two columns:
Column 1: The directory to frequency data (either in .geno file or freqSum format)
Column 2: (optional) If this is blank and Column 1 is a .geno file, the tool will automatically detect the .ind file in the same directory as Column 1 indicates. Otherwise, this should be the header of the frequency file indicated in Column 1 (separated by “,” If the frequency file has multiple columns). The header of the columns in frequency files could be “CHR”, “POS”, “ALT” or “REF”, corresponding to chromosome, position, alternative or reference, if those columns are not about frequency information. Otherwise, the header should have a bracket indicating the total allele count, followed by the column name. If the frequency information is part of freqSum format, we add “~num” or “~deno” on the column name, corresponding to “numerator” or “denominator”.

If there is a case when we need to define new groups by packing existing individuals/populations, we use a string separated by “:” to define this. On the left of the colon is the name of the new population; on the right of colon are the names of multiple existing populations separated by “,”.

In each of the left and right population, there could be either one group or multiple groups separated by “,”. Left population 1 and Right population 1 are mandatory, while Left population 2 and Right population 2 are optional (for RASD). The tool will conduct all combinations of statistics.

There are some optional filtering conditions:
We have "refpop_max_miss" to indicate the proportion of missingness in the reference group that can be tolerated.
We have “outgroup_max_hetero” to indicate the maximum heterozygosity in the outgroup that can be tolerated.

There are some options to adjust in calculating the statistical value if required:
For each left/right population, we have a “fillna” option to indicate how we treat the missing data.
We can set “outgroup_round_freq = True” to round the allele frequency of the outgroup to 0 or 1.

The description of code files and parameters is shown below.


## Overview of scripts

- `RAS_src.R`: private functions used for preprocessing the different types of input data.
- `get.blockf.R`: the computation tool to calculate block-wise statistical values, command line version.
- `get.blockf.function.R`: the computation tool to calculate block-wise statistical values, function version.
- `get_f_from_fblock.R`: functions to summerize block-wise statistical values.

## Arguments and parameters

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
