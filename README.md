# RAStools

RAStools computes block-based rare allele sharing statistics from population
allele-frequency and EIGENSTRAT-style genotype inputs.

## Basic use

```r
library(RAStools)

f_blocks <- get_blockf(
  base_dir = "C:/lei_huang/RAS_tools",
  frequency_files_list = "test_example/RAS_empirical/frequency_files_list.txt",
  pop_def_file = "test_example/RAS_empirical/pop_def.txt",
  pop_left1 = "HGDP01381,HGDP01382",
  pop_right1 = "CEU,FIN,GBR",
  pop_right2 = "IBS,TSI",
  jackknife = "CHR",
  is_ascertained = TRUE,
  asc_refpop = "EUR5",
  asc_freq_table = "test_example/RAS_empirical/asc_freq_table.txt",
  refpop_max_miss = 0.75,
  refpop_remove_homo = TRUE,
  asc_outgroup = "AFR_all",
  outgroup_max_hetero = 0
)

summary <- get_f_from_fblock(f_blocks, return_result = TRUE)
```

The repository includes a large local example dataset under `data/` and
`test_example/`. Those files are intentionally excluded from package builds so
the installable package stays small; pass your own manifest and data paths when
using an installed copy.
