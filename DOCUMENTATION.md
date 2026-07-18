# Rare Allele Sharing (RAS) Tool

## What this tool does

This repository contains an R implementation of a block-based rare allele sharing
workflow. It reads allele-frequency data for multiple populations, applies optional
site filters, computes a per-site contrast statistic, and then summarizes that
statistic by jackknife block.

The tool supports both:

- a non-ascertained mode, and
- an ascertained mode that bins sites by derived allele frequency in a reference
  population.

It was written to support the analysis described in the associated preprint on
population structure using rare variants.

## Main scripts

- `RAS_src.R`
  - Core helper functions for reading input metadata, loading frequencies, defining
    custom populations, and computing block-level statistics.
- `get.blockf.R`
  - Command-line entry point.
  - Parses options, loads data, computes block statistics, and writes the result.
- `get.blockf.function.R`
  - Function wrapper around the same workflow.
  - Intended for programmatic use, though it appears not to be fully synchronized
    with the command-line script.
- `get_f_from_fblock.R`
  - Post-processing helper.
  - Reads block outputs and computes the weighted mean and jackknife standard error.

## Inputs

### 1. Frequency file manifest

The primary input is a tab-delimited manifest listing frequency sources.

Each row contains:

- column 1: path to a frequency file
- column 2: optional header information

If column 2 is omitted and the file is `.geno`, the script looks for an `.ind`
file with the same basename. Otherwise it looks for a `.header` file, or interprets
the string in column 2 as a comma-separated header specification.

### 2. Supported frequency formats

The tool accepts several formats:

- `.geno` files in EIGENSTRAT format
- two-column `freqSum` files, where the two columns represent numerator and
  denominator
- one-column allele-count files
- one-column allele-frequency files

The column headers can encode the total allele count in parentheses.

### 3. Population definition file

Optional file for defining synthetic populations from existing ones.

Format:

- `NewPop:PopA,PopB,PopC`

This creates a new population named `NewPop` by combining the listed existing
populations.

### 4. Position file

Optional two-column file containing `CHR` and `POS`.

If present, it restricts the analysis to the listed sites.

### 5. Ascertainment frequency table

Optional two-column file with:

- `min_freq`
- `max_freq`

Used only in ascertained mode to define derived-frequency bins.

## Command-line interface

`get.blockf.R` accepts the following key options:

- `--frequencyFilesList`
  - required
  - path to the frequency file manifest
- `--positionFile`
  - optional
  - site subset file
- `--popDefFile`
  - optional
  - custom population definitions
- `--popLeft1`, `--popRight1`
  - required
  - primary left and right population groups
- `--popLeft2`, `--popRight2`
  - optional
  - secondary left and right population groups for RASD contrasts
- `--popLeftFillna`, `--popRightFillna`
  - optional
  - missing-data handling for left and right populations
- `--jackknife`
  - optional
  - block definition, either `"CHR"` or an integer window length
- `--ascertained`
  - optional flag
  - enables ascertained analysis
- `--ascRefpop`
  - optional
  - reference population for derived-frequency ascertainment
- `--ascFreqTable`
  - optional
  - ascertainment bin table
- `--ascOutgroup`
  - required
  - outgroup population
- `--refpopMaxMiss`
  - optional
  - maximum tolerated missingness in the ascertainment reference population
- `--refpopRemoveHomo`
  - optional flag
  - removes homozygous sites in the ascertainment reference population
- `--outgroupMaxHetero`
  - optional
  - maximum tolerated outgroup heterozygosity
- `--outgroupRoundFreq`
  - optional flag
  - rounds outgroup frequency to 0 or 1 before ascertained calculations
- `-o`, `--outFile`
  - optional
  - output path; if omitted, results are printed

## Workflow

### Step 1: Read frequency metadata

`get_frequency_files_info()` reads the manifest and expands it into per-column
metadata.

It determines:

- which files are genotype files versus frequency files,
- which columns belong to positions versus populations,
- the denominator for each population column,
- which populations are represented by one column or two columns.

### Step 2: Resolve population names

The user-supplied left/right population strings are parsed as either:

- comma-separated lists, or
- files containing one population per line.

Custom populations from the population-definition file are added to the valid
population universe.

### Step 3: Build population bookkeeping

The script constructs `full_info$pops`, a table containing:

- `Nr_columns`
  - number of input columns used for the population
- `total`
  - total allele count
- `keep_frac`
  - whether numerator/denominator data should be retained
- `lth`
  - cache marker used by the lazy loading logic

### Step 4: Load the outgroup and build the site mask

The outgroup is loaded first.

Sites are retained if they satisfy:

- non-missing outgroup frequency
- optional outgroup heterozygosity filter
- if ascertained:
  - optional reference-population missingness filter
  - optional reference-population homozygote exclusion

If a position file is supplied, the retained sites must also appear in that file.

### Step 5: Load all needed populations

`get_freq()` lazily loads the remaining populations only for the surviving sites.

It can also materialize custom populations by summing the numerators and
denominators of their component populations.

### Step 6: Compute per-site values

For a site `i`, the core statistic is:

\[
f_i = (O_i - L_i)(O_i - R_i)
\]

where:

- `O_i` is the outgroup frequency,
- `L_i` is the left frequency or contrast,
- `R_i` is the right frequency or contrast.

In non-ascertained mode, the tool may define:

\[
L_i = L2_i - L1_i
\quad\text{and}\quad
R_i = R2_i - R1_i
\]

or, if only one population is provided on a side:

\[
L_i = O_i - L1_i
\quad\text{and}\quad
R_i = O_i - R1_i
\]

### Step 7: Aggregate by jackknife block

Each site is assigned to a block using:

- chromosome-only blocking if `jackknife == "CHR"`
- otherwise a windowed block index based on chromosome and position

For each block, the tool computes:

- `Nr_sites_all`
- `Nr_sites_nonmissing`
- `fsum`
- `f = fsum / Nr_sites_nonmissing`

### Step 8: Ascertained branch

If `--ascertained` is enabled, the code:

- derives a frequency class from the reference population,
- bins sites using `ascFreqTable`,
- computes block values within each bin,
- optionally builds RASD contrasts from the resulting RAS tables.

## Output

The main output is a tab-delimited table containing block-wise results.

Typical columns include:

- `min_freq`
- `max_freq`
- `asc_outgroup`
- `pop_left1`
- `pop_left2`
- `pop_right1`
- `pop_right2`
- `block_id`
- `Nr_sites_all`
- `Nr_sites_nonmissing`
- `f`

The summarization helper `get_f_from_fblock()` can then combine block rows into a
single estimate with a jackknife standard error.

## Summary statistic

Given block values `f_b` and weights `n_b`, the final estimate is:

\[
\bar f = \frac{\sum_b n_b f_b}{\sum_b n_b}
\]

with jackknife standard error computed from the block deviations.

## External dependencies

The code uses:

- `dplyr`
- `purrr`
- `tidyr`
- `data.table`
- `readr`
- `optparse`
- base R packages, including `tools`

## Notes and caveats

- `get.blockf.R` sources `RAS_src.R` from a hard-coded absolute path, so the
  script is not portable as written.
- `get.blockf.function.R` appears to be an older or incomplete variant of the
  command-line workflow. It references helper functions that are not present in
  the files in this repository.
- `allSitesBlock` is parsed in the CLI and function signatures, but it is not
  clearly used in the main command-line workflow.

## Suggested reading order

If you are trying to understand the code from scratch, inspect files in this order:

1. `README.txt`
2. `get.blockf.R`
3. `RAS_src.R`
4. `get_f_from_fblock.R`

