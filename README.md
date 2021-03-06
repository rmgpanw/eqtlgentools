
<!-- README.md is generated from README.Rmd. Please edit that file -->

# eqtlgentools

<!-- badges: start -->
<!-- badges: end -->

The goal of eqtlgentools is to generate beta and standard error values
for the
[cis-eQTL](https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_cis)
dataset from the [eQTLGen Consortium](https://www.eqtlgen.org/). This
requires [allele
frequencies](https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_allele_frequency)
from a separate file. Files are publicly available from
[here](https://www.eqtlgen.org/cis-eqtls.html)

The calculations for converting Z-scores to betas and standard errors
are documented in this
[here](https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/SMR_formatted/README_cis_full_SMR).

## Installation

Install from [github](https://github.com/rmgpanw/eqtlgentools) with:

``` r
devtools::install_github("rmgpanw/eqtlgentools")
```

## Usage

Run the following code (replacing the file paths with your own) in R to
generate a new file with beta and standard error columns appended to the
main cis-eQTL file. Adjust the `chunk_size` parameter if needed.

``` r
append_beta_se_to_eqtlgen_cis_eqtl(
  eqtlgen_cis_eqtl_filepath = "PATH/TO/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
  eqtlgen_cis_eqtl_maf_filepath = "PATH/TO/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt",
  outfile_path = "PATH/TO/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded-beta-se.txt",
  chunk_size = 1000000
)
```

For a full list of available functions, please see
[`Reference`](https://rmgpanw.github.io/eqtlgentools/reference/index.html)
section on [package website](https://rmgpanw.github.io/eqtlgentools/).
