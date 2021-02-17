#' Append columns for beta and standard error values to a dataframe of eQTLGen
#' cis-eQTL data
#'
#' A wrapper around \code{\link{mutate_beta_se_from_zscore}}. Works specifically
#' for the eQTLGen
#' \href{https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_cis}{cis-eQTL}
#' and
#' \href{https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_allele_frequency}{cis_eQTL
#' mean allele frequency} datasets.
#'
#' @param eqtlgen_cis_eqtl_df Dataframe of eQTLGen cis-eQTL data
#' @param eqtlgen_cis_eqtl_maf_df eQTLGen cis-eQTL mean allele frequency data
#'
#' @return A dataframe with appended \code{beta} and \code{se} columns
#' @export
#'
mutate_beta_se_from_zscore_eqtlgen_cis_eqtl <- function(eqtlgen_cis_eqtl_df,
                                                        eqtlgen_cis_eqtl_maf_df) {
  eqtlgen_cis_eqtl_df %>%
    dplyr::left_join(eqtlgen_cis_eqtl_maf_df, by = "SNP") %>%
    mutate_beta_se_from_zscore(zscore = "Zscore",
                               maf = "AlleleB_all",
                               sample_n = "NrSamples",
                               beta_se = "beta")
}


#' Append columns for beta and standard error values to a dataframe
#'
#' Appends columns with beta and standard error values to a dataframe containing
#' columns for Z-scores, allele frequencies and sample sizes.
#'
#' @inheritParams beta_se_from_zscore
#'
#' @return A dataframe with appended \code{beta} and \code{se} columns
#' @export
#'
#' @seealso \code{\link{beta_se_from_zscore}}
mutate_beta_se_from_zscore <- function(df,
                                       zscore,
                                       maf,
                                       sample_n,
                                       beta_se) {
  # mutate beta column
  df$beta <- beta_se_from_zscore(
    df = df,
    zscore = zscore,
    maf = maf,
    sample_n = sample_n,
    beta_se = "beta"
  )

  # mutate se column
  df$se <- beta_se_from_zscore(
    df = df,
    zscore = zscore,
    maf = maf,
    sample_n = sample_n,
    beta_se = "se"
  )

  return(df)
}

#' Get beta or standard error from z-score and allele frequency
#'
#' Take a dataframe containing Z-scores, allele frequencies and sample sizes and
#' returns a vector of either beta or standard error values (specified by the
#' `beta_se` argument)
#'
#' @param df A dataframe containing the required input values
#' @param zscore Colname from \code{df} containing zscores
#' @param maf Colname from \code{df} containing allele frequencies
#' @param sample_n Colname from \code{df} containing sample size
#' @param beta_se Character. Must be either \code{"beta"} or \code{"se"}
#'
#' @return A numerical vector of either beta or standard error values
#' @export

beta_se_from_zscore <- function(df,
                                zscore,
                                maf,
                                sample_n,
                                beta_se) {
    # return beta
  if (beta_se == 'beta') {
    df[[zscore]] /
      sqrt(2 * df[[maf]] * (1 - df[[maf]]) * (df[[sample_n]] + df[[zscore]] ^ 2))

    # return se
  } else if (beta_se == "se") {
    1 /
      sqrt(2 * df[[maf]] * (1 - df[[maf]]) * (df[[sample_n]] + df[[zscore]] ^ 2))

    # error if invalid beta_se value
  } else {
    stop("Argument 'beta_se' must be either 'beta' or 'se'")
  }

}
