"""
A throwaway notebook to test merit.data_tools.calc_beta_from_z.

1. Load first 100 rows of raw eQTLgen cis-eqtl data and allele frequency data
2. Use merit to append beta and SE columns
3. Save to a file called eqtlgen_ciseqtl_head_100_merit.txt
"""

# %%
from merit.data_tools import calc_beta_from_z, get_unique_columns
from merit import constants as con
import pandas as pd
import numpy as np

# %%
# read raw data
eqtlgen_100 = pd.read_csv(
    "../tests/data/raw/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
    sep="\t",
    nrows=100)

allele_freq = pd.read_csv(
    "../tests/data/raw/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt",
    usecols=["AlleleB_all", "SNP"],
    sep="\t")

# %%
# rename allele frequency col in allele_freq
allele_freq = allele_freq.rename(columns={"AlleleB_all": "allele_freq"})

# set indices for merge
eqtlgen_100.set_index("SNP", inplace=True)
allele_freq.set_index("SNP", inplace=True)

# %%
# make my version of calc_beta_from_z


def calc_beta_from_z_aw(df, allele_freq="effect_allele_freq",
                        effect_size="effect_size",
                        no_of_samples="number_of_samples"):
    """
    Calculates standardized beta (correlation) coefficients from z scores.
    This uses the method from SMR

    Parameters
    ----------
    df : :obj:`pandas.DataFrame`
        The DataFrame containing the z_score effect size that we want to
         convert into a beta.
    allele_freq : :obj:`pandas.DataFrame` or str, optional, default: effect_allele_freq
        If this is a string then it is assumed that it is an allele frequency
        column of that name within df. If it is a DataFrame, it is assumed that
        it is a separate DataFrame of allele frequencies with the column
        `allele_freq` (which must not exist in df). Also, any allele frequency
        DataFrames must have an index that can be merged with df
    effect_size : str, optional, default: effect_size
        The name of the effect size column in df
    no_of_samples : str or int, optional, default: number_of_samples
        If it is a string, it is assumed that it is a number of samples column
        in df, with separate sample sized for each variant. If it is an integer
        it is assumed that it is an overall sample size for the cohort, that
        will be applied to each variant

    Returns
    -------
    df : :obj:`pandas.DataFrame`
        The DataFrame containing with the effect sizes represented as a beta
        and the standard errors calculated
    """
    allele_freq_col = "allele_freq"

    # FIXME
    ## no_of_samples_col = con.NUMBER_OF_SAMPLES.name
    no_of_samples_col = no_of_samples

    delete_no_of_samples = False
    delete_allele_freq = False
    unique_cols = get_unique_columns([no_of_samples_col,
                                      allele_freq_col], df.columns)

    if isinstance(no_of_samples, int):
        no_of_samples_col = unique_cols[no_of_samples_col]

        df[no_of_samples_col] = no_of_samples
        delete_no_of_samples = False

    if isinstance(allele_freq, pd.DataFrame):
        af = allele_freq[[allele_freq_col]]
        df.rename({allele_freq_col: unique_cols[allele_freq_col]}, axis=1,
                  inplace=True)
        allele_freq_col = unique_cols[allele_freq_col]
        df = df.merge(af, left_index=True, right_index=True)
        delete_allele_freq = False

    df[con.STANDARD_ERROR.name] = 1 / np.sqrt(
        2 * df[allele_freq_col] *
        (1 - df[allele_freq_col]) *
        (df[no_of_samples_col] + df[effect_size]**2)
    )

    df[effect_size] = df[effect_size] * df[con.STANDARD_ERROR.name]

    if delete_no_of_samples is True:
        df.drop(no_of_samples_col, axis=1, inplace=True)

    if delete_allele_freq is True:
        df.drop(allele_freq_col, axis=1, inplace=True)

    return df


# %%
# calculate beta/SE
merit_eqtlgen_100 = calc_beta_from_z_aw(
    df=eqtlgen_100,
    allele_freq=allele_freq,
    effect_size="Zscore",
    no_of_samples="NrSamples")

# %%
# save to .txt
merit_eqtlgen_100.to_csv("../tests/data/eqtlgen_ciseqtl_head_100_merit.txt",
                         sep="\t",
                         index=False)
