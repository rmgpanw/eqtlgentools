#' Append columns for beta and standard error values to the eQTLGen cis-eQTL
#' dataset.
#'
#' Calculates beta and standard error values for the eQTLGen eQTLGen
#' \href{https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_cis}{cis-eQTL}
#' dataset from Z-score, sample size and allele frequency (these values are in a
#' \href{https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_allele_frequency}{
#' separate file}) values.
#'
#' @param eqtlgen_cis_eqtl_filepath Character. Path to
#'   \href{https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_cis}{cis-eQTL
#'    file}
#' @param eqtlgen_cis_eqtl_maf_filepath Character. Path to
#'   \href{https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_allele_frequency}{cis_eQTL
#'    mean allele frequency file}
#' @inheritParams process_file_chunked
#' @inheritParams readr::read_delim_chunked
#'
#' @export
#' @seealso \code{\link{process_file_chunked}},
#'   \code{\link{mutate_beta_se_from_zscore_eqtlgen_cis_eqtl}}
append_beta_se_to_eqtlgen_cis_eqtl <- function(eqtlgen_cis_eqtl_filepath,
                                               eqtlgen_cis_eqtl_maf_filepath,
                                               outfile_path,
                                               chunk_size = 10000) {
  # read required cols from allele frequency file
  message("Reading allele frequency file")
  allele_freq <- data.table::fread(file = eqtlgen_cis_eqtl_maf_filepath,
                       select = c("AlleleB_all", "SNP"))

  # process file in chunks
  process_file_chunked(
    infile_path = eqtlgen_cis_eqtl_filepath,
    outfile_path = outfile_path,
    data_processing_function = mutate_beta_se_from_zscore_eqtlgen_cis_eqtl,
    chunk_size = chunk_size,
    eqtlgen_cis_eqtl_maf_df = allele_freq
  )
}

#' Apply a processing function to a file and output the result to a new file in
#' chunks
#'
#' A wrapper around \code{\link[readr]{read_delim_chunked}} and
#' \code{\link[readr]{write_delim}}. Allows a large file of tabular data to be
#' processed in chunks that fit into memory.
#'
#' @param infile_path Character. Path to infile.
#' @param outfile_path Character. Path to outfile. An error is raised if a file
#'   already exists at this location.
#' @param data_processing_function Function. A function that accepts a dataframe
#'   as its first argument and returns a (processed) dataframe also
#' @inheritParams readr::write_delim
#' @inheritParams readr::read_delim_chunked
#' @param ... Additional parameters are passed on to
#'   \code{data_processing_function}
#' @param verbose Time taken message printed after each chunk has processed. Default is \code{TRUE}.
#'
#' @export
process_file_chunked <- function(infile_path,
                                 outfile_path,
                                 data_processing_function,
                                 chunk_size = 10000,
                                 delim = "\t",
                                 col_types = NULL,
                                 # indexes = NULL,
                                 verbose = TRUE,
                                 append = FALSE,
                                 quote = "\"",
                                 escape_backslash = FALSE,
                                 escape_double = TRUE,
                                 col_names = TRUE,
                                 locale = readr::default_locale(),
                                 na = c("", "NA"),
                                 quoted_na = TRUE,
                                 comment = "",
                                 trim_ws = FALSE,
                                 skip = 0,
                                 guess_max = min(1000, chunk_size),
                                 progress = readr::show_progress(),
                                 skip_empty_rows = TRUE,
                                 ...) {

  start_time <- proc.time()

  # Error message if file already exists
  if (
    file.exists(outfile_path)) {
    stop("Error! ", outfile_path, " already exists.")
  }

  # Chunk processing function - applies data_processing_function then writes to outfile
  f <- function(x,
                pos) {
    time_taken <- proc.time() - start_time

    if (verbose == TRUE) {
      message(
        "Writing from line ",
        pos,
        ". Time taken: ",
        (time_taken[3] %/% 60),
        " minutes, ",
        (round(time_taken[3] %% 60)),
        " seconds"
      )
    }

    x <- data_processing_function(x, ...)

    # write to outfile
    readr::write_delim(
      x = x,
      file = outfile_path,
      delim = delim,
      na = "NA",
      append = TRUE,
      col_names = ifelse(pos == 1, TRUE, FALSE),
      quote_escape = "double",
      eol = "\n",
    )
  }

  # Read file and write to table in chunks
  message("Writing chunk to ", outfile_path)
  readr::read_delim_chunked(
    file = infile_path,
    callback = DataFrameCallback$new(f),
    chunk_size = chunk_size,
    delim = delim,
    col_types = col_types,
    quote = quote,
    escape_backslash = escape_backslash,
    escape_double = escape_double,
    col_names = col_names,
    locale = locale,
    na = na,
    quoted_na = quoted_na,
    comment = comment,
    trim_ws = trim_ws,
    skip = skip,
    guess_max = guess_max,
    progress = progress,
    skip_empty_rows = skip_empty_rows
  )

  # Completion message
  time_taken <- proc.time() - start_time
  message("Complete. Time taken: ",
          (time_taken[3] %/% 60),
          " minutes, ",
          (round(time_taken[3] %% 60)),
          " seconds")
}
