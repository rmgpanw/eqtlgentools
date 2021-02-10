library(eqtlgentools)

test_df <- data.frame(
  zscore = c(200.7534, 200.6568, 200.2654),
  maf = c(0.2554750, 0.2554561, 0.2543986),
  sample_n = c(30596, 30596, 30598)
)


test_that("Correct beta values are generated", {
  expect_equal(beta_se_from_zscore(test_df, "zscore", "maf", "sample_n", "beta"),
               c(1.222411, 1.222187, 1.222806),
               tolerance = 1e-06)
})

test_that("Correct standard error values are generated", {
  expect_equal(beta_se_from_zscore(test_df, "zscore", "maf", "sample_n", "se"),
               c(0.006089118, 0.006090932, 0.006105927),
               tolerance = 1e-07)
})
