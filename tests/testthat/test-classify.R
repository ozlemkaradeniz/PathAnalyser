library(testthat)
source_file("R/classification.R")
data("ER_sig")
data("ER_data_se1")

test_that("signature is not a dataframe", {
  expect_error(classify(1, ER_data_se),
               "Signature argument is not a dataframe.")
  expect_error(classify(TRUE, ER_data_se),
               "Signature argument is not a dataframe.")
  mat <- matrix(nrow = 1, ncol = 2)
  expect_error(classify(mat, ER_data_se),
               "Signature argument is not a dataframe.")
})

test_that("expression matrix is not a dataframe or matrix", {
  expect_error(classify(ER_sig_df, 1),
               "Expression dataset is not a matrix or dataframe object.")
  expect_error(classify(ER_sig_df, TRUE),
              "Expression dataset is not a matrix or dataframe object.")
})

test_that("Incorrect thresholds generate appropriate errors", {
  expect_error(classify(ER_sig_df, ER_data_se, up_thresh=c(0, 2),
                        dn_thresh=c(-0.4, 0.4)),
               "Up-regulated gene set thresholds are not between -1 and 1.")
  expect_error(classify(ER_sig_df, ER_data_se, up_thresh=c(-0.4, 0.4),
                        dn_thresh=c(-1.5, 0.5)),
               "Up-regulated gene set thresholds are not between -1 and 1.")
  expect_error(classify(ER_sig_df, ER_data_se, up_thresh=c(-0.4, 0.7)),
               "up_thresh argument provided but not dn_thresh argument.")
  expect_error(classify(ER_sig_df, ER_data_se, dn_thresh=c(-0.4, 0.7)),
               "dn_thresh argument provided but not up_thresh argument.")
})
