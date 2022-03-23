library(testthat)
source_file("R/classification.R")
data("ER_sig_df")
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
  expect_error(classify(ER_sig_df, ER_data_se, up_thresh=c(-0.4, 1.4),
                        dn_thresh=c(-0.4, 0.4)),
               "Gene set thresholds for classification must be between -1 and 1.")
  expect_error(classify(ER_sig_df, ER_data_se, up_thresh=c(-1.4, 0.4),
                        dn_thresh=c(-0.4, 0.4)),
               "Gene set thresholds for classification must be between -1 and 1.")
  expect_error(classify(ER_sig_df, ER_data_se, up_thresh=c(-0.4, 0.4),
                        dn_thresh=c(-0.4, 1.4)),
               "Gene set thresholds for classification must be between -1 and 1.")
  expect_error(classify(ER_sig_df, ER_data_se, up_thresh=c(-0.4, 0.4),
                        dn_thresh=c(-1.4, 0.4)),
               "Gene set thresholds for classification must be between -1 and 1.")
})
