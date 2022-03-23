source("../../R/calculate_accuracy.R", chdir = TRUE)
library(testthat)

test_that("parameter which are not matrix or dataframe", {
  expect_error(calculate_accuracy(1,1), "true_labels_df argument is not a dataframe or matrix.")
  df <- data.frame()
  expect_error(calculate_accuracy(df,1), "predicted_labels_df argument is not a dataframe or matrix.")
})
