test_that("classification method tests", {

data("ER_sig_df")
data("ER_data_se1")
load("normalized_ER_data_se1.rda")
load("ER_classes_df.default.rda")

classes_df.default_actual<- classify_GSVA_percent(ER_sig_df, normalized_ER_data_se1)

expect_equal(classes_df.default_actual, ER_classes_df.default)

})

