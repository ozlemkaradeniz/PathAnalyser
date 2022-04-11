test_that("classification method tests", {

data("ER_sig_df")
data("normalized_ER_data_se1")
data("ER_classes_df.default")

classes_df.default_actual<- classify_GSVA_percent(ER_sig_df, normalized_ER_data_se1)

expect_equal(classes_df.default_actual, ER_classes_df.default)

})
