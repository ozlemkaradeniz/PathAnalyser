test_that("classification method tests", {

data("ER_sig_df")
data("ER_data_se1")

normalized_se <- check_signature_vs_dataset(log_cpm_transformation(ER_data_se1), ER_sig_df)

classes_df.default<- classify_GSVA_percent(ER_sig_df, normalized_se)

table.default <- table(classes_df.default$class)
table.default <- as.numeric(table.default)

expect_identical(c(4,5,11), table.default)

})
