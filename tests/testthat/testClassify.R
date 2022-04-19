test_that("classification method tests", {

data("ER_sig_df")
data("ER_TCGA_RNAseq")
load("normalized_ER_TCGA_RNAseq.rda")
load("ER_classes_df.default_actual.rda")

ER_classes_df.default<- classify_GSVA_percent(ER_sig_df, normalized_ER_TCGA_RNAseq)

expect_equal(ER_classes_df.default_actual, ER_classes_df.default)

})



