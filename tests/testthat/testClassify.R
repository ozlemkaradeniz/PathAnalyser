test_that("classification method tests", {

data("ER_SET_sig")
data("ER_TCGA_RNAseq")
load("normalized_ER_TCGA_RNAseq.rda")
load("ER_classes_df.default_actual.rda")

ER_classes_df.default<- classify_gsva_percent(normalized_ER_TCGA_RNAseq, ER_SET_sig)

expect_equal(ER_classes_df.default_actual, ER_classes_df.default)

})



