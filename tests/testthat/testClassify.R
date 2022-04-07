test_that("classification method tests", {

ER_expression_set <- read_expression_data("../../raw_data/toy_data.txt")
er_sign_df <- read_signature_data("../../raw_data/ESR1_UP.v1._UP.grp", "../../raw_data/ESR1_DN.v1_DN.grp")

normalized_se <- check_signature_vs_dataset(log_cpm_transformation(ER_expression_set), er_sign_df)

classes_df.default<- classify_GSVA_percent(er_sign_df, normalized_se)

table.default <- table(classes_df.default$class)
table.default <- as.numeric(table.default)

expect_identical(c(4,5,11), table.default)

classes_df.percent <- classify_GSVA_percent(er_sign_df, normalized_se,
                                            percent_thresh=30)

table.percent <- table(classes_df.percent$class)
table.percent <- as.numeric(table.percent)

expect_identical(c(6,6,8), table.percent)

classes_df.absolute <- classify_GSVA_abs(er_sign_df, normalized_se,
                                up_thresh.low=-0.2, up_thresh.high=0.2,
                                dn_thresh.low=-0.2, dn_thresh.high=0.2)

table.absolute <- table(classes_df.absolute$class)
table.absolute <- as.numeric(table.absolute)

expect_identical(c(9,10,1), table.absolute)

})
