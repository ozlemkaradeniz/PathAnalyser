#' Calculation of accuracy by  confusion matrix
#' @description Creates confusion matrix and calculate accuracy of classification method
#' @author Ozlem Karadeniz \email{ozlem.karadeniz.283@cranfield.ac.uk}
#' @param true_labels_df actual labels dataframe or matrix
#' @param predicted_labels_df predicted labels dataframe or matrix generated in classification method
#'
#' @return confusion_matrix_HER2
#' @export
#'
#' @examples
#' calculate_accuracy(HER2_labelDF,Her2_pathway)

calculate_accuracy <-function(true_labels_df, predicted_labels_df){

        # check true_classes_df arg is dataframe or matrix
        if (is.matrix(true_labels_df)){
                true_labels_df <- as.data.frame(true_labels_df)
        }
        else if (!is.data.frame(true_labels_df)){
                stop("true_labels_df argument is not a dataframe or matrix.")
        }


        # check predicted_labels_df arg is dataframe or matrix
        if (is.matrix(predicted_labels_df)){
                predicted_labels_df <- as.data.frame(predicted_labels_df)
        }

        else if (!is.data.frame(predicted_labels_df)){
                stop("predicted_labels_df argument is not a dataframe or matrix.")
        }

        # HER2 status Equivocal, [Not Evaluated] and Indeterminate are all set to uncertain

        true_labels_df$HER2 <- sapply(true_labels_df$HER2,
                 function(x){gsub(pattern = "(Equivocal|\\[Not Evaluated\\]|Indeterminate)",
                                  replacement = "Uncertain", x)})

        for (i in 1:nrow(predicted_labels_df)) {
                predicted_labels_df[i,1] <- gsub("\\.", "-", predicted_labels_df[i,1])
        }


        df<-merge(x=true_labels_df, y=predicted_labels_df, by.x = "CaseID", by.y="sample")

        TP <- nrow(df[(df$class == "Active" & df$HER2 == "Positive"),]) # predicted positive, actual positive
        FP <- nrow(df[(df$class == "Active" & df$HER2 == "Negative"),]) # predicted positive, actual negative

        TN <- nrow(df[(df$class == "Inactive" & df$HER2 == "Negative"),]) # predicted negative, actual negative
        FN <- nrow(df[(df$class == "Inactive" & df$HER2 == "Positive"),]) # predicted negative, actual positive

        prd_uncert_act_positive = nrow(df[(df$class == "Uncertain" & df$HER2 == "Positive"),]) # predicted uncertain, actual positive
        prd_uncert_act_negative = nrow(df[(df$class == "Uncertain" & df$HER2 == "Negative"),]) # predicted uncertain, actual negative
        prd_positive_act_uncertain = nrow(df[(df$class == "Active" & df$HER2 == "Uncertain"),]) # predicted positive, actual uncertain
        prd_negative_act_uncertain = nrow(df[(df$class == "Inactive" & df$HER2 == "Uncertain"),]) # predicted negative, actual uncertain
        prd_uncert_act_uncertain = nrow(df[(df$class == "Uncertain" & df$HER2 == "Uncertain"),]) # predicted uncertain, actual uncertain

        # create confusion matrix
        matrix_data <- c(TP, FN, prd_uncert_act_positive, FP, TN, prd_uncert_act_negative, prd_positive_act_uncertain, prd_negative_act_uncertain, prd_uncert_act_uncertain)
        confusion_matrix_HER2<-matrix(matrix_data,nrow = 3,ncol = 3,
                                      dimnames = list(c("Prediction Positive ","Prediction Negative", "Prediction Uncertain"),
                                                      c("Actual Positive","Actual Negative" , "Actual Uncertain")))

        classified_samples_proportion <- (TP + FN + FP + TN) / sum(confusion_matrix_HER2) * 100
        accuracy_amongst_classified_samples <- (TP + TN) / (TP + FN + FP + TN) * 100

        print("Confusion Matrix ")
        print("--------------------------------------------------------------")
        print(confusion_matrix_HER2)

        print("Statistics in Confusion Matrix ")
        print("--------------------------------------------------------------")
        print(paste0("Proportion of classified samples: ", format(round(classified_samples_proportion, 2), nsmall = 2)))
        print(paste0("Accuracy amongst classified samples: " , format(round(accuracy_amongst_classified_samples, 2), nsmall = 2)))
        print(paste0("True Positive(TP): " , TP))
        print(paste0("True Negative(TN): " , TN))
        print(paste0("False Negative(FN): " , FN))
        print(paste0("False Positive(FP): " , FP))
        print("--------------------------------------------------------------")
        print(paste0("True Positive Rate(TPR)(sensitivity)(Recall): ",
                     format(round(TP / (TP + FN) * 100 , 2) , nsmall =2)))
        print(paste0("True Negative Rate(TNR)(specificity): ",
             format(round(TN / (TN + FP) * 100 , 2) , nsmall =2)))
        print(paste0("Precision (Positive predictive value): ",
                     format(round(TP / (TP + FP) * 100 , 2) , nsmall =2)))
        print(paste0("False Positive Rate(FPR): ",
             format(round(FP / (FP + TN) * 100 ,  2) , nsmall =2)))
        print(paste0("False Negative Rate(FNR): ",
             format(round(FN / (FN + TP) * 100 ,  2) , nsmall =2)))
        print("--------------------------------------------------------------")
        print(summary(confusion_matrix_HER2))

        return(confusion_matrix_HER2)
}

































