#' Classifies samples according to pathway activity using GSVA
#' Creates confusion matrix and calculate accuracy of classification method
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
        
        # create confusion matrix
        confusion_matrix_HER2<-matrix(0,nrow = 3,ncol = 3,
                                      dimnames = list(c("Actual_Pos","Actual_Neg","Actual_Uncertain"),
                                                      c("Prediction_Pos","Prediction_Neg", "Prediction_Uncertain")))
        df<-merge(x=true_labels_df, y=predicted_labels_df, by.x = "CaseID", by.y="sample")
        
        true_positive <- nrow(df[(df$class == "Positive" & df$HER2 == "Positive"),]) # predicted positive, actual positive
        false_positive <- nrow(df[(df$class == "Positive" & df$HER2 == "Negative"),]) # predicted positive, actual negative
        confusion_matrix_HER2[1,1] <- true_positive
        confusion_matrix_HER2[2,1] <- false_positive
        confusion_matrix_HER2[3,1] <- nrow(df[df$class == 'Positive',])-true_positive-false_positive
        

        true_negative <- nrow(df[(df$class == "Negative" & df$HER2 == "Negative"),]) # predicted negative, actual negative
        false_negative <- nrow(df[(df$class == "Negative" & df$HER2 == "Positive"),]) # predicted negative, actual positive
        confusion_matrix_HER2[2,2] <- true_negative
        confusion_matrix_HER2[1,2] <- false_negative
        confusion_matrix_HER2[3,2] <- nrow(df[df$class == 'Negative',])-true_negative-false_negative
        
        true_uncertain <- nrow(df[(df$class == "Uncertain" & df$HER2 != "Uncertain"),]) # predicted uncertain, actual uncertain
        confusion_matrix_HER2[3,3] <- true_uncertain
        confusion_matrix_HER2[1,3] <- nrow(df[(df$class == 'Uncertain' & df$HER2 == "Positive"),])  # predicted uncertain, actual positive
        false_uncertain <- confusion_matrix_HER2[1,3]
        confusion_matrix_HER2[2,3] <- nrow(df[(df$class == 'Uncertain' & df$HER2 == "Negative"),])  # predicted uncertain, actual negative
        false_uncertain <- false_uncertain + confusion_matrix_HER2[2,3]
        accuracy <- (true_positive + true_negative + true_uncertain) / nrow(df) * 100
        
        print("Confusion Matrix ")
        print("--------------------------------------------------------------")
        print(confusion_matrix_HER2)

        print("Accuracy of classification summary ")
        print("--------------------------------------------------------------")
        print(paste0("accuracy: ", format(round(accuracy, 2), nsmall = 2)))
        print(paste0("total number samples: " , nrow(df)))
        print(paste0("True Positive(TP): " , true_positive))
        print(paste0("True Negative(TN): " , true_negative))
        print(paste0("True Uncertain: " , true_uncertain))
        print(paste0("False Negative(FN): " , false_negative))
        print(paste0("False Positive(FP): " , false_positive))
        print(paste0("False Uncertain: " , false_uncertain))
        print("--------------------------------------------------------------")
        print(paste0("True Positive Rate(TPR): ", 
                     format(round(true_positive / (true_positive + false_negative) * 100 , 2) , nsmall =2)))
        print(paste0("True Negative Rate(TNR): ", 
             format(round(true_negative / (true_negative + false_positive) * 100 , 2) , nsmall =2)))
        print(paste0("False Positive Rate(FPR): ", 
             format(round(false_positive / (false_positive + true_negative) * 100 ,  2) , nsmall =2)))
        print(paste0("False Negative Rate(FNR): ", 
             format(round(false_negative / (false_negative + true_positive) * 100 ,  2) , nsmall =2)))
        print("--------------------------------------------------------------")
        print(summary(confusion_matrix_HER2))
        
        return(confusion_matrix_HER2)
} 

































