# Clear workspace
rm(list=ls()) 

install.packages("BiocManager")
BiocManager::install("GSVA")

# Installing
install.packages("readr")
# Loading
library("readr")

library(GSVA)

setwd("~/GroupProject1")


X <- read.table("TCGA_unannotated.txt", header=TRUE,sep = "\t")
X <- data.matrix(X)

up <- read_lines("ERBB2_UP.V1_UP.grp", skip = 3, n_max = -1L)
dn <- read_lines("ERBB2_UP.V1_DN.grp", skip = 3, n_max = -1L)

gs_list <- list("up"=up, "dn" = dn)

gsva.es <- gsva(X, gs_list, verbose=FALSE)
write.csv(gsva.es, "/Users/ozlemkaradeniz/GroupProject1/gsva_scores.csv")

gsva.es.t <- as.data.frame(t(gsva.es))

df_up <- data.frame(samples= rownames(gsva.es.t), up=gsva.es.t[,1])
df_up_sorted <- df_up[order(df_up$up, decreasing = TRUE),]

df_down <- data.frame(samples= rownames(gsva.es.t), down=gsva.es.t[,2])
df_down_sorted <- df_down[order(df_down$down, decreasing = TRUE),]

# compute quartiles
quantUp <- quantile(df_up_sorted$up)
quantDown <- quantile(df_down_sorted$down)


# pathway analyse
pathway <- data.frame(samples=df_up_sorted[,1], HER2=vector(mode="character", nrow(df_up_sorted)))

pathwayActivity <- lapply(pathway$samples, function(sample) {
        if (df_up_sorted[df_up_sorted$samples==sample, 2] > quantUp[length(quantUp)-1] && df_down_sorted[df_down_sorted$samples==sample, 2] < quantDown[2]) {
                "+"
        } else if (df_up_sorted[df_up_sorted$samples==sample, 2] < quantUp[2] && df_down_sorted[df_down_sorted$samples==sample, 2] > quantDown[length(quantUp)-1]){
                "-"
        } else {
                "?"
        }
})


pathwayActivity <- factor(pathwayActivity, levels=c("+", "-", "?"))

pathway$HER2 <- pathwayActivity

summary(pathway$HER2)
pathway

for (i in 1:nrow(pathway)) {
        pathway[i,1] <- gsub("\\.", "-", pathway[i,1])
}

labelDF <- read.table("Sample_lable.txt", header=TRUE,sep = "\t")

read_labels=function(labelDataFile){
        labelDF <- read.table(labelDataFile, header=TRUE,sep = "\t") 
        return(labelDF)
}

labelDF <- read_labels("Sample_lable.txt")

compare_classifications = function(true_classes_df, predicted_classes_df){
        # create confusion matrix
        confusion_matrix_HER2<-matrix(0,nrow = 3,ncol = 3,
                                      dimnames = list(c("Actual_Pos","Actual_Neg","Actual_Uncertain"),
                                                      c("Prediction_Pos","Prediction_Neg", "Prediction_Uncertain")))
        df<-merge(x=true_classes_df, y=predicted_classes_df, by.x = "CaseID", by.y="samples")
        
        true_positive <- nrow(df[(df$HER2.y == '+' & df$HER2.x == "Positive"),]) # predicted positive, actual positive
        false_positive <- nrow(df[(df$HER2.y == '+' & df$HER2.x == "Negative"),]) # predicted positive, actual negative
        confusion_matrix_HER2[1,1] <- true_positive
        confusion_matrix_HER2[2,1] <- false_positive
        confusion_matrix_HER2[3,1] <- nrow(df[df$HER2.y == '+',])-true_positive-false_positive
        
        true_negative <- nrow(df[(df$HER2.y == '-' & df$HER2.x == "Negative"),]) # predicted negative, actual negative
        false_negative <- nrow(df[(df$HER2.y == '-' & df$HER2.x == "Positive"),]) # predicted negative, actual positive
        confusion_matrix_HER2[2,2] <- true_negative
        confusion_matrix_HER2[1,2] <- false_negative
        confusion_matrix_HER2[3,2] <- nrow(df[df$HER2.y == '-',])-true_negative-false_negative
        
        true_uncertain <- nrow(df[(df$HER2.y == '?' & df$HER2.x != "Negative" & df$HER2.x != "Positive"),]) # predicted uncertain, actual uncertain
        confusion_matrix_HER2[3,3] <- true_uncertain
        confusion_matrix_HER2[1,3] <- nrow(df[(df$HER2.y == '?' & df$HER2.x == "Positive"),])  # predicted uncertain, actual positive
        confusion_matrix_HER2[2,3] <- nrow(df[(df$HER2.y == '?' & df$HER2.x == "Negative"),])  # predicted uncertain, actual negative
        
        accuracy <- (true_positive + true_negative + true_uncertain) / nrow(df) * 100
        print(paste0("accuracy: ", format(round(accuracy, 2), nsmall = 2)))
        print(paste0("total number samples: " , nrow(df)))
        print(paste0("true_positive: " , true_positive))
        print(paste0("true_negative: " , true_negative))
        print(paste0("true_uncertain: " , true_uncertain))
        print(paste0("sum in confusion_matrix_HER2: " , sum(confusion_matrix_HER2)))
        print(confusion_matrix_HER2)
} 

compare_classifications(labelDF,pathway)

































