# Clear workspace
rm(list=ls())

# Installing
#install.packages("BiocManager")
#BiocManager::install("GSVA")
#install.packages("readr")

# Loading
library("readr")
library(GSVA)
install.packages("~/git_project/PathAnalyser_0.0.0.9000.tar.gz", repos = NULL, type="source")
library("PathAnalyser")


data_se <- read.table("~/GroupProject1/TCGA_unannotated.txt", header=TRUE, sep = "\t", row.names=1)
data_se <- data_se[,1:200]
up <- read_lines("~/GroupProject1/ERBB2_UP.V1_UP.grp", skip = 3, n_max = -1L)
dn <- read_lines("~/GroupProject1/ERBB2_UP.V1_DN.grp", skip = 3, n_max = -1L)
sig_df <- data.frame(c(dn,up), expression= c(rep(-1, length(dn)), rep(1, length(up))))

readsig(up, dn)

"%notin%" <- Negate("%in%")
data_se_not_in_sig_df <- data_se[rownames(data_se) %notin% sig_df[,1],]
data_se_in_sig_df <- data_se[rownames(data_se) %in% sig_df[,1],]
data_se_subset <- rbind(data_se_not_in_sig_df[1:300,],data_se_in_sig_df)


data_se_subset<-readmatrix(data_se_subset)
evaluatematrix(data_se_subset)

normalized_se <- transform(data_se[])
data_se<- check_signature_vs_dataset(normalized_se, as.data.frame(up),as.data.frame(dn))

visualise_GSVA(sig_df, data_se)

classes_df <-classify(sig_df, data_se)

labelDF <- read.table("~/GroupProject1/Sample_lable.txt", header=TRUE,sep = "\t")

confusion_matrix_HER2<-calculate_accuracy(labelDF, classes_df)

classes_pca(normalized_se,classes_df, pathway_name = "ER")

