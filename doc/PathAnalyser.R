## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)

## ---- echo=F, fig.cap="PathAnalyser workflow"---------------------------------
knitr::include_graphics("pathway_workflow.png")

## ----eval=F-------------------------------------------------------------------
#  # If not already installed
#  install.packages("BiocManager")
#  
#  BiocManager::install(c("GSVA", "pROC", "edgeR", "reshape2", "ggplot2","limma",
#                         "VennDiagram", "NCmisc", "futile.logger"),
#                       dependencies = TRUE)
#  

## ----eval=F-------------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("ozlemkaradeniz/PathAnalyser")

## ----eval=FALSE, message=FALSE, warning=FALSE, cache=FALSE--------------------
#  library(utils)
#  install.packages("~/PathAnalyser", repo=NULL, type="source")

## ----message=FALSE, warning=FALSE, cache=FALSE--------------------------------
library(PathAnalyser)

## ----fig.cap="An example gene expression matrix text file"--------------------
knitr::include_graphics("expr_dataset_example.png")

## ----fig.cap="An example gene set file format containing either the up- or down-regulated gene sets of the gene signature"----
knitr::include_graphics("up_gene_signature.png")

## -----------------------------------------------------------------------------
data("ER_sig_df")
data("HER2_sig_df")

## -----------------------------------------------------------------------------
dim(ER_sig_df)
head(ER_sig_df)
# up-regulated genes are given an expression value of 1
ER_up_sig <- ER_sig_df[ER_sig_df$expression == 1 ,]
dim(ER_up_sig)
# down-regulated genes are given an expression value of -1
ER_dn_sig <- ER_sig_df[ER_sig_df$expression == -1 ,]
dim(ER_dn_sig)
# for more details
?ER_sig_df

## -----------------------------------------------------------------------------
dim(HER2_sig_df)
head(HER2_sig_df)
# up-regulated genes are given an expression value of 1
HER2_up_sig <- HER2_sig_df[HER2_sig_df$expression == 1 ,]
dim(HER2_up_sig)
# Down-regulated genes are given an expression value of -1
HER2_dn_sig <- HER2_sig_df[HER2_sig_df$expression == -1 ,]
dim(HER2_dn_sig)
# for more details
?HER2_sig_df

## -----------------------------------------------------------------------------
data("ER_data_se1")
data("HER2_data_se1")
# column names represent case (sample) IDs from TCGA

## -----------------------------------------------------------------------------
dim(ER_data_se1)
# Expression data for first 6 genes
head(ER_data_se1)

## -----------------------------------------------------------------------------
dim(HER2_data_se1)
# Expression data for first 6 genes
head(HER2_data_se1)

## ----eval=FALSE---------------------------------------------------------------
#  data_se <- read_expression_data("/Users/taniyapal/Documents/Group Project/TCGA_unannotated.txt")
#  data_se<-data_se[,1:200]
#  print(data_se[1:5, 1:5])
#  
#  

## ----eval=FALSE---------------------------------------------------------------
#  sig_df <- read_signature_data("/Users/taniyapal/Documents/Group Project/ESR1_UP.v1._UP.csv", "/Users/taniyapal/Documents/Group Project/ESR1_DN.v1_DN.csv")
#  head(sig_df)

## ----fig.wide=TRUE------------------------------------------------------------
# using toy data set as expression matrix
data("ER_data_se1")
data_se <- ER_data_se1 
normalized_se <- log_cpm_transformation(data_se)

## ----fig.wide=TRUE------------------------------------------------------------
normalized_se <- check_signature_vs_dataset(normalized_se, ER_sig_df)

## ---- echo=FALSE, fig.wide=TRUE, fig.cap="Overiew of algorithm diagram"-------
knitr::include_graphics("algorithm_diagram.png")

## ----plot, fig.cap="Density plots showing the distribution of GSVA scores for samples in logCPM-normalised `ER_data_se1` data set for the up-regulated (left) and down-regulated gene set (right) of the pathway signature using `gsva_scores_dist` function.", fig.wide=TRUE----
gsva_scores_dist(ER_sig_df, normalized_se)

## ----thresholds, fig.cap="GSVA scores distributions of the samples for the down-regulated (left) and up-regulated gene set (right) of the pathway gene signature, showing absolute thresholds of -0.2 and 0.2 for characterising low and high expression abundance of both gene sets respectively. The expression data set used is the logCPM-normalised `ER_data_se1` gene expression data set.", fig.wide=TRUE----
plot <- gsva_scores_dist(ER_sig_df, normalized_se)
# Add thresholds on plot
library(ggplot2)
data_threshs <- data.frame(Geneset=c("Up", "Down"), vline=c(-0.2, 0.2))
plot + geom_vline(xintercept=data_threshs$vline, linetype=2)

## ----relaxed, fig.cap="GSVA scores distributions of samples with more relaxed thresholds for assessing expression consistency of the down-regulated gene set (left) and the up-regulated gene set (right) with the pathway gene expression signature. Low and high expression abundance thresholds for the down-regulated gene set (left) are set to -0.1 and 0.1 respectively, while low and high expression abundance thresholds for the up-regulated gene set (right) are also set to -0.1 and 0.1 respectively. Data is generated from running GSVA on logCPM-normalised `ER_data_se1` expression dataset with ER signature (`ER_sig_df`).", fig.wide=TRUE----
# more relaxed thresholds, fewer uncertain labels
data_threshs <- data.frame(Geneset=c("Up", "Down"), 
                           vline=c(-0.1, 0.1))
plot + geom_vline(xintercept=data_threshs$vline, linetype=2)

## ----stringent, fig.cap="Density plots for GSVA scores distributions for samples with relatively more stringent thresholds for assessing consistency of the up-regulated gene set (left) and down-rgulated gene set (right) with the pathway signature. The low and high expression abundance thresholds are -0.4 and 0.4 for both up-regulated and down-regulated gene sets respectively. Data is generated from running GSVA on logCPM-normalised `ER_data_se1` expression dataset with ER signature (`ER_sig_df`).", fig.wide=TRUE----
# more stringent thresholds, greater uncertain labels
data_threshs <- data.frame(Geneset=c("Up", "Down"), vline=c(-0.4, 0.4))
plot + geom_vline(xintercept=data_threshs$vline, linetype=2)

## -----------------------------------------------------------------------------
classes_df <- classify_GSVA_abs(ER_sig_df, normalized_se, 
                                    up_thresh.low=-0.2, up_thresh.high=0.2, 
                                    dn_thresh.low=-0.2, dn_thresh.high=0.2)

## -----------------------------------------------------------------------------
# default percentile = 25% (quartile)
classes_df.percent <- classify_GSVA_percent(ER_sig_df, normalized_se)
# custom percentile = 30%
classes_df.percent <- classify_GSVA_percent(ER_sig_df, normalized_se, 
                                            percent_thresh=30)

## -----------------------------------------------------------------------------
classes_pca(normalized_se, classes_df, pathway_name = "ER")

## ----fig.wide=TRUE------------------------------------------------------------
true_labels_df <- read.table("Sample_labels.txt", sep="\t", header=T)
confusion_matrix <- calculate_accuracy(true_labels_df, classes_df, 
                                       pathway = "ER", display_statistics=TRUE, 
                                       display_roc_curve=TRUE)

## -----------------------------------------------------------------------------
# Load transcriptomic data set (gene expression matrix of samples)
data_se <- read_expression_data("../raw_data/toy_data.txt")
dim(data_se)
head(data_se)
# read signature data from the two individual gene set files for up-regulated
# and down-regulated gene sets
ER_sig <- read_signature_data("../raw_data/ESR1_UP.v1._UP.grp", "../raw_data/ESR1_DN.v1_DN.grp")
dim(ER_sig)

## ----fig.wide=T---------------------------------------------------------------
# Transforming matrix with log cpm transformation and sanity check of the transformation
normalized_se <- log_cpm_transformation(data_se)

## ----fig.wide=TRUE------------------------------------------------------------
normalized_se <- check_signature_vs_dataset(normalized_se, ER_sig)
dim(normalized_se)

## ----fig.wide=TRUE------------------------------------------------------------
gsva_scores_dist(ER_sig, normalized_se)


## -----------------------------------------------------------------------------
classes_df <- classify_GSVA_percent(ER_sig, normalized_se)
head(classes_df)

## -----------------------------------------------------------------------------
classes_df_50p <- classify_GSVA_percent(ER_sig, normalized_se, 
                                        percent_thresh=50)

## -----------------------------------------------------------------------------
classes_df_abs <- classify_GSVA_abs(ER_sig, normalized_se, up_thresh.low = -0.2, 
                                    up_thresh.high = 0.2, dn_thresh.low = -0.2, 
                                    dn_thresh.high = 0.2)

## -----------------------------------------------------------------------------
classes_pca(normalized_se, classes_df, pathway_name = "ER")

## ----fig.wide=TRUE------------------------------------------------------------
# read true pathway activity class labels for data set samples
true_labels_df <- read.table("Sample_labels.txt", 
                             header=TRUE, sep = "\t")
# assess accuracy and generate confusion matrix for classification
confusion_matrix <- calculate_accuracy(true_labels_df, classes_df, 
                                       pathway = "ER")
# detail breakdown of classification evaluation (accuracy, % classified etc)
confusion_matrix <- calculate_accuracy(true_labels_df, classes_df, 
                                       pathway = "ER", display_statistics = T)

## ----fig.wide=TRUE------------------------------------------------------------
# detail breakdown of classification evaluation (accuracy, % classified etc)
confusion_matrix <- calculate_accuracy(true_labels_df, classes_df, 
                                       pathway = "ER", display_statistics = T,
                                       display_roc_curve = T)

## ----cache=FALSE--------------------------------------------------------------
sessionInfo()

