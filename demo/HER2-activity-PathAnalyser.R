# Demo R script for using PathAnalyser functionality
#
# In this demo we will walk through the steps of classifying samples in breast
# cancer RNAseq data set by HER2 pathway activity.
# This script follows the example provided in the vignette, for further
# information please consult the vignette.
#
# Expression data information
#---------------------------------
# The gene expression data set featured in this demo is a subset of RNA-Seq data
# for 20,124 genes for 60 primary breast cancer (biopsy) samples obtained from
# TCGA online. This subset was taken from an initial data set of 1,100 primary
# breast tumor samples from TCGA.
#
# Gene signature information
#----------------------------
# The gene signature data files were obtained from the Sensitive to endocrine
# therapy genome index proposed by Symmans et al. (2011), which consists of a
# list of genes co-expressed with HER2 pathway.

# Installing package using devtools
#install.packages("devtools") # un-comment this line to install the devtools package
#devtools::install_github("ozlemkaradeniz/PathAnalyser") # un-comment this line to install the package

# Load the library into R
library(PathAnalyser)

#-------------------------------------------------------------------------------
# 1. Reading input gene expression and signature files
#-------------------------------------------------------------------------------
# read expression matrix data
data_mat <- read_expression_data("inst/extdata/HER2_toydata_RNAseq.txt")
# Column names represent sample names while row names represent gene names
head(data_mat)
# there are 20,124 genes for 60 samples
dim(data_mat)

# read the up-regulated and down-regulated gene signature files
sig_df <- read_signature("inst/extdata/SMID_BREAST_CANCER_ERBB2_UP.grp",
                         "inst/extdata/SMID_BREAST_CANCER_ERBB2_DN.grp"
  )
# first column represents gene names and the second column represents their
# expression value (1 for up-regulated expression when HER2 pathway is active and
# -1 for down-regulated expression when HER2 pathway is active)
# up-regulated genes are at the top of signature data frame
head(sig_df)
# down-regulated genes of the gene signature are at the bottom of the gene
# signature data frame
tail(sig_df)

# there are 156 genes in this HER2 signature obtained from the SET index
dim(sig_df)

#-------------------------------------------------------------------------------
# 2. QC and data pre-processing
#-------------------------------------------------------------------------------
# logCPM normalise RNAseq raw counts of expression data set
norm_data <- log_cpm_transform(data_mat)
# check the gene signature and gene signature have consistent genes
# filters out genes from expression data set that are not present in the gene
# signature and those that are not expressed in at least 10% of samples.
norm_data <- check_signature_vs_dataset(norm_data, sig_df)

#-------------------------------------------------------------------------------
# 3. Classification of samples based on pathway activity
#-------------------------------------------------------------------------------
# Samples can be classified according to a percentile threshold or an absolute
# threshold for GSVA scores

# Using a percentile threshold (default = 25% so quartile threshold essentially)
classes_df.perc25 <- classify_GSVA_percent(norm_data, sig_df)
# using a percentile threshold of 50%
# (At this percentile the number of uncertain classifications are reduced to
# their minimum, as only those samples that have consistent expression with only
# the up-gene set signature or only the down-gene set of the signature are
# classified as uncertain)
classes_df.perc50 <- classify_GSVA_percent(norm_data, sig_df,
                                           percent_thresh = 50)

#-------------------------------------------------------------------------------
# 4. Visualise sample classification
#-------------------------------------------------------------------------------
# To generate a PCA plot showing clustering of samples based on pathway
# classified labels:
classes_pca(norm_data, classes_df.perc25, pathway = "HER2")
# PCA plot for 50th percentile classification threshold for our data
classes_pca(norm_data, classes_df.perc50, pathway = "HER2")

#-------------------------------------------------------------------------------
# 6. Classification evaluation (optional step if true pathway activity classes
#    are available)
#-------------------------------------------------------------------------------
# To generate a confusion matrix for actual classes vs predicted classes:
confusion_mat.perc25 <- calculate_accuracy("inst/extdata/Sample_labels.txt",
                                           classes_df.perc25, pathway = "HER2")

# for more detailed classification evaluation metric info use the optional
# parameter: display_stats=TRUE
confusion_mat.perc25 <- calculate_accuracy("inst/extdata/Sample_labels.txt",
                                           classes_df.perc25, pathway="HER2",
                                           show_stats=T)

# classification evaluation metrics for 50th percentile threshold classification
confusion_mat.perc50 <- calculate_accuracy("inst/extdata/Sample_labels.txt",
                                           classes_df.perc50, pathway = "HER2")

# for more detailed classification evaluation metric info use the optional
# parameter: display_stats=TRUE
confusion_mat.perc50 <- calculate_accuracy("inst/extdata/Sample_labels.txt",
                                           classes_df.perc50, pathway="HER2",
                                           show_stats=T)

