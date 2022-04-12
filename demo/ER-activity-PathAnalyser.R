# Demo R script for using PathAnalyser functionality
#
# In this demo we will walk through the steps of classifying samples in breast
# cancer RNAseq data set by ER pathway activity.
# This script follows the example provided in the vignette, for further
# information please consult the vignette.
#
# Expression data information
#---------------------------------
# The gene expression data set featured in this demo is a subset of RNA-Seq data
# for 20,124 genes for 20 primary breast cancer (biopsy) samples obtained from
# TCGA online. This subset was taken from an initial data set of 1,100 primary
# breast tumor samples from TCGA.
#
# Gene signature information
#----------------------------
# The gene signature data files were obtained from the Sensitive to endocrine
# therapy genome index proposed by Symmans et al. (2011), which consists of a
# list of genes co-expressed with ER pathway.

# Installing package using devtools
#install.packages("devtools") # un-comment this line to install the devtools package
#devtools::install_github("ozlemkaradeniz/PathAnalyser") # un-comment this line to install the package

# Load the library into R
library(PathAnalyser)

#-------------------------------------------------------------------------------
# 1. Reading input gene expression and signature files
#-------------------------------------------------------------------------------
# read expression matrix data
data_se <- read_expression_data("raw_data/toy_data.txt")
# Column names represent sample names while row names represent gene names
head(data_se)
# there are 20,124 genes for 20 samples
dim(data_se)

# read the up-regulated and down-regulated gene signature files
sig_df <- read_signature_data("raw_data/ESR1_UP.v1._UP.grp",
                              "raw_data/ESR1_DN.v1_DN.grp")
# first column represents gene names and the second column represents their
# expression value (1 for up-regulated expression when ER pathway is active and
# -1 for down-regulated expression when ER pathway is active)
head(sig_df) # up-regulated genes are at the top of signature data frame
tail(sig_df) # down-regulated genes of the gene signature are at the bottom of
# the gene signature data frame

# there are 160 genes in this ER signature obtained from the SET index
dim(sig_df)

#-------------------------------------------------------------------------------
# 2. QC and data pre-processing
#-------------------------------------------------------------------------------
# logCPM normalise RNAseq raw counts of expression data set
normalized_se <- log_cpm_transformation(data_se)
# check the gene signature and gene signature have consistent genes
# filters out genes from expression data set that are not present in the gene
# signature and those that are not expressed in at least 10% of samples.
normalized_se <- check_signature_vs_dataset(normalized_se, sig_df)

#-------------------------------------------------------------------------------
# 3. Classification of samples based on pathway activity
#-------------------------------------------------------------------------------
# Samples can be classified according to a percentile threshold or an absolute
# threshold for GSVA scores

# Using a percentile threshold (default = 25% so quartile threshold essentially)
classes_df.perc <- classify_GSVA_percent(sig_df, normalized_se)
# using a percentile threshold of 50%
# (At this percentile the number of uncertain classifications are reduced to
# their minimum, as only those samples that have consistent expression with only
# the up-gene set signature or only the down-gene set of the signature are
# classified as uncertain)
classes_df.perc50 <- classify_GSVA_percent(sig_df, normalized_se,
                                          percent_thresh = 50)

#-------------------------------------------------------------------------------
# 4. Visualise sample classification
#-------------------------------------------------------------------------------
# To generate a PCA plot showing clustering of samples based on pathway
# classified labels:
classes_pca(normalized_se, classes_df.perc, pathway_name = "ER")
# PCA plot for 50th percentile classification threshold for our data
classes_pca(normalized_se, classes_df.perc50, pathway_name = "ER")

#-------------------------------------------------------------------------------
# 6. Classification evaluation (optional step if true pathway activity classes
#    are available)
#-------------------------------------------------------------------------------
# To generate a confusion matrix for actual classes vs predicted classes:
confusion_mat <- calculate_accuracy("raw_data/Sample_labels.txt", classes_df.perc, pathway="ER")

# for more detailed classification evaluation metric info and roc curve diagram
# use the optional parameter: display_stats=TRUE and display_roc_curve=TRUE
confusion_mat <- calculate_accuracy("raw_data/Sample_labels.txt", classes_df.perc, pathway="ER",
                                    display_statistics=T, display_roc_curve=T)

# classification evaluation metrics for 50th percentile threshold classification
confusion_mat <- calculate_accuracy("raw_data/Sample_labels.txt", classes_df.perc50, pathway="ER")

# for more detailed classification evaluation metric info and roc curve diagram
# use the optional parameter: display_stats=TRUE and display_roc_curve=TRUE
confusion_mat <- calculate_accuracy("raw_data/Sample_labels.txt", classes_df.perc50, pathway="ER",
                                    display_statistics=T, display_roc_curve=T)

