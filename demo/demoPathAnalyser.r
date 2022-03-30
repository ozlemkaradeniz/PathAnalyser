# Clear workspace
rm(list=ls())

install.packages("~/git_project/PathAnalyser_0.0.0.9000.tar.gz", repos = NULL, type="source")
library("PathAnalyser")

gene_expression_file <- "~/GroupProject/TCGA_unannotated.txt"
expression_mtx <- read_expression_data(gene_expression_file)

# 200 samples from the expression set is used due to the performance issues
expression_mtx <- expression_mtx[,1:200]

up_signature_file <- "~/GroupProject/ERBB2_UP.V1_UP.grp"
down_regulated_file <- "~/GroupProject/ERBB2_UP.V1_DN.grp"
signature_df <- read_signature_data(up_signature_file, up_signature_file)

#Expression gene sub-matrix is created below in order to overcome performance issues coming with the whole gene expression data matrix
#sub-matrix contains all genes in signature dataframe plus additional 300 genes
"%notin%" <- Negate ("%in%")
expression_mtx_not_in_sig_df <- expression_mtx[rownames(expression_mtx) %notin% signature_df[,1],]
expression_mtx_in_sig_df <- expression_mtx[rownames(expression_mtx) %in% signature_df[,1],]
expression_mtx_subset <- rbind(expression_mtx_not_in_sig_df[1:300,],expression_mtx_in_sig_df)

# Getting statistical information for the gene expression data matrix
evaluatematrix(expression_mtx_subset)

# cpm normalization is performed on the ene expression data matrix
normalized_expression_mtx <- transform(expression_mtx)
upregulated_genes_df <- signature_df[signature_df$expression == 1, ]
downregulated_genes_df <- signature_df[signature_df$expression == -1, ]
data_se<- check_signature_vs_dataset(normalized_expression_mtx, upregulated_genes_df, downregulated_genes_df)

visualise_GSVA(signature_df, expression_mtx)

classes_df <-classify(signature_df, expression_mtx)

labelDF <- read.table("~/GroupProject/Sample_lable.txt", header=TRUE,sep = "\t")

confusion_matrix_HER2<-calculate_accuracy(labelDF, classes_df)


classes_pca(normalized_se,classes_df, pathway_name = "ER")

