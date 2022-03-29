#' Reading gene expression data from file
#' 
#' @description  Reads gene expression data file from user,
#' removes Nas, if the first column consist of gene names,
#' gives the genenames to the rownames of the matrix and deletes
#' the first column. Before converting the first column to rownames
#' it checks for any duplicated gene present in the first column.
#' Following that the matrix is converted to a numeric data matrix.
#' It also returns a heatmap which helps the user
#' to identify duplicate samples (columns) from the
#' heatmap, which in turn will help the user to remove duplicated 
#' samples according to their choice.
#'
#' @author Tania Pal \email{taniya.pal.094@cranfield.ac.uk}
#' @param matrix  Full path of the gene expression data file name 
#'
#' @return Structured matrix containing gene symbols/IDs as rownames
#' Also returns a heatmap with identifiable duplicated samples for 
#' the convinience of the user.
#' and the Sample IDs as colnames.
#' @export
#'
#' @examples read_expression_data("/Users/taniyapal/Documents/Group Project/Matrix used for GSVA Analysis by Taniya.csv"))

read_input_file<- function(file_name){
  
  #getting the delimiter for the file whether it is "\t" or "," or " "
  #install.packages("reader")
  library(reader)
  delimiter=get.delim(file_name)
  input=read.delim(file_name, sep=delimiter)
  input=na.omit(input)
  #removing duplicated gene symbols from first column
  
  input=input[!duplicated(input[,1]),]
  
  #giving the gene symbols of the first column to rownames
  
  rownames(input)=input[,1]
  input=input[,-1]
  
  
  #removing NAs from gene expression values

  input=data.matrix(input)
  return(input)
 
  ###Deal with duplicated samples
  library(reshape2)
  data=melt(cor(input))
  ggplot(data=data, aes(x=Var1, y=Var2)) +
    geom_tile(aes(fill=value)) + 
    scale_fill_gradient2(limits=c(-1, 1)) + 
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Heatmap showing duplicated samples in gene expression matrix") + xlab(NULL) + ylab(NULL)
  
  
  
  
}


#' Reads up-regulated and down-regulated gene signatures from files,
#' formats it according to the needs of the package analysis,
#' and combines the up and down regulated signature files
#' into a single data frame
#' 
#' @description Reads up and down regulated signature files from user 
#' and structures it according to package requirements. It returns a 
#' dataframe. The first column of the dataframe  lists together
#' both up and down regulated signatures. The second column
#' signifies whether they are up(+1) or down(-1) regulated.
#' 
#' @author Tania Pal \email{taniya.pal.094@cranfield.ac.uk}
#' @param  up_sig   Up-regulated gene-set
#' @param  down_sig  Down-regulated gene-set
#'
#' @return A dataframe containing both up regulated and down 
#' regulated signature files and signifying their up or down 
#' expression with +1 and -1 respetively.
#' @export
#'
#' @examples read_signature_data("ESR1_UP.v1_UP.csv","ESR1_DN.v1_DN.csv" )

read_signature_data=function(up_sig_file, down_sig_file){
  
  #loading the required package
  require("readr")
  
  #reading the up regualted gene signature file 
  up <- read_lines(up_sig_file, skip = 3, n_max = -1L)
  
  #reading the down regulated gene signature file
  dn <- read_lines(down_sig_file, skip = 3, n_max = -1L)
  
  #vector combining both up and down regulated signatures
  list=c(up,dn)
  
  #combining up and downregulated signatures in
  #a column of datraframe
  sig_df <- data.frame("Signatures"=list, "Symbols representing expression"=c(rep(1, length(up)), rep(-1, length(dn))))
  
  return(sig_df)
}

#' Transformation of gene expression data with log cpm
#' transformation method
#' @description Plots the raw data before transformation in the form 
#' of a boxplot. Then, it normalizes the raw counts of gene
#' expression with the log cpm transformation
#' method and returns a boxplot of the gene expression 
#' matrix after transformation for sanity check of the
#' transformation. The user can check the distribution of 
#' the gene expression values after log cpm
#' transformation with the help  of a box plot.
#'
#' @author Tania Pal \email{taniya.pal.094@cranfield.ac.uk}
#' @param Gene expression matrix after being
#' pre processed by read_input_file function
#'
#' @return Two boxplots, one representing
#' distribution of gene expression values before transformation,
#' another one after transformation.
#' @export
#'
#' @examples transform_matrix(input)

transform_matrix=function(input){
  
  
  par(mfrow=c(2,1))
  #boxplot before transformation
  plot.before.transformation=boxplot(log(input+0.5), main="Plot before normalization", axes=F)
  
  library(edgeR)
  
  #counts per million (cpm) transformation of raw gene expression
  #values of matrix
  cpm.values=cpm(input)
  
  #Taking the log of the cpm transformed values
  #(log-cpm transformation)
  input=log(cpm.values)
  
  #boxplot after transformation
  plot.after.transformation=boxplot(log(cpm.values+0.5), main="Plot after normalization", axes=F)
 
  
 
}

