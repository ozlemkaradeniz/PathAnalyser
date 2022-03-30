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
#' @author Taniya Pal \email{taniya.pal.094@cranfield.ac.uk}
#' 
#' @param filename  Full path of the gene expression data file name 
#'
#' @return Structured matrix containing gene symbols/IDs as rownames
#' Also returns a heatmap with identifiable duplicated samples for 
#' the convinience of the user.
#' and the Sample IDs as colnames.
#' @export
#'
#' @examples read_expression_data("/Users/taniyapal/Documents/Group Project/Matrix used for GSVA Analysis by Taniya.csv"))

read_input_file<- function(file_name){
  
  #loading the required packages
  library(reader)
  
  #getting the delimiter for the file whether it is "\t" or "," or " "
  delimiter=get.delim(file_name)
  
  #reading the file provided by the user
  input=read.delim(file_name, sep=delimiter)
  
  #removing the NAs 
  input=na.omit(input)
  
  #removing duplicated gene symbols from first column
  input=input[!duplicated(input[,1]),]
  
  #giving the gene symbols of the first column to rownames
  if (typeof(input[,1])=="character")
    rownames(input)=input[,1]
    input=input[,-1]
  
  
  #converting the data frame in numeric matrix
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


#' Reads up-regulated and down-regulated gene signatures from files
#' 
#' @description Reads up and down regulated signature files from user 
#' and structures it according to package requirements. It returns a 
#' dataframe. The first column of the dataframe  lists together
#' both up and down regulated signatures. The second column
#' signifies whether they are up(+1) or down(-1) regulated.
#' 
#' @author Taniya Pal \email{taniya.pal.094@cranfield.ac.uk}
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
  
  #reading the up regulated gene signature file 
  up_sig <- read_lines(up_sig_file, skip = 3, n_max = -1L)
  
  #reading the down regulated gene signature file
  dn_sig<- read_lines(down_sig_file, skip = 3, n_max = -1L)
  
  #vector combining both up and down regulated signatures
  list=c(up_sig,dn_sig)
  
  #combining up and down regulated signatures in
  #a column of dataframe
  sig_df <- data.frame("Signatures"=list, "Symbols representing expression"=c(rep(1, length(up_sig)), rep(-1, length(dn_sig))))
  
  return(sig_df)
}

#' Transformation of gene expression data with log cpm 
#' transformation method
#' 
#' @description Plots the raw data before transformation in the form 
#' of a boxplot. Then, it normalizes the raw counts of gene
#' expression with the log cpm transformation
#' method and returns a boxplot of the gene expression 
#' matrix after transformation for sanity check of the
#' transformation. The user can check the distribution of 
#' the gene expression values after log cpm
#' transformation with the help  of a box plot.
#'
#' @author Taniya Pal \email{taniya.pal.094@cranfield.ac.uk}
#' 
#' @param input expression matrix after being
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


#'Checking overlap between matrix and signature gene symbols
#'
#' @description Plots a Venn diagram checking which gene symbols 
#' overlap between the signature file
#' and the expression matrix file.
#' The user gets an estimates of how many gene symbols overlap
#'  between the two files.
#'
#' @author Taniya Pal \email{taniya.pal.094@cranfield.ac.uk}
#'
#' @param input The expression matrix after being processed and 
#' returned by read_input_file function
#' @param sig_df The signature data frame after being
#' processed and returned by read_signature_data function
#' 
#' @return A Venn diagram representing the overlap of 
#' gene symbols between the two files
#'
#' @export
#' @examples check_matrix_sig_overlap(input, sig_df)


check_matrix_sig_overlap=function(input, sig_df){
  #loading the required library
  library(VennDiagram)
  
  #listing the gene symbols in expression matrix
  rownames=as.list(rownames(expdata))
  
  #listing the gene symbols in signature
  signatures=as.list(sig$Signatures)
  
  #finding the length of intersection between the two lists
  intersection=length(intersect(rownames, signatures))
  
  #Setting the color
  myCol=c("green", "yellow")
  
  draw.pairwise.venn(
    area1 = length(rownames), area2=length(signatures),
    category = c("Expression Matrix" , "Signature"),
    main = 'Checking overlap between matrix and gene signature symbols.png', lwd = 2,
    lty = 'blank',
    fill = myCol, cat.cex = 0.7,
    cat.fontface = c("bold", "bold"),
    cat.default.pos = "text",
    cat.pos = c(-27, 27),
    cat.dist = c(-0.055, 0.055),
    cat.fontfamily = c("sans", "sans"), cross.area=94, scaled=F)
}
