#' Reading gene expression data from file
#' @description  Reads gene expression data file, performs validation checks and converts to numeric data matrix
#'
#' @author Tania Pal \email{taniya.pal.094@cranfield.ac.uk}
#' @param matrix  Gene expression data file
#'
#' @return Structured matrix for the convenience of the package
#' @export
#'
#' @examples read_expression_data(file_name)

read_expression_data <- function(file_name){
  # #rownames(expdata)=expdata[,1]
  # expdata=expdata[-1,-1]
  # for (i in 1:ncol(expdata)){
  #   class(expdata[,i])="numeric"
  # }
  #
  # for(i in 1:ncol(expdata)){
  #   expdata[,i][is.na(expdata[,i])]<-mean(expdata[,i],na.rm=TRUE)
  # }
  # expdata=as.matrix(expdata)
  # return(expdata)
  file_name <- "~/GroupProject/TCGA_unannotated.txt"
  raw.df <- read.table(file_name, header=TRUE, sep = "\t")

  samples.df <- read.table(file_name, nrows=1)
  samples <- unlist(samples.df[1,])

  # Deal with duplicated samples (colnames)
  # if there is any duplicated sampled, sample name is renamed to be unique
  if(any(duplicated(samples))){
    print(paste("File",
                basename(file_name),
                "has duplicated samples. Duplicated sample names have been changed to unique names."))
    samples <- make.unique(samples)
  }

  # Deal with duplicated genes
  # if there is any duplicated gene, error message is raised and the program stops
  raw.df <- read.table(file_name, skip=1)
  genes <- raw.df[,1]

  # Check for duplicates
  if(any(duplicated(genes))){
    stop(paste("File",
               basename(file_name),
               "has duplicated genes"))
  }

  # Gene expression dataframe is created with sample names and gene names
  data.df <- raw.df[,-1]
  rownames(data.df) <- genes
  colnames(data.df) <- samples

  # dataframe is converted to numeric data matrix
  data.mx <- as.matrix(data.df)

  #error is raised and program stops if any problem occurs during coercion into numeric data matrix
  expr.mx <- tryCatch(matrix(as.numeric(data.mx),nrow=nrow(data.mx)),
                      warning=function(w){return="Fail"},
                      error=function(e){return="Fail"})

  # Note use of identical() instead of equality check
  if(identical(expr.mx,"Fail")){stop(paste("Can't convert some elements of",
                                           basename(file_name),
                                           "to numeric values"))}


  # if the numbers of NAs are not equal in dataframe and numeric data matrix, error is raised and program stops
  if( sum(is.na(data.df)) != sum(is.na(expr.mx)) ){
    stop(paste("Number of NAs in numeric data matrix is different from the original dataframe!,
               there might be error during conversion"))
  }

  rownames(expr.mx) <- rownames(data.mx)
  colnames(expr.mx) <- colnames(data.mx)

  return(expr.mx)
}

#' Reading up-regulated and down-regulated gene signatures from files
#' @description Reads up and down regulated signature files from user and structures it according to package requirements
#'
#' @author Tania Pal \email{taniya.pal.094@cranfield.ac.uk}
#' @param  up_sig   Up-regulated gene-set
#' @param  down_sig  Down-regulated gene-set
#'
#' @return Up and down regulated signature files in the form a list for the convenience of  the package
#' @export
#'
#' @examples read_signature_data(up_sig_file, down_sig_file)

read_signature_data=function(up_sig_file, down_sig_file){
  require("readr")
  up <- read_lines(up_sig_file, skip = 3, n_max = -1L)
  dn <- read_lines(down_sig_file, skip = 3, n_max = -1L)
  sig_df <- data.frame(c(dn,up), expression= c(rep(-1, length(dn)), rep(1, length(up))))
  return(sig_df)
}

#' Evaluation of gene expression data
#' @description Checks matrix characteristics and normalizes matrix with log cpm normalization
#'
#' @author Tania Pal \email{taniya.pal.094@cranfield.ac.uk}
#' @param expdata Gene exoression data matrix or dataframe
#'
#' @return Information about the matrix such as dimensions, range, mean and plots before and after normalization boxplots of matrix
#' @export
#'
#' @examples evaluatematrix(expdata)

evaluatematrix=function(expdata){
  #checking the dimensions
  numberofrows=nrow(expdata)
  numberofcolumns=ncol(expdata)
  dimensions=c(paste("number of rows:", numberofrows), (paste ("number of columns:", numberofcolumns)))

  #replacing the missing values(NA) with mean of the column
  for(i in 1:ncol(expdata)){
    expdata[is.na(expdata[,i]), i] = mean(expdata[,i], na.rm = TRUE)
  }

 #boxplot before normalization
  par(mfrow=c(2,1))
  plot.before.normalization=boxplot(expdata, col="blue", main="PLot before normalization")
  library(edgeR)
  normalized.expdata=cpm(expdata)
  log(normalized.expdata)
  for(i in 1:ncol(expdata)){
    normalized.expdata[,i][is.na(normalized.expdata[,i])]<-mean(normalized.expdata[,i],na.rm=TRUE)
  }
  plot.after.normalization=boxplot(normalized.expdata, col="blue", main="Plot after normalization")

  #statistical summary of the expression matrix
  library(matrixStats)
  matrix.mean=mean(expdata)
  min<-min(expdata)
  max<-max(expdata)
  #install.packages("matrixStats")
  matrix.range=c(paste("Minimum value:", min), paste("Maximum value:", max))
  matrix.stats=c(paste("Range of matrix::", matrix.range), paste("Mean of matrix::", matrix.mean))


  output= c(paste("Dimensions:" ,dimensions), paste("Statistics of matrix:",matrix.stats) )
  print(output)

}

