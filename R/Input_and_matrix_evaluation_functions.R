#' Reading gene expression data matrix
#' @description  Reads gene expression matrix from user and structures it according to page requirements
#'
#' @author Tania Pal \email{taniya.pal.094@cranfield.ac.uk}
#' @param matrix  Gene expression data matrix
#'
#' @return Structured matrix for the convenience of  the package
#' @export
#'
#' @examples readmatrix(matrix)

readmatrix=function(expdata){
  #rownames(expdata)=expdata[,1]
  expdata=expdata[-1,-1]
  for (i in 1:ncol(expdata)){
    class(expdata[,i])="numeric"
  }

  for(i in 1:ncol(expdata)){
    expdata[,i][is.na(expdata[,i])]<-mean(expdata[,i],na.rm=TRUE)
  }
  expdata=as.matrix(expdata)
  return(expdata)

}

#' Reading up-regulated and down-regulated gene signatures
#' @description Reads up and down regulated signature files from user and structures it according to package requirements
#'
#' @author Tania Pal \email{taniya.pal.094@cranfield.ac.uk}
#' @param  up_sig   Up-regulated gene-set
#' @param  down_sig  Down-regulated gene-set
#'
#' @return Up and down regulated signature files in the form a list for the convenience of  the package
#' @export
#'
#' @examples readsig(up_sig, down_sig)

readsig=function(up_sig, down_sig){
  upsig=as.list(up_sig)
  down_sig=as.list(down_sig)
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

