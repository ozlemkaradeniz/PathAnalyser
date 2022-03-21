
#' @param readmatrix Reads gene expression matrix from user and structures it according to paage requirements
#' @return Structured matrix for the convenience of  the package
#' @examples readmatrix()
#' @export
readmatrix=function(matrix){
  rownames(expdata)=expdata[,1]
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
#' @param readsig Reads up and down regulated signature files from user and structures it according to package requirements
#' @return Up and down regulated signature files in the form a list for the convenience of  the package
#' @examples readsig()
#' @export

readsig=function(up_sig, down_sig){
  upsig=as.list(up_sig)
  down_sig=as.list(down_sig)
}
#' @param evaluatematrix Checks matrix characteristics and normalizes matrix with log cpm normalization
#' @return Information about the matrix such as dimensions, range, mean and plots before and after normalization boxplots of matrix
#' @examples evaluatematrix()
#' @export

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
  #install.packages("matrixStats")
  matrix.range=c(paste("Minimum value:", min), paste("Maximum value:", max))
  matrix.stats=c(paste("Range of matrix::", matrix.range), paste("Mean of matrix::", matrix.mean))
  
  
  output= c(paste("Dimensions:" ,dimensions), paste("Statistics of matrix:",matrix.stats) )
  print(output)
  
}

