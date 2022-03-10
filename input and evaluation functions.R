#read gene expression matrix from user

readmatrix=function(expdata, rowData, colData){
  #making sure the expression data is in the form of a matrix
  expdata=as.matrix(expdata)
  #converting the first column of matrix to rownames
  rownames(expdata)=expdata[,1]
  expdata=expdata[,-1]
  #creating Summarized Experiment type of Data
  input.es=SummarizedExperiment(assays=expdata, rowData=rowData, colData=colData)
  return(input_es)
  
  
}
readsig=function(up_sig, down_sig){
  #coverting the gene signature files as list
  upsig=as.list(up_sig)
  downsig=as.list(down_sig)
  return(up_sig)
  return(down_sig)
  
}

#evaluation of the gene expression matrix

evaluatematrix=function(expdata){
  #checking the dimensions
  numberofrows=nrow(expdata)
  numberofcolumns=ncol(expdata)
  dimensions=c(paste("number of rows:", numberofrows), (paste ("number of columns:", numberofcolumns)))
  return(dimensions)
  #replacing the missing values(NA) with mean of the column
  for(i in 1:ncol(expdata)){
    expdata[is.na(expdata[,i]), i] = mean(expdata[,i], na.rm = TRUE)
  }
  return(expdata)
  #plotting a boxplot of the gene expression matrix
  plot=boxplot(expdata[,1:ncol(expdata)])
  outliers=plot$out
  return(plot)
  return(outliers)
  
  #statistical summary of the expression matrix
  matrix.mean=colMeans(expdata)
  install.packages("matrixStats")
  library(matrixStats)
  matrix.range=colRanges(expdata)
  matrix.summary=data.frame(matrix.mean, matrix.range)
  return(matrix.summary)
}