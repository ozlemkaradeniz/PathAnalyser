
#' Transformation of the raw count data by using  log CPM (counts per million) along
#' with sanity check for transformation
#' @description  Transformation of the raw count data by using  log CPM (counts per million)
#' by calling cpm method in edgeR library. Along with performing transformation 
#' of raw data,the function plots the distribution of raw data before transformation 
#' in the form of a boxplot. Then, it normalizes the raw counts of gene
#' expression with the log cpm transformation method and returns a boxplot 
#' of the gene expression matrix after transformation for sanity check of the
#' transformation. The user can check the distribution of the gene expression values 
#' after log cpm transformation with the help of the box plot.
#'
#' @author Rishabh Kaushik and Taniya Pal \email{rishabh.kaushik.126@cranfield.ac.uk, taniya.pal.094@cranfiled.ac.uk}
#' @param data_es  Gene expression matrix with gene IDs/symbols as row names
#' and Sample IDs as column names (return value of read_expression matrix function)
#'
#' @return Transformed gene expression data matrix
#' @export
#'
#' @examples log_cpm_transformation(formatted_matrix)
#'

log_cpm_transformation <-function (data_es){
  require(edgeR)
  
  #box plot before transformation
  boxplot(log(data_es+0.5), main="Plot before log cpm transformation", axes=F)

  data_es_log_cpm<-cpm(data_es, log=TRUE)

  #box plot after transformation
  boxplot(data_es_log_cpm, main="Plot after log cpm transformation", axes=F)

  
  return(data_es_log_cpm)

}


