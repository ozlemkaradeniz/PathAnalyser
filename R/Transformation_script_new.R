
#' Transformation of the raw count data by using  log CPM (counts per million)
#' @description  Transformation of the raw count data by using  log CPM (counts per million) by calling
#' cpm method in edgeR library
#'
#' @author Rishabh Kaushik and Taniya Pal \email{rishabh.kaushik.126@cranfield.ac.uk, taniya.pal.094@cranfiled.ac.uk}
#' @param data_es  Gene expression matrix with gene IDs/symbols as row names
#' and Sample IDs as column names
#'
#' @return Transformed gene expression data matrix
#' @export
#'
#' @examples log_cpm(data_es)
#'

log_cpm <-function (data_es){
  require(edgeR)

  #box plot before transformation
  plot.before.transformation=boxplot(log(data_es+0.5), main="Plot before log cpm transformation", axes=F)

  data_es_log_cpm<-cpm(data_es, log=TRUE)

  #box plot after transformation
  plot.after.transformation=boxplot(log(data_es_log_cpm+0.5), main="Plot after log cpm transformation", axes=F)


  return(data_es_log_cpm)

}


