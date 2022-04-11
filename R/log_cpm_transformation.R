
#' Log CPM transformation of RNA-seq raw count data by using log CPM
#' @description  Performs a logCPM transformation of RNA-seq raw count data
#' using the cpm method from edgeR library. Along with performing transformation
#' of raw data,the function plots the distribution of raw data before the
#' transformation in the form of a box plot. Then, it normalizes the raw counts
#' of gene expression with the log cpm transformation method and returns a boxplot
#' of the gene expression matrix after transformation for sanity check of the
#' transformation. The user can check the distribution of the gene expression values
#' after log cpm transformation with the help of the box plot.
#'
#' @author Rishabh Kaushik and Taniya Pal \email{rishabh.kaushik.126@cranfield.ac.uk, taniya.pal.094@cranfiled.ac.uk}
#' @param data_es  Gene expression matrix with gene IDs/symbols as row names
#' and Sample IDs as column names (return value of read_expression matrix function)
#'
#' @return Transformed gene expression data matrix
#' @importFrom edgeR cpm
#' @importFrom graphics boxplot
#' @export
#'
#' @examples
#' \dontrun{log_cpm_transformation(formatted_matrix)}
log_cpm_transformation <-function (data_es){
  #box plot before transformation
  boxplot(log(data_es+0.5), main="Plot before log cpm transformation", xlab=" ",
          ylab="Log Raw Counts", xaxt="n")

  data_es_log_cpm<-cpm(data_es, log=TRUE)

  #box plot after transformation
  boxplot(data_es_log_cpm, main="Plot after log cpm transformation", xlab=" ",
          ylab="Log Counts Per Million (CPM)", xaxt="n")

  return(data_es_log_cpm)

}


