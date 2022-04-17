#' Log CPM transformation of RNA-seq raw count data by using log CPM
#' @description  Performs a logCPM transformation of RNA-seq raw count data
#' using the counts per million (CPM) method from edgeR library. In addition to
#' the log CPM transformation, the function plots two boxplots as sanity check
#' for logCPM transformation of a gene expression matrix. The first boxplot
#' displays the distribution of the raw counts for each sample in thousands,
#' while the second boxplot shows the distribution of the logCPM normalised gene
#' expression matrix.
#'
#' @author Rishabh Kaushik and Taniya Pal \email{rishabh.kaushik.126@cranfield.ac.uk, taniya.pal.094@cranfiled.ac.uk}
#' @param data_es  An unnormalised gene expression matrix containing raw RNA-seq
#' counts with gene IDs/symbols as row names and sample IDs as column names
#'
#' @return logCPM transformed gene expression data matrix
#' @importFrom edgeR cpm
#' @importFrom graphics boxplot
#' @export
#'
#' @examples
#' \dontrun{log_cpm_transformation(formatted_matrix)}
log_cpm_transform <- function(data_es){
  #box plot before transformation
  boxplot(data_es / 1000, main="Plot before log cpm transformation", xlab=" ",
          ylab="Raw counts per thousand", xaxt="n")

  data_es_log_cpm<-cpm(data_es, log=TRUE)

  #box plot after transformation
  boxplot(data_es_log_cpm, main="Plot after log cpm transformation", xlab=" ",
          ylab="Log Counts Per Million (CPM)", xaxt="n")

  return(data_es_log_cpm)

}


