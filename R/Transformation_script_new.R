
#' Normalisation of the raw count data by using  log CPM (counts per million)
#' @description  Normalisation of the raw count data by using  log CPM (counts per million) by calling
#'  cmp method in edgeR library
#'
#' @author Rishabh Kaushik \email{rishabh.kaushik.126@cranfield.ac.uk}
#' @param data_es  Gene expression data matrix from TCGA
#'
#'
#' @return Normalized Gene expression data matrix
#' @export
#'
#' @examples log_cpm(data_es)
#'

log_cpm <-function (data_es){
   require(edgeR)

    data_es_log_cpm<-cpm(data_es, log=TRUE)
    return(data_es_log_cpm)

}


