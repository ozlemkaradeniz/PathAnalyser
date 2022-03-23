#' Transforms gene expression dataset according to CPM read count normalization
#' @description Transforms gene expression dataset according to CPM read count normalization,
#'              calls cmp method in edgeR library
#'
#' @author Rishabh Kaushik \email{rishabh.kaushik.126@cranfield.ac.uk}
#' @param data_es  Gene expression data matrix
#'
#'
#' @return Normalized Gene expression data matrix
#' @export
#'
#' @examples transform(data_es)
#'

transform <-function (data_es){
   require(edgeR)

    data_es_log_cpm<-cpm(data_es, log=TRUE)
    return(data_es_log_cpm)

}


