
#' pca plot for normalized gene expression dataset
#' @description plots pca for normalized gene expression dataset with the predicted labels
#'
#' @author Yi-Hsuan Lee \email{yi-hsuan.lee@cranfield.ac.uk}
#' @param normalized_data_filtered  Normalized Gene expression data matrix
#' @param predicted_labels_df  dataframe which labels for each sample predicted by classification algorithm
#' @param pathway_name  pathway name ER or HER2, default is ER
#'
#'
#' @return
#' @export
#'
#' @examples classes_pca(normalized_data_filtered, predicted_labels_df, 'ER')
#'

classes_pca <-
  function(normalized_data_filtered,
           predicted_labels_df,
           pathway_name = "ER") {

    library(ggfortify)

    data_filtered<-t(normalized_data_filtered)
    row.names(predicted_labels_df)<-predicted_labels_df$sample
    Merged<-merge(predicted_labels_df,data_filtered,by='row.names',all=T)
    Merged<-Merged[,-1]
    pca_res <- prcomp(Merged[3:ncol(Merged)], scale. = TRUE)
    autoplot(pca_res, data = Merged, colour = "class",main=pathway_name)

  }
