
#' pca plot for normalized gene expression dataset
#' @description plots pca for normalized gene expression dataset with the predicted labels
#'
#' @author Yi-Hsuan Lee \email{yi-hsuan.lee@cranfield.ac.uk}
#' @param normalized_data_filtered  Normalized Gene expression data matrix
#' @param predicted_labels_df  dataframe which labels for each sample predicted by classification algorithm
#' @param pathway_name  pathway name, default is "ER"
#' @importFrom plotly plot_ly layout %>%
#'
#' @return
#' @export
#'
#' @examples classes_pca(normalized_data_filtered, predicted_labels_df, 'ER')
#'
classes_pca <- function(normalized_data_filtered, predicted_labels_df,
           pathway_name = "ER") {

    data_filtered<-t(normalized_data_filtered)
    row.names(predicted_labels_df)<-predicted_labels_df$sample
    Merged<-merge(predicted_labels_df,data_filtered,by='row.names',all=T)
    Merged<-Merged[,-1]

    prin_comp <- prcomp(Merged[3:ncol(Merged)], rank. = 2)
    components <- prin_comp[["x"]]
    components <- data.frame(components)
    components <- cbind(components, Merged$class)
    row.names(components)<-Merged$sample
    pca<-plot_ly(components,type = "scatter",x = ~PC1, y = ~PC2,text=rownames(components),
            mode="markers",color = components$`Merged$class`,marker=list(size=11))%>%
      layout(title=pathway_name,
             xaxis=list(title="PC1"),
             yaxis=list(title="PC2"))
    pca
  }

