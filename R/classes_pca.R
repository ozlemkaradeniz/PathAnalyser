#' PCA plot visualising pathway-based classification of samples in a dataset
#' @description This function generates a PCA plot showing the clustering of
#' samples in the normalised expression data set, coloured by predicted pathway
#' activity labels: active (green), inactive (orange) and uncertain (blue). A
#' ideal pathway-based classification would generate a PCA plot showing tight
#' clustering of samples in each activity class and less overlap between classes.
#' @author Yi-Hsuan Lee \email{yi-hsuan.lee@cranfield.ac.uk}
#' @param normalized_data_filtered  A (logCPM) normalized gene expression data
#' matrix, with row names consisting of the HUGO gene symbols and column names
#' corresponding to the name / ID of each sample in the dataset
#' @param predicted_labels_df  a data frame containing the pathway
#' activity labels predicted by classification algorithm for each sample in the
#' dataset. The first column is called "sample" which contains the sample names
#' in the data set and the second column is called "class" containing the
#' corresponding predicted pathway activity of the sample.
#' @param pathway  pathway name, default is "Predicted Pathway Activity"
#' @importFrom plotly plot_ly layout %>%
#' @importFrom stats prcomp
#'
#' @export
#'
#' @examples
#' \dontrun{classes_pca(normalized_data_filtered, predicted_labels_df, 'ER')}
classes_pca <- function(normalized_data_filtered, predicted_labels_df,
           pathway = "Pathway Activity") {
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
      layout(title=pathway,
             xaxis=list(title="PC1"),
             yaxis=list(title="PC2"))
    pca
  }

