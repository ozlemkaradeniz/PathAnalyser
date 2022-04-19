#' PCA plot visualising pathway-based classification of samples in a dataset
#' @description This function generates a PCA plot showing the clustering of
#' samples in the normalised expression data set, coloured by predicted pathway
#' activity labels: active (green), inactive (orange) and uncertain (blue). An
#' ideal pathway-based classification would generate a PCA plot showing tight
#' clustering of samples in each activity class and less overlap between classes.
#' @author Yi-Hsuan Lee \email{yi-hsuan.lee@cranfield.ac.uk}
#' @param norm_data  A (logCPM) normalized gene expression data
#' matrix, with row names consisting of the HUGO gene symbols and column names
#' corresponding to the name / ID of each sample in the dataset
#' @param predicted_labels_df A data frame containing the pathway
#' activity labels predicted by classification algorithm for each sample in the
#' dataset. The first column is called "sample" which contains the sample names
#' in the data set and the second column is called "class" containing the
#' corresponding predicted pathway activity of the sample.
#' @param pathway  The pathway name used in the title of the plot, default is
#' a placeholder "Pathway Activity".
#' @importFrom plotly plot_ly layout %>%
#' @importFrom stats prcomp
#'
#' @export
#'
#' @examples
#' \dontrun{classes_pca(norm_data, predicted_labels_df, 'ER')}
classes_pca <- function(norm_data, predicted_labels_df,
           pathway = "Pathway Activity") {
    # check normalised gene expression matrix (norm_data)
    if (!is.matrix(norm_data)) {
      stop("Normalised expression data set (norm_data) is not a matrix.")
    }

    if (!is.numeric(norm_data)){
      stop("Normalised expression data set (norm_data) contains non-numerical data.")
    }

    if (!is.data.frame(predicted_labels_df)) {
      stop("predicted_labels_df is not a data frame.")
    }

    data_filtered <- t(norm_data)
    row.names(predicted_labels_df) <- tryCatch(predicted_labels_df$sample,
                                               error=function(e){
                                                 stop('Predicted labels data frame does not have "sample" column.')
    })
    merged <- merge(predicted_labels_df, data_filtered, by='row.names', all=T)
    merged <- merged[,-1]
    prin_comp <- tryCatch(prcomp(merged[3:ncol(merged)], rank. = 2),
                          error=function(e){
                            stop("Sample names/IDs in the normalised expression matrix are incosistent with those in the predicted labels data frame.")
                          })
    components <- prin_comp[["x"]]
    components <- data.frame(components)
    components <- cbind(components, merged$class)
    row.names(components)<-merged$sample
    pca <- plot_ly(components,type = "scatter",x = ~PC1, y = ~PC2,text=rownames(components),
            mode="markers",color = components$`merged$class`,marker=list(size=11))%>%
      layout(title=pathway,
             xaxis=list(title="PC1"),
             yaxis=list(title="PC2"))
    pca
  }

