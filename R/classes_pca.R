

classes_pca <-
  function(normalized_data_filtered,
           predicted_labels_df,
           pathway_name = "ER") {

    library(ggfortify)

    data_filtered<-t(normalized_data_filtered)

    Merged<-merge(predicted_labels_df,data_filtered,by='row.names',all=T)
    Merged<-Merged[,-2]
    pca_res <- prcomp(Merged[3:ncol(Merged)], scale. = TRUE)
    autoplot(pca_res, data = Merged, colour = "classes",main=pathway_name)

  }


