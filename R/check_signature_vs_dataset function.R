

check_HR2_signature_vs_dataset <-
  function(data_es,sig_up_df, sig_dn_df) {
    # filter gene signature in expression matrix
    up_data <- data_es[rownames(data_es) %in% sig_up_df[, 1],]
    dn_data <- data_es[rownames(data_es) %in% sig_dn_df[, 1],]
    filtered <- rbind(up_data, dn_data)
    filtered_mat <- data.matrix(filtered)
    gene_count <- nrow(filtered_mat)
    cat("Number of feature present in expression dataset:", gene_count)
    if (gene_count != 0) {
      count_list <- list(0)
      count_list$minimun <- 0
      count_list$maximun <- 0
      count_list$average <- 0
      for (i in 1:nrow(filtered_mat)) {
        gene_case_count <- list() # get gene count from each case
        
        gene_case_count <- filtered_mat[i,]
        count_list$minimun[i] <- min(gene_case_count)
        count_list$maximun[i] <- max(gene_case_count)
        count_list$average[i] <- mean(gene_case_count)
        cat(
          "\n",
          row.names(filtered_mat)[i],
          ":\nMinimun count:",
          count_list$minimun[i],
          "\nMaximun count:",
          count_list$maximun[i],
          "\nAverage count:",
          count_list$average[i],
          "\n"
        )
        
      }
      return(filtered_mat)
    }
  }

