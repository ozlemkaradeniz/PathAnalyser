
#' Validity check in gene signatures and gene expression datasets
#' @description performs validity check against gene signatures and gene expression dataset and
#'          filters gene expression dataset accordingly
#'
#' @author Yi-Hsuan Lee \email{yi-hsuan.lee@cranfield.ac.uk}
#' @param data_norm  Normalized Gene expression data matrix
#' @param sig_df  signature gene-set
#'
#'
#'
#' @return Filtered Normalized Gene expression data matrix
#' @import ggplot2
#' @export
#'
#' @examples \dontrun{check_signature_vs_dataset(data_norm, sig_df)}
check_signature_vs_dataset <-function(data_norm, sig_df) {
    sig_df<-as.data.frame(sig_df)
    sig_df<-as.data.frame(sig_df[,1])
    # filter gene signature in expression matrix
    filtered <- data_norm[rownames(data_norm) %in% sig_df[, 1],]
    delet_gene_list <- list()

    # calculate each gene present or absent in each case
    # store which row append in delet_gene list
    if (nrow(filtered) != 0) {
      for (i in 1:nrow(filtered)) {
        absent <- 0
        present <- 0
        for (j in 1:ncol(filtered)) {
          if (filtered[i, j] != 0) {
            present <- present + 1
          }
          else{
            absent <- absent + 1
          }
        }
        present_proportion <- present / (present + absent)
        if (present_proportion <= 0.1) {
          delet_gene_list <- c(delet_gene_list, rownames(filtered)[i])
        }
      }
      # delet present proportion less than 0.1
      if (length(delet_gene_list)!=0) {
        for (i in 1:length(delet_gene_list)) {
          filtered <-
            filtered[!rownames(filtered) %in% delet_gene_list[[i]][1],]
        }
      }
      filtered_mat <- data.matrix(filtered)
      gene_count <- nrow(filtered_mat)
      cat("Number of feature present in expression dataset:",
          gene_count,"\n")
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

        }




         df_plot<- data.frame( legend= rep(c("average"), each = gene_count),
                      gene = rep(1:gene_count, 1))
        df_plot$counts<-""
        for(i in 1:gene_count){
          df_plot[i,3]<-count_list$average[i]
        }
        df_plot$counts <- as.numeric(as.vector(df_plot$counts))

        p<-ggplot(data=df_plot, aes(x=gene, y=counts)) +
          geom_bar(stat="identity")+labs(title = "Mean normalised counts per gene")
      print(p)
        return(filtered_mat)
      }

    }
    else{
      cat("No gene present in sigunature")
    }
  }

