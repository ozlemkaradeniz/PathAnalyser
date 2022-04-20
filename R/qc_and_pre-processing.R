#' Log CPM transformation of RNA-seq raw count data by using log CPM
#' @description  Performs a logCPM transformation of RNA-seq raw count data
#' using the counts per million (CPM) method from the edgeR library. In addition
#' to the log CPM transformation, the function can also plot a boxplot as sanity
#' check for logCPM transformation of the gene expression matrix. The first
#' series of boxplots display the distribution of the raw counts for each sample
#' in thousands, while the second series of boxplots show the distribution of
#' the logCPM normalised gene expression matrix.
#'
#' @author Rishabh Kaushik and Taniya Pal
#' \email{rishabh.kaushik.126@@cranfield.ac.uk, taniya.pal.094@@cranfiled.ac.uk}
#' @param dataset  An unnormalised gene expression matrix containing raw RNA-seq
#' (integer) counts with gene symbols as row names and sample IDs as column
#' names.
#' @param boxplot Optional argument that is a boolean (TRUE or FALSE) indicating
#' whether boxplots displaying before and after log CPM transformation should be
#' displayed. (Default=TRUE)
#'
#' @return A logCPM transformed gene expression data matrix with gene symbols as
#' row names and samples names / IDs as column names.
#' @importFrom edgeR cpm
#' @importFrom graphics boxplot
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{log_cpm_transform(data.matrix)}
log_cpm_transform <- function(dataset, boxplot = TRUE) {
  samples <- counts <- NULL
  if (!is.matrix(dataset)) {
    stop("Dataset provided is not a matrix.")
  }

  if (!is.numeric(dataset)) {
    stop("The gene expression matrix contains non-numerical values.")
  }

  if (any(dataset %% 1 != 0)) {
    stop("The unnormalised expression matrix does not contain integer counts.")
  }

  # normalise dataset using log CPM
  norm_data <- cpm(dataset, log = TRUE)

  # generate boxplots if boxplots remains selected
  if (boxplot == TRUE) {
    # box plots of samples before log CPM transformation
    scaled_data <- dataset / 1000
    scaled_data <- as.data.frame(scaled_data)
    tScaled <- t(scaled_data)
    melt_data <- melt(
      tScaled,
      id.vars = samples,
      varnames = c("samples", "genes"),
      value.name = "counts"
    )
    plot1 <- ggplot(melt_data, aes(x = samples, y = counts)) +
      geom_boxplot(fill = "slateblue") +
      ggtitle("Plot before log CPM transformation") +
      ylab("Raw Counts Per Thousand") +
      theme_bw() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )
    plot(plot1)

    norm_df <- as.data.frame(norm_data)
    tNorm_df <- t(norm_df)
    melt_norm <-
      melt(
        tNorm_df,
        id.vars = samples,
        varnames = c("samples", "genes"),
        value.name = "counts"
      )
    # box plots after transformation
    plot2 <- ggplot(melt_norm, aes(x = samples, y = counts)) +
      geom_boxplot(fill = "slateblue", alpha = 0.2) +
      ylab("Log Counts Per Million (CPM)") +
      ggtitle("Plot after log CPM transformation") +
      theme_bw() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )

    plot(plot2)
  }
  return(norm_data)
}

#' Validity check in gene signatures and gene expression datasets
#' @description This function performs a validity check against the gene
#' signature and gene expression data set and filters genes that are absent in
#' expression data set and /or are not expressed in at least 10% of the total
#' number of samples. A bar plot of mean-normalised counts for each gene is also
#' displayed.
#'
#' @author Yi-Hsuan Lee \email{yi-hsuan.lee@cranfield.ac.uk}
#' @param norm_data  Normalized gene expression data matrix
#' @param sig_df A signature data frame containing two columns: the first column
#' contains a list of gene (represented by gene symbols) that are differentially
#' expressed in a given pathway and their relative expression values are given
#' in the second column, with 1 representing up-regulated genes and -1
#' representing down-regulated genes.
#'
#' @return A filtered gene expression data matrix with gene symbols as row names
#' and sample names / IDs as column names
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{check_signature_vs_dataset(norm_data, sig_df)}
check_signature_vs_dataset <- function(norm_data, sig_df, barplot=TRUE) {
  # bind gene and counts variables locally to the function
  gene <- counts <- NULL

  # check expression data matrix
  if (!is.matrix(norm_data)) {
    stop("Gene expression matrix (norm_data) is not a matrix.")
  }

  if (!is.numeric(norm_data)) {
    stop("Gene expression matrix (norm_data) contains non-numerical data.")
  }

  # check signature data frame is in an appropriate format
  check_sig_df(sig_df)

  # filter gene signature in expression matrix
  filtered <- norm_data[rownames(norm_data) %in% sig_df[, 1],]
  delete_gene_list <- list()

  # calculate each gene present or absent in each sample/case
  # store which row append in delete_gene list
  if (nrow(filtered) != 0) {
    for (i in 1:nrow(filtered)) {
      absent <- 0
      present <- 0
      for (j in 1:ncol(filtered)) {
        if (filtered[i, j] >= 0) {
          present <- present + 1
        }
        else{
          absent <- absent + 1
        }
      }
      present_proportion <- present / (present + absent)
      if (present_proportion <= 0.1) {
        delete_gene_list <- c(delete_gene_list, rownames(filtered)[i])
      }
    }
    # delete gene if present proportion less than 0.1
    if (length(delete_gene_list) != 0) {
      for (i in 1:length(delete_gene_list)) {
        filtered <-
          filtered[!rownames(filtered) %in% delete_gene_list[[i]][1],]
      }
    }
    filtered_mat <- data.matrix(filtered)
    gene_count <- nrow(filtered_mat)
    cat("Number of feature present in expression dataset:",
        gene_count,
        "\n")
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

      if (barplot) {
        df_plot <-
          data.frame(legend = rep(c("average"), each = gene_count),
                     gene = rep(1:gene_count, 1))
        df_plot$counts <- ""
        for (i in 1:gene_count) {
          df_plot[i, 3] <- count_list$average[i]
        }
        df_plot$counts <- as.numeric(as.vector(df_plot$counts))

        p <- ggplot(data = df_plot, aes(x = gene, y = counts)) +
          geom_bar(stat = "identity", fill="#0072B2") +
          ggtitle("Mean normalised counts per gene") +
          ylab("Count") + theme_bw() +
          scale_x_continuous(expand=c(0, 0)) +
          theme(plot.title = element_text(hjust = 0.5),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank())
        print(p)
      }

      return(filtered_mat)
    }

  }
  else{
    cat("No genes in the gene expression martix are present in the gene signature.")
  }
}
