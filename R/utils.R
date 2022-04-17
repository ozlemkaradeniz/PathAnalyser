#' Check for duplicated genes
#'
#' @description Checks for duplicated gene names and finds the mean expression
#' values or counts for the duplicated gene for each sample.
#' @param dataset A data frame containing the gene expression data with the first
#' column containing gene names and subsequent columns named after each sample
#' in the data set.
#'
#' @return A data frame containing the expression data without duplicate gene
#' names with the first column containing gene names and subsequent columns
#' named after each sample in the data set.
#' @keywords internal
#' @noRd
#' @examples \dontrun{dedup_genes <- duplicate_genes(dataset)}
duplicate_genes <- function(dataset){
  genes <- dataset[,1]
  # extract duplicated genes from data set
  dup_genes <- dataset[duplicated(genes), 1]
  if (length(dup_genes) != 0) {
    warning(
      sprintf('Duplicate gene(s) "%s" found, ave. expression values for each sample for the duplicated gene(s) will be used.',
              paste(dup_genes, collapse=", ")))
    uniq_genes <- dataset[!duplicated(genes), 1]
    uniq_genes_df <- dataset[!(genes %in% dup_genes),]
    dup_genes_df <- dataset[genes %in% dup_genes,]
    # for each duplicate gene find mean expression values / counts for each
    # sample
    for (gene in dup_genes) {
      dup_gene_df <- dataset[genes== gene,]
      dedup_gene <- colMeans(dup_gene_df[2:ncol(dup_gene_df)])
      uniq_genes_df[nrow(uniq_genes_df) + 1, 2:ncol(dataset)] <- dedup_gene
      uniq_genes_df[nrow(uniq_genes_df), 1] <- gene
    }
    return(uniq_genes_df)
  }
  # if there are no duplicates return original data frame
  return(dataset)

}


#' Remove duplicate samples from an gene expression dataset
#'
#' @param dataset A data frame containing gene expression data with the first
#' column containing gene names and subsequent columns are named after each
#' sample in the data set.
#'
#' @return A data frame containing the data set without duplicate samples, with
#' the first column containing the gene names and subsequent columns named after
#' each sample in the data set.
#' @noRd
#' @examples \dontrun{dataset <- duplicate_samples(dataset)}
duplicate_samples <- function(dataset) {
  # removing duplicated samples from gene expression data frame
  if (any(duplicated(colnames(dataset)))) {
    dup_ind <- which(duplicated(colnames(dataset)), arr.ind = T)
    dup_samples <- unique(colnames(dataset)[dup_ind])
    warning(
      sprintf('Duplicate sample(s) "%s" found, removing duplicate sample(s) from data set.',
                    paste(dup_samples, collapse=", ")))
    dataset <- dataset[,unique(colnames(dataset))]
  }
  return(dataset)
}
