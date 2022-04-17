#' Reading gene expression data from file
#'
#' @description  Reads gene expression matrix data file from user,
#' removes NAs, if the first column consist of gene names,
#' gives the gene names to the row names of the matrix and deletes
#' the first column. Before converting the first column to row names
#' it checks for any duplicated gene present in the first column.
#' Following that the matrix is converted to a numeric data matrix.
#'
#'
#' @author Taniya Pal \email{taniya.pal.094@@cranfield.ac.uk}
#'
#' @param file  Path for gene expression matrix file, which is either tab or
#' comma value separated. The first column must contain gene names and the
#' subsequent column names should be the sample name or ID.
#'
#' @return Structured matrix containing gene symbols/IDs as row names and
#' sample IDs as column names and the Sample IDs as column names.
#' @importFrom reader get.delim
#' @importFrom utils read.delim
#' @importFrom stats na.omit
#'
#' @export
#'
#' @examples
#' \dontrun{read_expr_data("/Users/taniyapal/Documents/Group Project/TCGA_unannotated.txt")}
read_expression_data <- function(file){
  if (!file.exists(file)) {
    stop(sprintf('File "%s" not found.', file))
  }
  # determine the delimiter for the file whether it is "\t" or "," or " "
  delimiter <- tryCatch(get.delim(file),
                        warning=function(w){
                          stop("Gene expression data file is not delimited.")
                        })
  # reading the file provided by the user
  dataset <- read.delim(file, sep=delimiter, row.names=NULL, check.names = F)
  # removing genes with NAs as GSVA doesn't consider those genes
  dataset <- na.omit(dataset)
  # If duplicate genes present take mean of duplicate genes
  dataset <- duplicate_genes(dataset)

  # if first column name is not a sample name read that column as row names
  if (colnames(dataset)[1] == "") {
    row.names(dataset) <- dataset[, 1]
    dataset <- dataset[,-1]
  }
  # check for duplicated samples
  dataset <- duplicate_samples(dataset)
  # set gene names as first column
  row.names(dataset) <- dataset[, 1]
  dataset <- dataset[,-1]
  # converting the data frame to numeric matrix
  dataset <- data.matrix(dataset)
  # check for non-numeric characters in the gene
  data_matrix <- as.numeric(dataset)
  if (sum(is.na(data_matrix)) != 0) {
    stop("Some expression values in expression data set are non-numerical.")
  }
  return(data_matrix)
}


#' Reads up-regulated and down-regulated gene signatures from files
#'
#' @description Reads up and down regulated signature files from user
#' and structures it according to package requirements. It returns a
#' dataframe. The first column of the dataframe  lists together
#' both up and down regulated signatures. The second column
#' signifies whether they are up(+1) or down(-1) regulated.
#'
#' @author Taniya Pal \email{taniya.pal.094@cranfield.ac.uk}
#' @param  up_sig_file   Up-regulated gene-set format file
#' @param  down_sig_file  Down-regulated gene-set format file
#'
#' @return A data frame containing both up regulated and down
#' regulated signature files and signifying their up or down
#' expression with +1 and -1 respectively.
#' @importFrom utils read.delim
#' @export
#'
#' @examples
#' \dontrun{read_sign_data("ESR1_UP.v1._UP.csv","ESR1_DN.v1_DN.csv" )}
read_signature <- function(up_sig_file, down_sig_file){
  # reading the up regulated gene signature file
  up_sig <- read.delim(up_sig_file, comment.char="#", sep="\n")

  # reading the down regulated gene signature file
  dn_sig <- read.delim(down_sig_file, comment.char="#", sep="\n")

  # vector combining both up and down regulated signatures
  genes <- c(up_sig[,1], dn_sig[,1])
  # generate expression values of 1 for up-regulated and -1 for down-regulated
  expression <- c(rep(1, nrow(up_sig)), rep(-1, nrow(dn_sig)))

  # combining up and down regulated signatures in single data frame
  sig_df <- data.frame("gene"=genes, "expression"=expression)

  return(sig_df)
}
