#' Reading gene expression data from file
#'
#' @description  Reads gene expression matrix data file which is either
#' tab/comma/white-space value separated text file. After reading in the input
#' file, any gene rows that contain NAs are removed and the data set is screened
#' for duplicates. Duplicate genes are reduced to single gene row entries, with
#' each value corresponding to the mean expression value for each sample for the
#' duplicated gene entry. If duplicate samples are detected then the only the
#' first sample of the duplicates are retained. Finally, the data frame is
#' converted to a numerical matrix, where row names represent gene names and
#' columns represent sample names or IDs.
#'
#'
#' @author Taniya Pal \email{taniya.pal.094@@cranfield.ac.uk}
#'
#' @param file  Path for gene expression matrix file, which is either tab /
#' comma / white-space value separated. The first column must contain gene names
#' and the subsequent column names should be the sample name or ID.
#'
#' @return A numerical matrix containing gene symbols/IDs as row names and
#' sample IDs as column names.
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
  dataset <- read.delim(file, sep=delimiter, row.names=NULL, check.names = F,
                        header=T)

  # removing genes with NAs as GSVA doesn't consider those genes
  dataset <- na.omit(dataset)
  # # If duplicate genes present take mean of duplicate genes
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
  data_matrix <- as.matrix(dataset)
  # check for non-numeric characters in the gene expression matrix
  tryCatch(sum(is.na(matrix(as.numeric(data_matrix), nrow=nrow(dataset)))) != 0,
           warning=function(w){
     stop("Some expression values in expression data set are non-numerical.")
  })
  rownames(data_matrix) <- rownames(dataset)
  colnames(data_matrix) <- colnames(dataset)
  return(data_matrix)
}


#' Reads up-regulated and down-regulated gene signatures from gene set files
#'
#' @description Reads up and down regulated signature files provided either in
#' gene set file format (.grp) or gene matrix transposed format (.gmt) creating
#' a data frame which the first column named "gene" containing the gene names
#' and the second column called "expression" containing the corresponding
#' expression value for a gene in the gene signature, where 1 signifies
#' up-regulation and -1 represents down-regulation of the gene in the gene
#' signature.
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
  # read files
  up_sig <- read_sig_file(up_sig_file)
  dn_sig <- read_sig_file(down_sig_file)
  # vector combining both up and down regulated signatures
  genes <- c(up_sig[,1], dn_sig[,1])
  # generate expression values of 1 for up-regulated and -1 for down-regulated
  expression <- c(rep(1, nrow(up_sig)), rep(-1, nrow(dn_sig)))

  # combining up and down regulated signatures in single data frame
  sig_df <- data.frame("gene"=genes, "expression"=expression)

  return(sig_df)
}

#' Read gene set file constituting one of the two gene-set of the gene signature
#'
#' @description Reads data from gene set format files (.grp) or gene matrix
#' transposed format files (.gmt).
#' @param sig_file A geneset file containing either the up-regulated or
#' down-regulated gene set for the gene signature to be used in classification.
#' The file can be in gene set file format, gene matrix transposed format or a
#' tab or comma value separated file.
#' @importFrom utils read.delim
#' @return a data frame containing gene signature information
#'
#' @noRd
#' @examples \dontrun{up_sig <- read_sig_file("geneset.grp")}
read_sig_file <- function(sig_file){
  sig <- NULL
  # if file is in gene set format (.grp)
  if (grepl("\\.grp$", sig_file)) {
    # reading the up regulated gene signature file
    sig <- tryCatch(read.delim(sig_file, comment.char="#", sep="\n"),
                       warning=function(w){
                         stop(sprintf('File "%s" is an incorrectly formated gene set format file (.grp)',
                              sig_file))
                       })
  # if file is in gene matrix transposed format (.gmt)
  } else if (grepl("\\.gmt$", sig_file)){
    lines <- readLines(sig_file, warn=FALSE)
    geneset <- unlist(strsplit(lines, "\t"))
    # the data frame should contain all NAs when converting to numeric, if it
    # only contains characters / words, if there are numbers i.e. some non-NAs
    # then there are number present which means the file contains numerical data
    if (suppressWarnings(!all(is.na(as.numeric(geneset))))) {
      stop(sprintf('File "%s" contains numerical data and is incorrectly formatted for gene matrix transposed format file (.gmt)',
           sig_file))
    }
    geneset_name <- geneset[1]
    sig <- as.data.frame(geneset[-c(1,2)])
    colnames(sig) <- geneset_name
    # or a delimited text file
  } else {
    stop(sprintf('Unable to read signature data from file "%s"', sig_file))
  }
  return(sig)
}
