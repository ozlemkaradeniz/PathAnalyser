#' Reading gene expression data from file
#'
#' @description  Reads gene expression data file from user,
#' removes NAs, if the first column consist of gene names,
#' gives the gene names to the row names of the matrix and deletes
#' the first column. Before converting the first column to row names
#' it checks for any duplicated gene present in the first column.
#' Following that the matrix is converted to a numeric data matrix.
#'
#'
#' @author Taniya Pal \email{taniya.pal.094@cranfield.ac.uk}
#'
#' @param file_name  Full path of the gene expression data file name
#'
#' @return Structured matrix containing gene symbols/IDs as row names and
#' sample IDs as column names and the Sample IDs as column names.
#' @importFrom reader get.delim
#' @importFrom utils read.delim
#' @importFrom stats na.omit
#'
#' @export
#'
#' @examples read_expression_data("/Users/taniyapal/Documents/Group Project/TCGA_unannotated.txt")

read_expression_data <- function(file_name){
  # getting the delimiter for the file whether it is "\t" or "," or " "
  delimiter <- get.delim(file_name)

  # reading the file provided by the user
  data_se <- read.delim(file_name, sep=delimiter, check.names=F)

  # removing the NAs
  data_se <- na.omit(data_se)
  # removing duplicated samples from gene expression data frame
  data_se <- data_se[,unique(colnames(data_se))]

  #converting the data frame in numeric matrix
  data_se <- data.matrix(data_se)
  return(data_se)
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
#' @param  up_sig   Up-regulated gene-set
#' @param  down_sig  Down-regulated gene-set
#'
#' @return A data frame containing both up regulated and down
#' regulated signature files and signifying their up or down
#' expression with +1 and -1 respectively.
#' @importFrom utils read.delim
#' @export
#'
#' @examples read_signature_data("ESR1_UP.v1._UP.csv","ESR1_DN.v1_DN.csv" )

read_signature_data <- function(up_sig_file, down_sig_file){
  # reading the up regulated gene signature file
  up_sig <- read.delim(up_sig_file, comment="#", sep="\n")

  # reading the down regulated gene signature file
  dn_sig <- read.delim(down_sig_file, comment="#", sep="\n")

  # vector combining both up and down regulated signatures
  genes <- c(up_sig[,1], dn_sig[,1])
  # generate expression values of 1 for up-regulated and -1 for down-regulated
  expression <- c(rep(1, nrow(up_sig)), rep(-1, nrow(dn_sig)))

  # combining up and down regulated signatures in single data frame
  sig_df <- data.frame("genes"=genes, "expression"=expression)

  return(sig_df)
}

#'Reading both expression data and signature files,also checking overlap
#'between matrix and signature gene symbols
#'
#' @description Reads both expression data and up and down regulated signature
#' data, if check_matrix_sig_overlap parameter in the function is true, returns
#' a Venn diagram representating which gene symbols
#' overlap between the signature file
#' and the expression matrix file.
#' The user gets an estimates of how many gene symbols overlap
#' between the two files.
#'
#' @author Taniya Pal \email{taniya.pal.094@cranfield.ac.uk}
#'
#' @param input_filename The expression matrix file name
#' @param up_sig_filename The up regulated signature data file name
#' @param dn_sig_filename The down regulated signature data file name
#' @param check_matrix_sig_overlap Boolean, if true returns a Venn diagram plotting
#' overlap between matrix and gene signature file symbols/IDs.
#'
#' @return A list with 3 objcets:
#' 1) formatted expression matrix with gene IDs/symbols as row names and
#' samples IDs as column names
#' 2) a data frame, whose first column has combined up and down regulated signature
#' IDs/symbols, second column has representation of up regulated expression as +1
#' and down regulated expression as -1
#' 3) a Venn diagram representing the overlap of gene symbols/IDs between the
#' two files
#' @import VennDiagram
#'
#' @export
#' @examples read_input_file("TCGA_unannotated.txt","ESR1_UP.v1_UP.csv",
#' "ESR1_DN.v1_DN.csv", check_mtarix_sig_overlap=T)

read_input_file=function(input_filename, up_sig_filename, dn_sig_filename, check_matrix_sig_overlap=F){

  #reading gene expression matrix
  matrix=read_expression_data(input_filename)

  #reading up and down regulated signature files
  sig_df=read_signature_data(up_sig_filename, dn_sig_filename)



  #if the user wants to display the Venn Diagram
  if (check_matrix_sig_overlap==T){

    #listing the gene symbols in expression matrix
    rownames=as.list(rownames(matrix))

    #listing the gene symbols in signature
    signatures=as.list(sig_df$Signatures)

    #finding the length of intersection between the two lists
    intersection=length(intersect(rownames, signatures))

    #Setting the color
    myCol=c("green", "yellow")


    plot=draw.pairwise.venn(
      area1 = length(rownames), area2=length(signatures),
      category = c("Expression Matrix" , "Signature"),
      main = 'Checking overlap between matrix and gene signature symbols', lwd = 2,
      lty = 'blank',
      fill = myCol, cat.cex = 0.7,
      cat.fontface = c("bold", "bold"),
      cat.default.pos = "text",
      cat.pos = c(-27, 27),
      cat.dist = c(-0.055, 0.055),
      cat.fontfamily = c("sans", "sans"), cross.area=94, scaled=F)


    output=list(matrix, sig_df, plot)
    return(output)

  }

}
