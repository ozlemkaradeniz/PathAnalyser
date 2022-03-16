#' Classifies samples according to pathway activity using GSVA
#' @author Anisha Thind \email{a.thind@@cranfield.ac.uk}
#' @param sig_df Signature dataframe
#' @param data_se Expression dataset matrix
#' @param thresh Threshold between 0-1 for pathway classification (default: 0.25)
#'
#' @return classes_df
#' @export
#'
#' @examples
#' # Default threshold
#' classify(HER2_sig_df, HER_data_se)
#' # Threshold = 0.05 (top 5% of up-regulated and bottom 5% of down-regulated genes)
#' classify(HER2_sig_df, HER_data_se, thresh=0.05)
classify <- function(sig_df, data_se, thresh=0.25) {
  require(GSVA)
  # check signature arg is dataframe
  if (!is.data.frame(sig_df)) {
    stop("Signature argument is not a dataframe.")
  }

  # check input is input matrix
  if (is.data.frame(data_se)) {
    data_se <- as.matrix(data_se)
  } else if (!is.matrix(data_se)) {
    stop("Expression dataset is not a matrix or dataframe object.")
  }

  # create list of gene sets (up-regulated and then down regulated) for gsva
  sigs <- list()
  sigs$up <- sig_df[sig_df$expression == 1, 1]
  sigs$dn <- sig_df[sig_df$expression == -1, 1]

  # run GSVA on expression data using 2 gene sets (up-regulated and
  # down-regulated) of the signature
  if (typeof(data_se) == "integer") {
    # If the data is of type integer, these are raw counts and follow Poisson
    scores <- gsva(data_se, sigs, kcdf="Poisson", verbose=F)
  } else {
    scores <- gsva(data_se, sigs, verbose=F)
  }

  # transpose the matrix
  scores <- t(scores)
  # create an up and down regulated matrix
  scores_up <- as.data.frame(scores[,1], row.names=rownames(scores))
  colnames(scores_up) <- "up"

  scores_dn <- as.data.frame(scores[,2])
  colnames(scores_dn) <- "dn"

  # sort GSVA scores in descending order for the up-regulated
  # and down-regulated set
  sorted_up <- scores_up[order(scores_up$up, decreasing=T), ,drop=F]
  sorted_dn <- scores_dn[order(scores_dn$dn, decreasing=T), ,drop=F]

  # compute quantiles according to user-specified threshold
  if (thresh > 0 && thresh < 1) {
    thresh_up <- quantile(sorted_up$up, probs=seq(0, 1, thresh))
    thresh_dn <- quantile(sorted_dn$dn, probs=seq(0, 1, thresh))
  } else {
    stop("Threshold given is not a decimal between 0 and 1.")
  }


  # HER 2 (pathway)
  classes_df <- data.frame(sample=rownames(sorted_up),
                          class=vector(mode="character", nrow(sorted_up)))
  row.names(classes_df) <- classes_df$samples
  # loop through genes if a gene has high score for up-regulated and has
  # low score for down-regulated genes, then that sample is positive for pathway
  # if reverse then negative, else uncertain
  classes_list <- lapply(rownames(sorted_up), function(i) {
   if (sorted_up[i, 1] > thresh_up[length(thresh_up) - 1] &&
       sorted_dn[i, 1] < thresh_dn[2]) {
      classes_df[i, 2] <- "Positive"
    } else if (sorted_up[i, 1] < thresh_up[2] &&
               sorted_dn[i, 1] > thresh_dn[length(thresh_dn) - 1]){
      classes_df[i, 2] <- "Negative"
    } else {
     classes_df[i, 2] <- "Uncertain"
    }
  })

  # convert list to vector
  classes <- unlist(classes_list)
  classes <- factor(classes, levels=c("Positive", "Negative", "Uncertain"))

  # add vector to pathway df
  classes_df$class <- classes
  cat("Summary of sample classification based on pathway activity:\n")
  cat("--------------------------------------------------------------\n")
  cat(sprintf("Total number of samples: %d\n", nrow(classes_df)))
  print(summary(classes_df$class))
  return(classes_df)
}
