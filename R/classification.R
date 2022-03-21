#' Classifies samples according to pathway activity using Gene Set Variation Analysis (GSVA)
#' @author Anisha Thind \email{a.thind@@cranfield.ac.uk}
#' @param sig_df Signature data frame
#' @param data_se Expression data set matrix
#' @param up_thresh Numerical vector of thresholds for the up-regulated gene set
#' @param dn_thresh Numerical vector of thresholds for the down-regulated gene set
#'
#' @return classes_df
#' @export
#'
#' @examples
#' # Default thresholds for up-regulated and down-regulated genesets
#' classify(HER2_sig_df, HER_data_se)
#' # up-regulated gene set threshold: bottom = -0.4, top = 0.4
#' # down-regulated gene set thresholds: bottom = -0.4, top = 0.4
#' classify(HER2_sig_df, HER_data_se, up_thresh=c(-0.4, 0.4), dn_thresh=c(-0.4, 0.4))
classify <- function(sig_df, data_se, up_thresh, dn_thresh) {
  require(GSVA)
  # check signature arg is data frame
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


  # check thresholds
  if (!missing(up_thresh) && !missing(dn_thresh)) {
    # ensure thresholds are numbers between
    stopifnot("Up-regulated gene set thresholds are not between -1 and 1."=all(up_thresh >= -1 && up_thresh <= 1),
              "Down-regulated gene set thresholds are not between -1 and 1."=all(dn_thresh >= -1 && dn_thresh <= 1))
  } else {
    # compute quantiles
    up_thresh <- quantile(sorted_up$up, probs=c(0.25, 0.75))
    dn_thresh <- quantile(sorted_dn$dn, probs=c(0.25, 0.75))
  }


  # HER 2 (pathway)
  classes_df <- data.frame(sample=rownames(sorted_up),
                          class=vector(mode="character", nrow(sorted_up)))
  row.names(classes_df) <- classes_df$samples
  # loop through genes if a gene has high score for up-regulated and has
  # low score for down-regulated genes, then that sample is positive for pathway
  # if reverse then negative, else uncertain
  classes_list <- lapply(rownames(sorted_up), function(i) {
   if (sorted_up[i, 1] > up_thresh[2] &&
       sorted_dn[i, 1] < dn_thresh[1]) {
      classes_df[i, 2] <- "Active"
    } else if (sorted_up[i, 1] < up_thresh[2] &&
               sorted_dn[i, 1] > dn_thresh[1]){
      classes_df[i, 2] <- "Inactive"
    } else {
     classes_df[i, 2] <- "Uncertain"
    }
  })

  # convert list to vector
  classes <- unlist(classes_list)
  classes <- factor(classes, levels=c("Active", "Inactive", "Uncertain"))

  # add vector to pathway df
  classes_df$class <- classes
  cat("Summary of sample classification based on pathway activity:\n")
  cat("--------------------------------------------------------------\n")
  cat(sprintf("Total number of samples: %d\n", nrow(classes_df)))
  print(summary(classes_df$class))
  return(classes_df)
}

#' Plots GSVA scores distribution for up-regulated and down-regulated sets of a gene signature
#'
#' @param sig_df gene signature given as a data frame
#' @param data_se gene expression data set provided as a matrix or data frame
#'
#' @return
#' @export
#'
#' @examples
visualiseGSVA <- function(sig_df, data_se) {
  require(GSVA)
  require(ggplot2)
  require(reshape2)

  # check signature arg is data frame
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
  colnames(scores_up) <- "Up"
  scores_dn <- as.data.frame(scores[,2])
  colnames(scores_dn) <- "Down"

  # sort GSVA scores in descending order for the up-regulated
  # and down-regulated set
  sorted_up <- scores_up[order(scores_up$Up, decreasing=T), ,drop=F]
  sorted_dn <- scores_dn[order(scores_dn$Down, decreasing=T), ,drop=F]
  sortedScores <- cbind(sorted_up, sorted_dn)
  meltScores <- melt(sortedScores,
                     id.vars=NULL,
                     measures.vars=c("Up", "Down"),
                     variable.name = "Geneset",
                     value.name="Score")

  # density plot
  ggplot(data=meltScores, aes(x=Score, fill=Geneset)) +
    geom_density() +
    xlab("GSVA Score") +
    ylab("Density") +
    scale_fill_discrete(name="Gene set") +
    theme_bw() +
    scale_y_continuous(expand=c(0, 0)) +
    scale_x_continuous(breaks=seq(-1, 1, 0.2)) +
    facet_grid(~ Geneset)
}


