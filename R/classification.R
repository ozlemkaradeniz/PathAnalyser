#' Classification Using Absolute GSVA Score Thresholds
#' @description Classifies samples according to pathway activity using the GSVA
#' algorithm with absolute GSVA scores thresholds for assessing expression
#' consistency with the up-regulated  and down-regulated gene-sets of the gene
#' signature.
#' @details The absolute threshold GSVA function classifies samples by pathway
#' activity using
#' @author Anisha Thind \email{a.thind@@cranfield.ac.uk}
#' @param sig_df Gene expression signature for a specific pathway given as data
#' frame containing a list of genes that are the most differentially expressed
#' when a specific pathway is active. Up-regulated genes are given an
#' expression of 1, while down-regulated genes are given an expression of -1.
#' @param data_se Normalized expression data set matrix given as a matrix or
#' data frame comprising the expression levels of genes for each sample in a
#' data set. Gene expression matrices can be from RNASeq or microarray
#' transcriptomic data sets.
#' @param up_thresh.low Number denoting the absolute GSVA score threshold for the
#' categorizing a sample as inconsistent expression with the up-regulated
#' gene-set from the gene signature.
#' @param up_thresh.high Number denoting the absolute GSVA score threshold for
#' categorizing a sample as having consistent expression with the up-regulated
#' gene set from the gene signature.
#' @param dn_thresh.low Number denoting the absolute GSVA score threshold for
#' the categorizing a sample as consistent expression with the down-regulated
#' gene-set from the gene signature.
#' @param dn_thresh.high Number denoting the absolute GSVA score threshold for
#' the categorizing a sample as inconsistent expression with the down-regulated
#' gene-set from the gene signature.
#'
#' @return A data frame containing a list of samples and their classified
#' pathway activity (Active, Inactive or Uncertain).
#' @export
#'
#' @examples
#' # Default thresholds for up-regulated and down-regulated gene-sets
#' classes_df <- classify_GSVA_abs(sig_df, data_se, up_thresh.low=-0.25,
#'      up_thresh.high=0.25, dn_thresh.low=-0.25, dn_thresh.high=0.35)
classify_GSVA_abs <- function(sig_df, data_se, up_thresh.low,
                              up_thresh.high,
                              dn_thresh.low,
                              dn_thresh.high){
  # check thresholds
  if (!is.numeric(up_thresh.high) && !is.numeric(dn_thresh.high) &&
      !is.numeric(up_thresh.low) && !is.numeric(dn_thresh.low)
      ) {
    stop("All thresholds provided must be numbers.")
    # reverse order if the low threshold is larger than the high threshold
  } else if (up_thresh.low > up_thresh.high) {
    stop("The high expression threshold for up-regulated gene-set
    (up_thresh.high) must be higher than the low threshold for the up-regulated
    gene-set (up_thresh.low)")
  } else if (dn_thresh.low > dn_thresh.high) {
    stop("The high expression threshold for down-regulated gene-set
    (dn_thresh.high) must be higher than the low threshold for down-regulated
    gene-set (dn_thresh.low)")
  }
  scores <- .run_GSVA(sig_df, data_se)
  classes_df <- .classify(scores, up_thresh.low, up_thresh.high, dn_thresh.low,
                         dn_thresh.high)
  return(classes_df)
}

#' Sample classification according to pathway activity using a percentile
#' threshold for assessing expression consistency with a gene signature.
#' @description Classify samples by ranking samples by abundance of expression
#' calculated via the GSVA algorithm. Samples are assessed for expression
#' consistency with both the up-regulated and down-regulated gene-sets using
#' percentile thresholds during the pathway activity sample classification.
#'
#' @author Anisha Thind \email{a.thind@@cranfield.ac.uk}
#' @param sig_df Gene expression signature for a specific pathway given as data
#' frame containing a list of genes that are the most differentially expressed
#' when a specific pathway is active. Up-regulated genes are given an
#' expression of 1, while down-regulated genes are given an expression of -1.
#' @param data_se Normalized gene expression matrix or data frame containing
#' expression values for samples in a data set (RNA Seq read counts or
#' microarray data).
#' @param up_thresh Percentile threshold of samples for checking consistency of
#' gene expression of a sample with first the up-regulated and then
#' down-regulated gene-set of the gene signature (default= 25% (quartile)).
#'
#' @return data frame containing predicted pathway activity classes for each
#' sample in the user provided data set (Active, Inactive or Uncertain).
#' @export
#'
#' @examples
#' # default using quartile threshold (25th percentile)
#' classes_df <- classify_GSVA_percent(ER_sig_df, ER_data_se1)
#' # custom percentile threshold e.g. 30th percentile
#' classes_df <- classify_GSVA_percent(ER_sig_df, ER_data_se1,
#'        thresh_percent=30)
classify_GSVA_percent <- function(sig_df, data_se, thresh_percent=25){
  # check threshold is a number between 0-100%
  thresh <- thresh_percent / 100
  scores <- .run_GSVA(sig_df, data_se)
  # compute percentiles
  tryCatch(up_thresh <- quantile(scores$Up[,1], c(thresh, 1-thresh)),
           error=function(e){
    stop("Percentile threshold given is not a number between 0 and 100.")
  })
  dn_thresh <- quantile(scores$Down[,1], c(thresh, 1-thresh))
  classes_df <- .classify(scores, up_thresh.low=up_thresh[1], up_thresh.high=up_thresh[2],
           dn_thresh.low=dn_thresh[1], dn_thresh.high=dn_thresh[2])
  return(classes_df)
}

#' GSVA score density distribution plot
#' @description Plots GSVA scores distribution for up-regulated and
#' down-regulated gene-sets of a provided gene expression signature
#' @author Anisha Thind \email{a.thind@@cranfield.ac.uk}
#' @param sig_df Normalized gene expression signature for a specific pathway
#' given as data frame containing a list of genes that are the most
#' differentially expressed when a specific pathway is active. Up-regulated
#' genes are given an expression of 1, while down-regulated genes are given an
#' expression of -1.
#' @param data_se Expression data set matrix given as a matrix or data frame
#' comprising the expression levels of genes in a set of samples. Gene
#' expression matrices can be from RNASeq
#' or microarray transcriptomic data sets.
#'
#' @return density plot displaying distribution of GSVA scores obtained for
#' the samples using up-regulated and down-regulated gene-sets from the gene
#' signature
#' @export
#'
#' @examples
#' gsva_scores_dist(ER_sig_df, ER_data_se)
gsva_scores_dist <- function(sig_df, data_se) {
  require(ggplot2)
  require(reshape2)

  # run GSVA using data provided
  gsva_scores <- .run_GSVA(sig_df, data_se)

  sortedScores <- cbind(gsva_scores$Up, gsva_scores$Down)

  # reshape data for plotting
  meltScores <- melt(sortedScores,
                     id.vars=NULL,
                     measures.vars=c("Up", "Down"),
                     variable.name = "Geneset",
                     value.name="Score")

  levels(meltScores$Geneset) <- c("Up-regulated Gene Signature",
                                  "Down-regulated Gene Signature")

  # density plot
  plot <- ggplot(data=meltScores, aes(x=Score)) +
    geom_density(fill="#69b3a2", alpha=0.8) +
    xlab("GSVA Score") +
    ylab("Density") +
    theme_bw() +
    scale_y_continuous(expand=c(0, 0)) +
    scale_x_continuous(breaks=seq(-1, 1, 0.2)) +
    facet_grid(~ levels(Geneset))
  return(plot)
}


#' Performs GSVA algorithm using both up-regulated and down-regulated gene-sets
#' of the gene signature
#' @description Generates GSVA scores for each sample generated by running GSVA
#' algorithm on the two gene-sets of the gene signature provided.
#' @author Anisha Thind \email{a.thind@@cranfield.ac.uk}
#' @param sig_df Gene expression signature for a specific pathway given as data
#' frame containing a list of genes that are the most differentially expressed
#' when a specific pathway is active. Up-regulated genes are given an
#' expression of 1, while down-regulated genes are given an expression of -1.
#' @param data_se Normalized expression data set matrix given as a matrix or
#' data frame comprising the expression levels of genes for each sample in a
#' data set. Gene expression matrices can be from RNASeq or microarray
#' transcriptomic data sets.
#' @keywords internal
#'
#' @return GSVA scores in the form of a list, containing two data frames: GSVA
#' scores for up-regulated gene-set for each sample, and GSVA scores for the
#' down-regulated gene-set GSVA for each sample
#'
#' @examples
#' run_GSVA(ER_sig_df, ER_data_se1)
.run_GSVA <- function(sig_df, data_se){
  require(GSVA)
  # check signature arg is data frame
  if (!is.data.frame(sig_df)) {
    stop("Signature argument is not a dataframe.")
  }

  if (all(abs(sig_df$expression) != 1)){
    stop("Expression values of signature data frame are not -1 or 1 and thus incompatible.")
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
  if (typeof(data_se) == "double" && all(data_se %% 1 == 0)) {
    # If the data is of type integer, these are raw counts and follow Poisson
    scores <- gsva(data_se, sigs, kcdf="Poisson", verbose=F)
  } else if (typeof(data_se) == "double") {
    scores <- gsva(data_se, sigs, verbose=F)
  } else {
    stop("Expression dataset contains non-numerical data. Try pre-processing dataset before classification.")
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
  return(list("Up"=sorted_up, "Down"=sorted_dn))
}

#' Generic classify function for classifying samples based on GSVA scores
#'
#' @param gsva_scores List containing a data frame (Up) of samples and their
#' corresponding GSVA scores when running GSVA on the up-regulated gene-set, and
#' a dataframe (Down) with GSVA scores for the same samples when running GSVA
#' on the down-regulated gene-set.
#' @param up_thresh Numerical vector containing the low and high threshold for
#' classification when checking consistency with only the up-regulated gene-set
#' of the gene signature
#' @param dn_thresh Numerical vector containing the low and high threshold for
#' classification when checking consistency with only the down-regulated
#' gene-set of the gene signature
#' @keywords internal
#'
#' @return data frame containing the pathway activity class labels for each
#' sample i.e. "Active", "Inactive" or "Uncertain".
#'
#' @examples
#' classes_df <- classify(gsva_scores, up_thresh, dn_thresh)
.classify <- function(gsva_scores, up_thresh.low, up_thresh.high, dn_thresh.low,
                     dn_thresh.high) {
  if (missing(up_thresh.high) || missing(dn_thresh.high) ||
      missing(up_thresh.low) || missing(dn_thresh.low)) {
    stop("Function requires thresholds for classifying samples using GSVA scores
    generated using the up-regulated and down-regulated gene-set of the gene
         signature.")
  }
  # create data frame containing samples and classes
  classes_df <- data.frame(sample=rownames(gsva_scores$Up),
                           class=vector(mode="character", nrow(gsva_scores$Up)))
  row.names(classes_df) <- classes_df$samples
  classes <- c()
  # Loop through sample rows which are ordered in descending order and checks
  # consistency with up-regulated and down-regulated parts of the gene signature
  classes_list <- sapply(rownames(gsva_scores$Up), function(i) {
    if ((gsva_scores$Up[i, 1] >= up_thresh.high) &&
        (gsva_scores$Down[i, 1] <= dn_thresh.low)) {
      classes[i] <- "Active"
    } else if ((gsva_scores$Up[i, 1] <= up_thresh.low) &&
               (gsva_scores$Down[i, 1] >= dn_thresh.high)) {
      classes[i] <- "Inactive"
    } else {
      classes[i] <- "Uncertain"
    }
    return(classes)
  })

  # convert list to vector
  classes <- factor(classes_list, levels=c("Active", "Inactive", "Uncertain"))

  # add vector to pathway df
  classes_df$class <- classes
  cat("Summary of sample classification based on pathway activity:\n")
  cat("--------------------------------------------------------------\n")
  cat("Number of samples in each pathway activity class:\n")
  print(table(classes_df$class))
  cat(sprintf("\nTotal number of samples: %d", nrow(classes_df)))
  classified <- classes_df[classes_df$class %in% c("Active", "Inactive"),]
  cat(sprintf("\nTotal number of samples classified: %d\n", nrow(classified)))
  return(classes_df)
}
