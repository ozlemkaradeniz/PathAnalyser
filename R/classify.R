# Classifier function for classifying samples according to pathway activity
# @author AT

# load libraries
require(GSVA)

# classifies samples according to pathway activity using GSVA
classify <- function(sigDf, dataSe, thresh=0.05) {
  # create list of gene sets (up-regulated and then down regulated) for gsva
  sigs <- list()
  sigs$up <- sigDf[sigDf$expression == 1, 1]
  sigs$dn <- sigDf[sigDf$expression == -1, 1]

  # run GSVA on expression data using 2 gene sets (up-regulated and
  # down-regulated) of the signature
  if (type(dataSe) == "integer") {
    # If the data is of type integer, these are raw counts and follow Poisson
    scores <- gsva(dataSe, sigs, kcdf="Poisson")
  } else {
    scores <- gsva(dataSe, sigs)
  }


  # transpose the matrix
  scores <- t(scores)
  # create an up and down regulated matrix
  scoresUp <- as.data.frame(scores[,1], row.names=rownames(scores))
  colnames(scoresUp) <- "up"

  scoresDn <- as.data.frame(scores[,2])
  colnames(scoresDn) <- "dn"


  # sort GSVA scores in descending order for the uporegulated
  # and down-regulated set
  sortedUp <- scoresUp[order(scoresUp$up, decreasing=T), ,drop=F]
  sortedDn <- scoresDn[order(scoresDn$dn, decreasing=T), ,drop=F]

  # compute quantiles according to user-specified threshold
  threshUp <- quantile(sortedUp$up, probs=seq(0, 1, thresh))
  threshDn <- quantile(sortedDn$dn, probs=seq(0, 1, thresh))

  # HER 2 (pathway)
  classesDf <- data.frame(sample=rownames(sortedUp),
                          class=vector(mode="character", nrow(sortedUp)))
  row.names(classesDf) <- classesDf$samples

  # loop through genes if a gene has high score for up-regulated and has
  # low score for down-regulated genes, then that sample is positive for pathway
  # if reverse then negative, else uncertain
  classesList <- lapply(rownames(sortedUp), function(i) {
    if (sortedUp[i, 1] > threshUp[4] && sortedDn[i, 1] < threshDn[2]) {
      classesDf[i, 2] <- "Positive"
    } else if (sortedUp[i, 1] < threshUp[2] && sortedDn[i, 1] > threshDn[4]){
      classesDf[i, 2] <- "Negative"
    } else {
      classesDf[i, 2] <- "Uncertain"
    }
  })

  # convert list to vector
  classes <- unlist(classesList)
  classes <- factor(classes, levels=c("Positive", "Negative", "Uncertain"))

  # add vector to pathway df
  classesDf$class <- classes

  cat("Summary of sample classification based on pathway activity:\n")
  print(summary(classesDf$class))
  return(classesDf)
}
