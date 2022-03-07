
#installing the required packages
#BiocManager::install("GSVA")

# loading the required packages
require(GSVA)
require(DESeq2)

# clear workspace
#rm(list=ls())

# turn off any open graphics windows
#graphics.off()



# load data
#data <- read.table("TCGA_unannotated.txt", sep="\t", header=T)

#ask the type of expression data from user, whether "microarray" or "RNA-Seq"
LOOP <- 1
while (LOOP == 1) {
  cat("\nEnter type of expression data (Microarray(M)/RNA-Seq(R): ")
  datatype <- scan("", what=double(0), nlines=1) 
  if (length(datatype)=0){
    LOOP=0 #terminate loop
    cat("Loop terminated\n")
    
  } elif (datatype=="M"){
    microarray=datatype
    
  } elif (datatype=="R"){
    microarray=datatype
  }

#read the gene expression table from the user
  
#ask the number of files for gene signature data provided by the user  
  
LOOP=2
while(LOOP==2){
  cat ("Enter the number of signature files you are going to provide: ")
  n=scan("", what=double(0),nlines=1)
  if (length(n)=0){
    LOOP=0 #terminate loop
    cat("Loop terminated\n")
    
  } elif (n>4){
    cat("Limit for number of signature files exceeded!\n")
  
}
  

    
read_expression_data=function(gene_exp_data){
  
  filepath=path(gene_exp_data)
  if (endsWith("filepath", ".tsv"))==TRUE{
    exp_data=readLines(gene_exp_data)
    exp_data=as.data.frame(exp_data)
    
  } else{
    message("Gene Expression file is not a .tsv file!")
  }


#read the gene signature data from the user
  
read_signature=function(up_reg_sig, down_reg_sig){
    filepath1=path(up_reg_sig)
    filepath2=path(down_reg_sig)
    if (endsWith("filepath1", ".txt")==T && endsWith("filepath2", ".txt")==T){
      up_sig_data=readlines(up_reg_sig)
      down_sig_data=readlines(down_reg_sig)
      up_sig_data=as.data.frame(up_sig_data)
      down_sig_data=as.data.frame(down_sig_data)
      
    } else {
      message("The gene signature data are not in text file format!")
    }

  }
#Quality and Preprocessing
QualityCheck=function{
  #Replacing outliers with NA
  for ( i in col(gene_exp_data){
    boxplot(i)
    outliers <- boxplot(i, plot = FALSE)$out
    df[i %in% outliers, colnames(gene_exp_data)] = NA
  }
  #Replacing NA with mean value of the column
  for(i in 1:ncol(gene_exp_data)){
    gene_exp_data[is.na(gene_exp_data[,i]), i] <- mean(data[,i], na.rm = TRUE)
  }
  
}


##Microarray

#scaling
norm.rma2.quantile <- threestep(gene_exp_data, background.method="RMA.2", normalize.method="quantile",summary.method="median.polish")
colnames(geneSigUp) <- "gene"
colnames(geneSigDn) <- colnames(geneSigUp)

filtered_up_data <- data[rownames(data) %in% geneSigUp[,1],]
filtered_dn_data <- data[rownames(data) %in% geneSigDn[,1],]
filtered_HER2 <- rbind(filtered_up_data, filtered_dn_data)
filtered_HER2_mat <- data.matrix(filtered_HER2)

# variance stabilising transformation 
matNorm <- vst(matData)

signature <- rbind(geneSigUp, geneSigDn)

# run each gene set (upregulated and then down regulated)
sigs <- list()
sigs$up <- geneSigUp$gene
sigs$dn <- geneSigDn$gene

scores <- gsva(filtered_HER2_mat, sigs)
# transpose the matrix
tScores <- t(scores)
# create an up and down regulated matrix
scoresUp <- as.data.frame(tScores[,1], row.names=rownames(tScores))
colnames(scoresUp) <- "up"

scoresDn <- as.data.frame(tScores[,2])
colnames(scoresDn) <- "dn"
scoresUp
scoresDn

# sort rows
sortedUp <- scoresUp[order(scoresUp$up, decreasing=T), ,drop=F]
sortedDn <- scoresDn[order(scoresDn$dn, decreasing=T), ,drop=F]

# compute quartiles
quantUp <- quantile(sortedUp$up)
quantDn <- quantile(sortedDn$dn)


# HER 2 (pathway)
pathway <- data.frame(samples=rownames(sortedUp), HER2=vector(mode="character", nrow(sortedUp)))
row.names(pathway) <- pathway$samples

# loop through genes if a gene has high score for upregulated and has 
# low score in down regulation, then that sample is positive
pathwayActivity <- lapply(rownames(sortedUp), function(i) {
  # if the sample is in the upper quartile of the GSVA score for upregulated geneset
  # and in the 
  if (sortedUp[i, 1] > quantUp[4] && sortedDn[i, 1] < quantDn[2]) {
    pathway[i, 2] <- "high"
    # quantUp[2] = 25%, quantUp[4] = 75%
  } else if (sortedUp[i, 1] < quantUp[2] && sortedDn[i, 1] > quantDn[4]){
    pathway[i, 2] <- "low"
  } else {
    pathway[i, 2] <- "uncertain"
  }
})

# convert list to vector
activity <- unlist(pathwayActivity)
activity <- factor(activity, levels=c("high", "low", "uncertain"))

# add vector to pathway df
pathway$HER2 <- activity

summary(pathway$HER2)
for (i in 1:nrow(pathway)) {
  pathway[i,1] <- gsub("\\.", "-", pathway[i,1])
}

rownames(pathway)

#write.table(pathway, file="pathwayActivity.txt", row.names=F, quote=F)


# assess normality
resNorm <- gsva(matNorm, sigs)
samples <- c(1:50)
#plot(samples, res[1,1:50], type='l')
#lines(samples, resNorm[1, 1:50], col="red")
