# PathAnalyser

  PathAnalyser is an user-friendly R package that provides functionality for assessing ER and 
  HER2 pathway activity in breast cancer transcriptomic datasets. The package enables 
  classification of samples in transciptomic datasets according to pathway activity by 
  using a gene expression signature associated with a pathway, a list of genes, which 
  vary in expression depending on the pathway activity. These transcriptional 
  signatures have been shown to have clinical predictive value, for example, ER 
  and HER2 gene signatures have been associated with molecular sub-types, 
  prognosis and treatment reponse in breast cancer. Although this package was 
  original developed to assess ER and HER2 pathway activity in breast cancer 
  transcriptomic datasets, the package functionality could be applied to 
  transcriptomic datasets and gene signatures outside the context of breast 
  cancer. In this vignette, we will describe how to use the PathAnalyser package 
  with microarray and RNA-seq expression transcriptomic datasets and gene 
  expression signatures associated with a specific pathway activity.
  
## PathAnalyser workflow

![PathAnalyser workflow](./vignettes/pathway_workflow.png?raw=true)


## Overview of PathAnalyser functionality


The PathAnalyser package functionality can be sub-divided into 5 broad 
categories corresponding to a typical workflow for assessing pathway activity of 
samples:

1.  Input functions for expression data and gene signatures
2.  Quality control and pre-processing of input data
3.  Classification of samples by pathway activity
4.  Evaluation of classification
5.  Visualization

The features unique to PathAnalyser package are:

1) It takes into account multiple sample for the analysis. This is an addition to GSVA package features, where a single sample is
taken into account for analysis.

2) It classifies the IDs/symbols into active and inactive pathway. This classification is based on the two matrices representating up and down regulated expression of genes respectively. The matrices consists of GSVA scores derived after GSVA analysis. It also returns a prediction accuracy score and ROC curve for verification of the accuray of the prediction.

3) It returns a PCA plot for visualization of the categories (Active/Inactive) the gene symbols/IDs have been classified into.

## Installation

### Dependencies

PathAnalyser needs the following:
- **R** (tested on version 3.6)
- **An internet connection**
- **The following R libraries:** (The number is the version tested during development)
```` 
   VennDiagram (1.7.1)     futile.logger (1.4.3) 
   reader (1.0.6)          NCmisc (1.1.6)        
   ggplot2 (3.3.5)         reshape2 (1.4.4)
   edgeR (3.34.1)          limma (3.48.3)        
   pROC (1.18.0)           GSVA (1.40.1)
````
**note:** The package is platform-independent; it runs multiple operating systems.

To install the dependencies you can use the following command in R :
````
# If not already installed
install.packages("BiocManager")

BiocManager::install(c("GSVA", "pROC", "edgeR", "reshape2", "ggplot2","limma", "reader", "VennDiagram", "NCmisc", "futile.logger"),                           dependencies = TRUE)

````

### Install PathAnalyser with devtools

The easiest way to get PathAnalyser is to install it directly from R using “devtools”:
````
install.packages("devtools")
library(devtools)
install_github(repo = "https://github.com/a-thind/PathAnalyser", dependencies = TRUE)
library(PathAnalyser)
````

### Install PathAnalyser from source

Alternatively you can clone the GitHub repository:
````
git clone git@github.com:a-thind/PathAnalyser.git
````
Then type the following in R
````
library(utils)
install.packages("./PathAnalyser/", repos = NULL, type = "source")
````

## Preparing the input

### Gene expression datasets
PathAnalyser also contains two built-in gene expression matrices each containing
RNA-seq raw read counts for primary breast cancer samples obtained from 20 
individuals (cases). Data for these matrices were obtained from The Cancer 
Genome Atlas ([TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)).
Each expression matrix contains 20,124 genes.
```{r}
data("ER_data_se1")
data("HER2_data_se1")

### ER dataset 1
The ER data set `ER_data_se1` contains RNASeq raw read counts for 20 primary 
breast cancer samples, 10 of which have ER pathway activity (ER+) and 10 which 
have inactive ER pathway activity (ER-).
```{r}
dim(ER_data_se1)
# Expression data for first 6 genes
head(ER_data_se1)
```

### HER2 dataset 1
Similarly, the HER2 data set `HER2_data_se1` contains RNASeq raw read counts for 
20 primary breast cancer samples, 10 of which have HER2 (ERBB2) pathway activity 
(HER2+) and 10 which have inactive HER2 pathway activity (HER2-).

```{r}
dim(HER2_data_se1)
# Expression data for first 6 genes
head(HER2_data_se1)
```

For more information you can visit vignete file.
https://github.com/a-thind/PathAnalyser/blob/main/vignettes/PathAnalyser.Rmd





