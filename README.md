# Summary

  PathAnalyser is an user-friendly R package that provides functionality for assessing ER and 
  HER2 pathway activity in breast cancer transcriptomic datasets. 

# Table of Contents

- [Summary](#summary)
- [PathAnalyser workflow](#pathanalyser-workflow)
- [Overview of PathAnalyser functionality](#overview-of-pathanalyser-functionality)
- [Installation](#installation)
    - [Dependencies](#dependencies)
    - [Install PathAnalyser with devtools](#install-pathanalyser-with-devtools)
    - [Install PathAnalyser from source](#install-pathanalyser-from-source)
- [If you wish to know more](#if-you-wish-to-know-more)
  
# PathAnalyser workflow

![PathAnalyser workflow](./vignettes/pathway_workflow.png?raw=true)

# Overview of PathAnalyser functionality

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

# Installation

## Dependencies

PathAnalyser needs the following:
- **R** (tested on version 4.1.1)
- **An internet connection**
- **The following R libraries:** (The number is the version tested during development)

```` 
   VennDiagram (1.7.1)     futile.logger (1.4.3) 
   reader (1.0.6)          NCmisc (1.1.6)        
   ggplot2 (3.3.5)         reshape2 (1.4.4)
   edgeR (3.34.1)          limma (3.48.3)        
   pROC (1.18.0)           GSVA (1.40.1)
````
**note:** The package is platform-independent; it runs on multiple operating systems.

To install the dependencies you can use the following command in R :

````
# If not already installed
install.packages("BiocManager")

BiocManager::install(c("GSVA", "pROC", "edgeR", "reshape2", "ggplot2","limma", "reader", "VennDiagram", "NCmisc", "futile.logger"),                           dependencies = TRUE)
````

## Install PathAnalyser with devtools

The easiest way to get PathAnalyser is to install it directly from R using “devtools”:

````
install.packages("devtools")
library(devtools)
install_github(repo = "https://github.com/ozlemkaradeniz/PathAnalyser", dependencies = TRUE)
library(PathAnalyser)
````

## Install PathAnalyser from source

Alternatively you can clone the GitHub repository:

````
git clone git@github.com:ozlemkaradeniz/PathAnalyser.git
````

Then type the following in R

````
library(utils)
install.packages("./PathAnalyser/", repos = NULL, type = "source")
````

# If you wish to know more

Look in the vignette here:
http://ozlemkaradeniz.github.io/PathAnalyser/






