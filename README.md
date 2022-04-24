# <img src="pathanalyser.png" alt="PathAnalyser">

# Summary

PathAnalyser is a flexible and user-friendly R package that provides functionality for assessing ER and HER2 pathway activity in breast cancer transcriptomic datasets using a gene expression signature. Typically, gene signatures can be broadly classified into the following three categories:

1.  **gene-sets** - list of genes without information regarding the strength or direction of association with a phenotype
2.  **weighted gene lists** - lists of genes with numerical weights representing strength and direction of association with a phenotype
3.  **gene-signature containing only direction of association** with the phenotype (i.e. a list of up-regulated and down-regulated genes).

Several currently available packages / algorithms classify samples using gene-sets such as GSVA and GSEA or weighted gene lists such as the PAM50 algorithm. However, despite the third type of signature (direction-associated gene signatures) being reported by numerous publications, there is currently *no software tools available for classifying samples based on these signatures*. PathAnalyser addresses this need by providing functionality for classifying samples by pathway activity using these highly reported and widely available gene signatures.

# Table of Contents

-   [Summary](#summary)
-   [PathAnalyser workflow](#pathanalyser-workflow)
-   [Installation](#installation)
    -   [Dependencies](#dependencies)
    -   [Install PathAnalyser from source](#install-pathanalyser-from-source)
-   [Input File Formats](#input-file-formats)
-   [Quick Start](#quick-start)
    -   [Load PathAnalyser](#load-pathanalyser)
    -   [Reading data from input files](#reading-data-from-input-files)
        -   [Gene expression data set file](#gene-expression-data-set-file)
        -   [Gene signature files](#gene-signature-files)
    -   [QC and pre-processing of gene expression dataset](#qc-and-pre-processing-of-gene-expression-dataset)
    -   [Classification based on pathway activity](#classification-based-on-pathway-activity)
    -   [Visualising classification](#visualising-classification)
    -   [Classification evaluation with true labels (optional)](#classification-evaluation-with-true-labels-optional)
-   [Acessing help](#accessing-help)
-   [Questions, bug reports or issues](#questions-bug-reports-or-issues)
-   [If you wish to know more](#if-you-wish-to-know-more)

# PathAnalyser workflow

The typical workflow for using the package is outlined below: <img src="vignettes/workflow_flowchart.png" width="600"/>

**Note**: Assessment of classification is optional and can only occur if true pathway class labels (e.g. "positive", "negative" or "uncertain") are available for the transcriptomic dataset.

# Installation

## Dependencies

PathAnalyser needs the following:

-   **R** (tested on version 4.1.1)
-   **GSVA** Bioconductor package (1.42.0) required by classification algorithm
-   **The following R libraries:** (The number is the version tested during development)

<!-- -->

       GSVA (1.40.1)           reader (1.0.6)
       edgeR (3.34.1)          ggplot2 (3.3.5)
       limma (3.48.3)          reshape2 (1.4.4)
       plotly (4.10.0)       
       

**Note:** The package is platform-independent; it was developed and runs on multiple operating systems (Windows, MacOS, Linux).

All dependencies should be installed together with the PathAnalyser package, however, they can be installed separately. To install all required CRAN dependencies of PathAnalyser, type the following in R:

```{r}
install.packages(c("ggfortify", "ggplot2", "plotly", "reader", "reshape2"))
```

All Bioconductor dependencies can be installed by typing the following in R:

```{r}
# if not previously installed
install.packages("BiocManager") 
BiocManager::install(c("GSVA","edgeR", "limma"))
```

## Install PathAnalyser from source

You can download the latest source tarball file (ending in .tar.gz) from the [latest release section](https://github.com/ozlemkaradeniz/PathAnalyser/releases) on the [PathAnalyser GitHub repository page](https://github.com/ozlemkaradeniz/PathAnalyser).

Then to install this local source package type the following in R:

```{r}
library(utils)
install.packages("PathAnalyser_1.0.0.tar.gz", repos = NULL, type = "source", dependencies=TRUE)
```

For instructions for more advanced installations (e.g. for developers), please consult the vignette.

# Input file formats

PathAnalyser can read two types of input data files:

1.  A **gene expression data set file** containing a table with sample names or IDs as column names, and gene names as row names<br/> An example gene expression data set file:

    <img src="vignettes/expr_dataset_example.png" alt="example expression dataset" width="400"/>

2.  **Gene signature files** (for up-regulated gene-set and another for down-regulated gene-set) which must be provided in either the [gene-set file format (GRP)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GRP:_Gene_set_file_format_.28.2A.grp.29) (shown below) or [gene matrix transposed file format (GMT)](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29). Currently, PathAnalyser requires 2 gene signatures files (one for the up-regulated and one for the down-regulated gene set) for a gene signature of a given pathway.<br/> An example gene signature file for the up-regulated gene-set of a gene signature:

    <img src="vignettes/HER2_up.png" alt="example gene signature file" width="600"/>

# Quick Start

## Load PathAnalyser

Once the package is installed, to start using PathAnalyser simply load the PathAnalyser package in R:

```{r}
library(PathAnalyser)
```

## Reading data from input files

### Gene expression data set file

To read a gene expression data set file (tab or comma value separated files i.e. files with extension .tsv/.csv/.txt) into a matrix format for classification, type the following in R:

```{r}
data_mat <- read_expression_data("gene_expr.txt")
```

### Gene signature files

To read the two signature files comprising the up-regulated gene set and down-regulated gene set of the gene signature, type the following in R using `up_sig_file` and `down_sig_file` parameters for the up-regulated and down-regulated gene-set of the signature respectively:

```{r}
sig_df <- read_signature(up_sig_file="up_gene_sig.grp", down_sig_file="down_gene_sig.grp") 
```

## QC and pre-processing of gene expression dataset

The classification functions of PathAnalyser require the expression dataset to be normalised.<br/>

**For unnormalised RNA-seq data only**, perform logCPM transformation on an unnormalised (raw count-containing) gene expression matrix (`data_mat`) generated by `read_expression_data` using `log_cpm_transform`:

```{r}
norm_data <- log_cpm_transform(data_set)
```

**Note:** Microarray datasets must be normalised prior to performing classification using PathAnalyser, as the package currently does not contain functionality for normalising microarray datasets.<br/>

For further quality control and data pre-processing including filtering genes from the gene expression matrix that are not present in the gene signature data frame, or those genes lacking expression values in \< 10% of the total number of samples call the `check_signature_vs_dataset` with the logCPM transformed gene expression matrix (`norm_data`) and gene signature data frame (`sig_df`):

```{r}
norm_data <- check_signature_vs_dataset(norm_data, sig_df)
```

**Note:** Genes with multiple names (typically found in microarray data sets), often delimited with "///", will be dropped from the gene expression matrix, regardless of the presence of one of these gene names in the gene signature.

## Classification based on pathway activity

Pathway-based classification using percentile thresholds, can be performed by using the classify_gsva_percent function with a normalised gene expression matrix and gene signature data frame:

```{r}
# Using default percentile threshold (quartile = 25%)
norm_data <- classify_gsva_percent(norm_data, sig_df)
```

A custom percentile threshold can be provided by the user for tuning the pathway-based classification, by adding the `percent_thresh` parameter:

```{r}
# Using a 50th percentile threshold (50%)
classes_df <- classify_gsva_percent(norm_data, sig_df, percent_thresh=50)
```

The generated output (`classes_df`) of the classification function is a data frame containing samples names as the first column and the predicted activity class for a given pathway as the second column ("Active", "Inactive", "Uncertain").

## Visualising classification

An interactive PCA plot for visualising the pathway-based classification of samples can be achieved by using the `classes_PCA` function with the normalised expression matrix (`norm_data`), the data frame produced by the `classify_gsva_percent` function (`classes_df`) and the pathway of interest:

```{r}
classes_PCA(norm_data, classes_df, pathway = "ER")
```

![An example PCA plot for samples classified by ER pathway activity](vignettes/pca_plot.png)

## Classification evaluation with true labels (optional)

**If true pathway class labels are available for the classified dataset**, users can obtain evaluation metrics for the classification such as accuracy, sensitivity, specificity etc using the `calculate_accuracy` function with the following arguments:

-   `true_labels`: a data frame containing sample names as the first column and true pathway class labels ("positive", "negative", "uncertain") in a column named after the pathway of interest

-   `predicted_labels`: a data frame containing the same sample names in the first column, named "sample", and predicted pathway class labels ("Active", "Inactive" or "Uncertain") as the second column which is named "class"

-   `pathway`: pathway of interest (**Note:** this pathway name must be identical to the name of the column containing the true pathway activity labels.

    For example:

```{r}
confusion_matrix <- calculate_accuracy("Sample_labels.txt", classes_df, 
                                           pathway = "ER")
```

For further examples of using PathAnalyser in pathway-based classification analysis, please refer to the demo script (under /demo folder) and use the provided supplementary data, or read the vignette.

# Accessing help

To access help pages for any of the functions or built-in data provided by PathAnalyser, prefix the name of the function or data set with a question mark, e.g. to get additional information on the `read_signature` function, type the following in R:

```{r}
?read_signature
```

# Questions, bug reports or issues

For any questions, feature requests, bug reports or issues regarding the latest version of PathAnalyser, please use the "[issues](https://github.com/ozlemkaradeniz/PathAnalyser/issues)" tab located at the top-left of the GitHub repository page.

# If you wish to know more

If the PathAnalyser GitHub repository is public, look in the vignette here: <http://ozlemkaradeniz.github.io/PathAnalyser/>

Alternatively, the vignette can be accessed within R, by typing the following (after PathAnalyser has been successfully installed):

```{r}
browseVignettes(package = "PathAnalyser")
```
