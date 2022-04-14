# PathAnalyser

# Summary

  PathAnalyser is a flexible and user-friendly R package that provides functionality for assessing ER and 
  HER2 pathway activity in breast cancer transcriptomic datasets by using a gene expression signature. Unlike 
  other available pathway assessment packages, which do not distinguish between up-regulated and down-regulated gene sets, 
  PathAnalyser classifies samples as "active" or "inactive" for a pathway, only if there is a consensus in the
  evidence of expression consistency with both the up-regulated and down-regulated gene sets of the gene signature.
 

# Table of Contents

- [Summary](#summary)
- [PathAnalyser workflow](#pathanalyser-workflow)
- [Installation](#installation)
    - [Dependencies](#dependencies)
    - [Install PathAnalyser from source](#install-pathanalyser-from-source)
- [Typical usage](#typical-usage)
- [If you wish to know more](#if-you-wish-to-know-more)
  
# PathAnalyser workflow
The typical workflow for using the package is outlined below:

<img src="vignettes/workflow_flowchart.png" width="600">

**Note**: Assessment of classification is optional and can only occur if true pathway class labels ("Active", "Inactive", "Uncertain")
are available for the transcriptomic dataset.

# Installation

## Dependencies

PathAnalyser needs the following:
- **R** (tested on version 4.1.1)
- **An internet connection for downloading the package bundle**
- **The following R libraries:** (The number is the version tested during development)

```` 
   VennDiagram (1.7.1)     futile.logger (1.4.3) 
   reader (1.0.6)          NCmisc (1.1.6)        
   ggplot2 (3.3.5)         reshape2 (1.4.4)
   edgeR (3.34.1)          limma (3.48.3)        
   pROC (1.18.0)           GSVA (1.40.1)
````
**Note:** The package is platform-independent; it runs on multiple operating systems.

All dependencies should be installed together with the PathAnalyser package,
however, they can be installed separately. To install all required CRAN 
dependencies of PathAnalyser, type the following in R:
```{r eval=F}
install.packages(c("ggfortify", "ggplot2", "glue", "lifecycle", "cli", "plotly",
                   "reader", "pROC", "reshape2", "rlang", "VennDiagram", "withr"
                   ))

```
All Bioconductor dependencies can be installed by typing the following in R:
```{r eval=F}
BiocManager::install(c("edgeR", "limma"))
```

## Install PathAnalyser from source

You can download the latest source version from the latest release section by clicking on this [link](https://github.com/ozlemkaradeniz/PathAnalyser/releases).

Then to install this local source package type the following in R:

````
library(utils)
install.packages("PathAnalyser", repos = NULL, type = "source")
````
# Typical usage
For typical usage, please read the demo script and use the provided supplementary data.

# If you wish to know more

Look in the vignette here:
http://ozlemkaradeniz.github.io/PathAnalyser/






