# scInt
scInt: accurate integration of multiple heterogeneous single-cell RNA-seq data sets by learning contrastive biological variation

scInt is a integration tool of exploring single-cell transcriptome in biological issues, especially in the comparison of single-cell samples from different biological conditions. If you have any suggestions or problems on this package, please contact Y. Zhou (yangz@stu.hit.edu.cn), or submit through this Github page directly.
# Installation
To run scInt R package, install from Github through ``devtools`` directly:
```R
install.packages('devtools')
library(devtools)
devtools::install_github("Jinsl-lab/scInt")
```
# Usage
## Data Integration
For integration task, scInt takes a list of gene by cell gene expression matrices as inputs and outputs the integrated low-dimensional data.

```R
library(scInt)
data <- create.scInt(dataset)
data <- normalize.Data(data)
data <- ident.cellIdentity(data)
data <- compute.Similarity(data)
data <- run.Integrate(data)
```

## Reference-based mapping
For mapping task, scInt takes a list of gene by cell gene expression matrices as inputs and outputs the integrated low-dimensional data.

```R
data <- create.scInt(dataset)
data <- normalize.Data(data)
data <- ident.cellIdentity(data)
data <- compute.Similarity(data)
data <- run.Integrate(data)
```

## Tutorials
For more details and full workflow of scInt inculing downstream analysis, please check [here]().

* For data Integration and reference-based mapping examples using simulated data sets, please see [here]().

* For data integration and subsequent downstream analysis on real data set, please see [here]().

# Reproducibility and data sets in the manuscript
The source codes and jupyter notebook scripts to reproduce the results in the manuscript are available on the [Github page](https://github.com/Jinsl-lab/scInt_reproducibility). The full information of used data sets can be find in the manuscript. The data sets are available for free in public data bases, and can also be downloaded in the [Zenodo repository]().

# Dependencies
scInt has been successfully installed and used on Windows, Linux and Mac OS. The dependencies including: Biobase, dplyr, irlba, Matrix, matrixStats, methods, rfunctions, Rcpp, RcppEigen, RcppArmadillo, RSpectra, Seurat (>= 4.1.0), uwot.
