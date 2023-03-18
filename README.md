# scInt
scInt: accurate integration of multiple heterogeneous single-cell RNA-seq data sets by learning contrastive biological variation

scInt is a integration tool of exploring single-cell transcriptome in biological issues, especially in the comparison of single-cell samples from different biological conditions. If you have any suggestions or problems on this package, please contact Y. Zhou (yangz@stu.hit.edu.cn), or submit through this GitHub page directly.

![Overview](https://github.com/Jinsl-lab/scInt/blob/main/tutorial/Overview.jpg)
# Installation
To run scInt R package, install from GitHub through ``devtools`` directly:
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

Selection of key parameters T and lambda:

```R
data <- compute.Similarity(data, select.T = TRUE, cand_Ts = seq(0.55, 0.8, 0.05))
data <- run.Integrate(data, select.lambda = TRUE, cand_lambdas = c(0.1, 1, 5, 10, 20))
```

## Reference-based mapping
```R
# reference integration
reference.data <- create.scInt(reference.dataset)
reference.data <- normalize.Data(reference.data)
reference.data <- ident.cellIdentity(reference.data)
reference.data <- compute.Similarity(reference.data)
reference.data <- run.Integrate(reference.data)

# mapping
query.data <- create.scInt(query.dataset)
query.data <- normalize.Data(query.data)
query.data <- ident.cellIdentity(query.data, vargenes = intersect(reference.data@vargenes, rownames(query.data@norm.data[[1]])))
query.data <- compute.Similarity(query.data)
query.data <- run.Integrate(query.data)
query.data <- run.Mapping(reference.data, query.data)
```

## Tutorials

* For data integration using simulated and real data sets, please see [here](https://github.com/Jinsl-lab/scInt/blob/main/tutorial/Simulation%26real_datasets.ipynb).

* For reference-based mapping on real data sets, please see [here](https://github.com/Jinsl-lab/scInt/blob/main/tutorial/Human_pancreas_mapping.ipynb).

# Reproducibility and data sets in the manuscript
The source codes and jupyter notebook scripts to reproduce the results in the manuscript are available on the [Github page](https://github.com/Jinsl-lab/scInt_reproducibility). The full information of used data sets can be find in the manuscript. The data sets are available for free in public data bases, and can also be downloaded in the [Zenodo repository](https://zenodo.org/record/7198744#.ZAvbO3ZByw6).

# Dependencies
scInt has been successfully installed and used on Windows, Linux and Mac OS (R version >= 4.0.2). The dependencies including: Biobase, dplyr, irlba, Matrix, matrixStats, methods, rfunctions, Rcpp, RcppEigen, RcppArmadillo, RSpectra, Seurat (>= 4.1.0), uwot.
