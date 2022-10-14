# scInt: accurate alignment of multiple heterogeneous single-cell RNA-seq datasets by learning contrastive biological variation
scInt is a integration tool of exploring single-cell transcriptome in biological issues, especially in the comparison of single-cell samples from different biological conditions. If you have any suggestions or problems on this package, please contact Y. Zhou (20b912010@stu.hit.edu.cn).
## Installation
To run scInt, install from this github directly:

```R
devtools::install_github("Jinsl-lab/scInt")
```

## Usage
scInt takes a list of gene by cell gene expression matrices as inputs and outputs the integrated low-dimensional data.

```R
data <- create.scInt(dataset)
data <- normalize.Data(data)
data <- ident.cellIdentity(data)
data <- compute.Similarity(data)
data <- run.Integrate(data)
```
