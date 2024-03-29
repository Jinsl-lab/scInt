% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scInt.R
\name{filter.cpca}
\alias{filter.cpca}
\title{Filter cell-cell similarities for each dataset i in cPCA space}
\usage{
filter.cpca(modi, norms, i, dims, T_th)
}
\arguments{
\item{modi}{The list of cell-cell similarities}

\item{norms}{The list of normalized matrices.}

\item{i}{The dataset i.}

\item{dims}{The dimension of cPCA.}

\item{T_th}{The (candidate) similarity threshold(s) for filtering cells in the cPCA space.}
}
\value{
List of filtered cell-cell similarities.
}
\description{
Filter cell-cell similarities for each dataset i in cPCA space
}
