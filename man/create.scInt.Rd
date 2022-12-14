% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scInt.R
\name{create.scInt}
\alias{create.scInt}
\title{Create a new scInt object}
\usage{
create.scInt(
  raw.data,
  meta = NULL,
  filter = TRUE,
  min.features = 300,
  min.cells = 3,
  modal = c("single", "multi"),
  do.sparse = TRUE
)
}
\arguments{
\item{raw.data}{List of raw data matrices (gene by cell).}

\item{meta}{The meta data of input datasets. (default NULL)}

\item{filter}{Whether filter genes and cells. (default TRUE)}

\item{min.features}{Filter the cells expressed features less than min.features. (default 300)}

\item{min.cells}{Filter the features expressed in cells less than min.cells. (default 3)}

\item{modal}{Dataset modal.}

\item{do.sparse}{Whether use sparse format ('dgCMatirx' class). (default TRUE)}
}
\value{
\code{scInt} object.
}
\description{
This function initializes a scInt object of multiple single cell dataset, for further integration and analysis.
It takes a list of multiple gene-cell expression matrices as inputs.
}
