% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellTree.R
\docType{data}
\name{HSMM_lda_model}
\alias{HSMM_lda_model}
\title{Pre-computed LDA model for HSMMSingleCell data}
\format{An object of class \code{list} of length 3.}
\usage{
HSMM_lda_model
}
\value{
LDA model obtained by running \code{\link{compute.lda}} on the \pkg{HSMMSingleCell} package's data. This demo model can be used with all functions in this package that require a fitted LDA model, such as \code{\link{compute.backbone.tree}}.

See examples for \code{\link{compute.lda}} for more details.
}
\description{
Pre-computed LDA model for HSMMSingleCell data
}
\seealso{
\link{compute.lda}
}
\keyword{data}

