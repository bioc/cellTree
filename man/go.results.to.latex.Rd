% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellTree.R
\name{go.results.to.latex}
\alias{go.results.to.latex}
\title{LaTeX output for Gene Ontology Enrichment results}
\usage{
go.results.to.latex(go.results, tex.file.name = "go_terms.tex",
  topic.colors = rainbow(length(go.results$all)))
}
\arguments{
\item{go.results}{GO Enrichment result list object, such as returned by \code{\link{compute.go.enrichment}}.}

\item{tex.file.name}{String (optional). If not \code{NULL}: name of the file to save to.}

\item{topic.colors}{RGB colour vector (optional). Colors to use for each topic.}
}
\value{
GO enrichment results in LaTeX format
}
\description{
Outputs or writes result tables of Gene Ontology enrichment testing in LaTeX format.
}
\examples{
# Load pre-computed LDA model for skeletal myoblast RNA-Seq data from HSMMSingleCell package:
data(HSMM_lda_model)

\donttest{
# Load GO mapping database for 'homo sapiens':
library(org.Hs.eg.db)
# Compute GO enrichment sets for each topic:
go.results = compute.go.enrichment(HSMM_lda_model, org.Hs.eg.db, bonferroni.correct=TRUE)

# Output LaTeX tables for GO results
latex.table = go.results.to.latex(go.results, tex.file.name='go_results.tex')

# [Optional] compile LaTeX to PDF (need to have LaTeX binaries installed):
library('tools')
texi2pdf('go_results.tex')
library('Biobase')
openPDF('go_results.pdf')
}
}

