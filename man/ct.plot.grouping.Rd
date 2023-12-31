% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellTree.R
\name{ct.plot.grouping}
\alias{ct.plot.grouping}
\title{Plot cell tree with grouping information}
\usage{
ct.plot.grouping(tree, file.output = NULL, show.labels = FALSE,
  force.recompute.layout = FALSE, height = 20, width = 10,
  vertebrae.distance = 0, backbone.vertex.size = 0, vert.vertex.size = 0)
}
\arguments{
\item{tree}{An \link[igraph]{igraph} tree, as returned by \code{\link{compute.backbone.tree}}}

\item{file.output}{String (optional). Path of a file where the plot should be saved in PDF format (rendered to screen if omitted).}

\item{show.labels}{Boolean (optional). Whether to write each cell's row number next to its vertex.}

\item{force.recompute.layout}{Boolean (optional). If set to \code{TRUE}, recomputes the graph's layout coordinates even when present.}

\item{height}{Numeric (optional). Height and width (in inches) of the plot.}

\item{width}{Numeric (optional). Height and width (in inches) of the plot.}

\item{vertebrae.distance}{Numeric (optional). If non-zero: forces a specific plotting distance (in pixels) between backbone cells and related peripheral cells (`vertebrae').}

\item{backbone.vertex.size}{Numeric (optional). Diameter (in pixels) of backbone and vertebrae cell vertices.}

\item{vert.vertex.size}{Numeric (optional). Diameter (in pixels) of backbone and vertebrae cell vertices.}
}
\value{
An updated \link[igraph]{igraph} object with \code{x} and \code{y} vertex coordinate attributes.
}
\description{
Plots a backbone tree (or MST) that was computed with \code{\link{compute.backbone.tree}}, displaying each cell's grouping.
}
\examples{
# Load pre-computed LDA model for skeletal myoblast RNA-Seq data from HSMMSingleCell package:
data(HSMM_lda_model)

# Recover sampling time (in days) for each cell:
library(HSMMSingleCell)
data(HSMM_sample_sheet)
days.factor = HSMM_sample_sheet$Hours
days = as.numeric(levels(days.factor))[days.factor]

# Compute near-optimal backbone tree:
b.tree = compute.backbone.tree(HSMM_lda_model, days)
# Plot resulting tree with sampling time as a vertex group colour:
ct.plot.grouping(b.tree)
}

