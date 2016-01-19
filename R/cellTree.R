#' Inference and visualisation of Single-Cell RNA-seq Data data as a hierarchical tree structure
#' 
#' This packages computes a Latent Dirichlet Allocation (LDA) model of single-cell RNA-seq data and build a compact tree modelling the relationship between individual cells over time or space.
#'
#' A typical use-case will require you to run \link{compute.lda} on your expression data, to fit an LDA model, followed by \link{compute.backbone.tree} to generate a tree structure from the LDA model.
#'
#' Plotting functions \link{ct.plot.grouping} and \link{ct.plot.topics} will then help you plot the tree structure with the desired annotations, while function \link{cell.ordering.table} will output an ordered table of cells, ranked by their position in the tree.
#'
#' To get further information on each topic, you can run Gene Ontology enrichment tests, using \link{compute.go.enrichment}, plot the result tables as a graph using \link{ct.plot.go.dag} or render it with LaTeX, using \link{go.results.to.latex}.
#'
"_PACKAGE"

#' Pre-computed LDA model for HSMMSingleCell data
#' 
#' @return LDA model obtained by running \code{\link{compute.lda}} on the \pkg{HSMMSingleCell} package's data. This demo model can be used with all functions in this package that require a fitted LDA model, such as \code{\link{compute.backbone.tree}}.
#'
#' See examples for \code{\link{compute.lda}} for more details.
#' 
#' @keywords data
#' @seealso \link{compute.lda}
"HSMM_lda_model"

#' LDA model inference
#'
#' This function fits a Latent Dirichlet Allocation (LDA) to single-cell RNA-seq data.
#'
#' @param data A matrix of (non-negative) RNA-seq expression levels where each row is a gene and each column is the cell sequenced.
#' @param method LDA inference method to use. Can be any unique prefix of `maptpx', `Gibbs' or `VEM' (defaults to `maptpx')
#' @param k.topics Integer (optional). Number of topics to fit in the model. If \code{method} is `maptpx', \code{k.topics} can be a vector of possible topic numbers and the the best model (evaluated on Bayes factor vs a null single topic model) will be returned.
#' @param log.scale Boolean (optional). Whether the data should be log-scaled.
#' @param sd.filter Numeric or \code{FALSE} (optional). Standard-deviation threshold below which genes should be removed from the data (no filtering if set to \code{FALSE}).
#' @param tot.iter,tol Numeric parameters (optional) forwarded to the chosen LDA inference method's contol class.
#' @return A LDA model fitted for \code{data}, of class \link[topicmodels]{LDA-class} (for methods 'Gibbs' or 'VEM') or \link[maptpx]{topics} (for 'maptpx')
#' @details Latent Dirichlet allocation (LDA) is a generative model that allows sets of observations to be explained by unobserved groups (topics) that explain why some parts of the data are similar [Blei, 2003]. Each topic is modelled as a (Dirichlet) distribution over observations and each set of observations is also modelled as a (Dirichlet) distribution over topics.
#' In lieu of the traditional NLP context of word occurence counts in documents, our model uses RNA-seq observation counts in single cells.
#' Three separate LDA inference methods can be used at the moment:
#' \itemize{
#'		\item{Gibbs} uses Collapsed Gibbs Sampling method (implemented by Xuan-Hieu Phan and co-authors in the \pkg{topicmodels} package [Phan, 2008]) to infer the parameters of the Dirichlet distributions for a given number of topics. It gives high accuracy but is very time-consuming to run on large number of cells and genes.
#'		\item{VEM} uses Variational Expectation-Maximisation (as described in [Hoffman, 2010]). This method tends to converge faster than Gibbs collapsed sampling, albeit with lower accuracy.
#'		\item{maptpx} uses the method described in [Taddy, 2011] and implemented in package \pkg{maptpx} to estimate the parameters of the topic model for increasing number of topics (using previous estimates as a starting point for larger topic numbers). The best model (/number of topics) is selected based on Bayes factor over the Null model. Although potentially less accurate, this method provides the fastest way to train and select from a large number of models, when the number of topics is not well known.
#' }
#' 
#' When in doubt, the function can be ran with its default parameter values and should produce a usable LDA model in reasonable time (using the `maptpx' inference method). The model can be further refined for a specific number of topics with slower methods.
#' While larger models (using large number of topics) might fit the data well, there is a high risk of overfitting and it is recommended to use the smallest possible number of topics that still explains the observations well. Anecdotally, a typical number of topics for cell differentiation data (from pluripotent to fully specialised) would seem to be around 4 or 5. 
#	
#' @seealso \link[topicmodels]{LDA}, \link[maptpx]{topics}, \link[topicmodels]{LDA_Gibbscontrol-class}, \link[topicmodels]{CTM_VEMcontrol-class}
#' @example examples/compute.lda.example.R
#' @references
#' \itemize{
#' 	\item Blei, Ng, and Jordan. ``Latent dirichlet allocation.'' the Journal of machine Learning research 3 (2003): 993-1022.
#'    \item Hoffman, Blei and Bach (2010). ``Online Learning for Latent Dirichlet Allocation.'' In J Lafferty, CKI Williams, J Shawe-Taylor, R Zemel, A Culotta (eds.), Advances in Neural Information Processing Systems 23, pp. 856-864. MIT Press, Cambridge, MA.
#'    \item Hornik and Grün. ``topicmodels: An R package for fitting topic models.'' Journal of Statistical Software 40.13 (2011): 1-30.
#'		\item Phan, Nguyen and Horiguchi. ``Learning to classify short and sparse text & web with hidden topics from large-scale data collections.'' Proceedings of the 17th international conference on World Wide Web. ACM, 2008.
#'		\item Taddy. ``On estimation and selection for topic models.'' arXiv preprint arXiv:1109.4518 (2011).
#' }
#'
#' @import slam
#' @export

# Disabling following imports to avoid collisions:
# import topicmodels
# import maptpx

compute.lda <- function(data, method="maptpx", k.topics=if(method=="maptpx") 2:15 else 4, log.scale=TRUE, sd.filter=0.5, tot.iter=if(method=="Gibbs") 200 else 1000000, tol = if(method=="maptpx") 0.05 else 10^-5) {
	.library.or.stop('slam')

	methods = c("maptpx", "VEM", "Gibbs")
	idx = pmatch(method, methods)
	if(is.na(idx)) stop(paste0("Unknown LDA inference method. Must be unique prefix of one of: ", paste0(methods, collapse=", ")))
	method = methods[idx]
	cat(paste0('Computing LDA model using: ', method, '\n'))

	data = .normalise.data(data, log.scale, sd.filter)
   
	if(method == "maptpx") {
		.library.or.stop('maptpx', "to use the 'maptpx' LDA inference method")

		# data.int = apply(data*100, 2, as.integer)
		# rownames(data.int) = rownames(data)
		# data.sparse = as.simple_triplet_matrix(t(data.int))
		data.sparse = slam::as.simple_triplet_matrix(t(data))
		
		lda.results = maptpx::topics(data.sparse, K=k.topics, tol=tol)
		
		cat(paste0("Selected k = ", lda.results$K, " topics\n"))
	} else {
		.library.or.stop('topicmodels', "to use the 'Gibbs' or 'VEM' LDA inference methods")
		if(length(k.topics) > 1) stop(paste0("At the moment, inference method '", method, "' only supports a single topic number value."))
		
		data.int = apply(data*100, 2, as.integer)
		rownames(data.int) = rownames(data)
		data.sparse = as.simple_triplet_matrix(t(data.int))
	
		if(method == "Gibbs")
			lda.results = topicmodels::LDA(data.sparse, vocab=rownames(data.int), k = k.topics, method = "Gibbs", control = list(keep = 10, verbose = 10, iter = tot.iter, burnin = 30))
		else
			lda.results = topicmodels::LDA(data.sparse, vocab=rownames(data.int), k = k.topics, method = "VEM", control = list(keep = 10, verbose = 10, var=list(iter.max=tot.iter, tol=tol), em=list(iter.max=tot.iter, tol=tol), alpha=0.1))
		
		cat(paste0("Model fit for k = ", lda.results@k, " topics\n"))
		
	}
		
	return (lda.results)
}

#' Cell Pairwise-Distance Matrix
#'
#' Computes the pairwise distance between cells, based on the topic histograms form a fitted LDA model.
#'
#' @inheritParams compute.backbone.tree
#' @return A square matrix of pairwise distance between cells in the input model. 
#' @details Distances between histograms are computed using the Chi-square distance
#' 
#' @example examples/get.cell.dists.example.R
#'
#' @export

get.cell.dists <- function(lda.results) {
	if(class(lda.results) != 'topics' & !is.list(lda.results)) {
		.library.or.stop('topicmodels', "to use this type of lda.results object")

		topic.distribs = lda.results@gamma
	} else {
		.library.or.stop('maptpx', "to use this type of lda.results object")

		topic.distribs = lda.results$omega
	}
	
	dists = apply(topic.distribs, 1, function(x.1) apply(topic.distribs, 1, function(x.2) .chi.square.dist(x.1, x.2)))
	
	return(dists)
}

#' Backbone Tree construction
#'
#' Builds a `backbone tree' from a fitted LDA model.
#'
#' @param lda.results A fitted LDA model, as returned by \code{\link{compute.lda}}
#' @param grouping An (optional) vector of labels for each cell in the \code{lda.results} object. E.g. a sampling times (numeric) or tissue categories.
#' @param start.group.label If a \code{grouping} parameter is provided, you can optionally specify the starting group. If no \code{start.group.label} is specified and the \code{grouping} vector is numeric, the lowest value will automatically be selected. Otherwise, the group with lowest mean-squared-distance between cells is selected.
#' @param absolute.width Numeric (optional). Distance threshold below which a cell vertex is considered to be attached to a backbone vertex (see paper for more details). By default, this threshold is computed dynamically, based on the distance distribution for each branch.
#' @param width.scale.factor Numeric (optional). A scaling factor for the dynamically-computed distance threshold (ignored if \code{absolute.width} is provided). Higher values will result in less branches in the backbone tree, while lower values might lead to a large number of backbone branches.
#' @param outlier.tolerance.factor Numeric (optional). Proportion of vertices, out of the total number of vertices divided by the total number of branches, that can be left at the end of the backbone tree-building algorithm.
#' @param rooting.method String (optional). Method used to root the backbone tree. Must be either NULL or one of: `longest.path', `center.start.group' or `average.start.group'. `longest.path` picks one end of the longest shortest-path between two vertices. `center.start.group' picks the vertex in the starting group with lowest mean-square-distance to the others. `average.start.group' creates a new artificial vertex, as the average of all cells in the starting group. If no value is provided, the best method is picked based on the type of grouping and start group information available.
#' @param only.mst If \code{TRUE}, returns a simple rooted minimum-spanning tree, instead of a backbone tree.
#' @param grouping.colors (Optional) vector of RGB colors to be used for each grouping.
#' @param merge.sequential.backbone (Optional) whether to merge sequential backbone vertices that are close enough. This will produce a more compact backbone tree, but at the cost of extra computing time.
#' @return A \link[igraph]{igraph} object with either a minimum rooted spanning-tree (if \code{only.mst} is \code{TRUE}) or a quasi-optimal backbone tree connecting all input cells. Cell topic distribution, distances and branch order are added as vertex/edge/graph attributes.
#' @details In order to easily visualise the structural and temporal relationship between cells, we introduced a special type of tree structure dubbed `backbone tree', defined as such:
#'
#' Considering a set of vertices \eqn{V} and a distance function over all pairs of vertices: \eqn{d: V \times V \rightarrow \strong{R}^+}{d: V × V -> R+}, we call \emph{backbone tree} a graph, \eqn{T} with backbone \eqn{B}, such that:
#' \itemize{
#'   \item \eqn{T} is a tree with set of vertices \eqn{V} and edges \eqn{E}.
#'   \item \eqn{B} is a tree with set of vertices \eqn{V_B \subseteq V}{V_B in V} and edges \eqn{E_B \subseteq E}{E_B in E}.
#'   \item All `vertebrae' vertices of \eqn{T}: \eqn{v \in V \setminus V_B}{v in V \ V_B} are connected by a single edge to the closest vertex in the set of backbone vertices \eqn{v^*_B \in V_B}{v*_B in V_B}. I.e: \eqn{v^*_B = argmin_{v_B \in V_B} d(v_B, v)}{v*_B = argmin_{v_B in V_B} d(v_B, v)}.
#'   \item For all vertices in \eqn{V \setminus V_B}{V \ V_B} are less than distance \eqn{\delta} to a vertex in the backbone tree \eqn{B}: \eqn{\forall v \in V \setminus V_B, \exists v_B \in V_B}{for all v in V \ V_B, there is a v_B in V_B} such that \eqn{d(v, v_b) \le \delta}{d(v, v_b) < \delta}.
#'}
#' In this instance, we relax the last condition to cover only `most' non-backbone vertices, allowing for a variable proportion of outliers at distance \eqn{> \delta}{> \delta} from any vertices in \eqn{V_B}.
#' 
#' We can then define the `optimal' backbone tree to be a backbone tree such that the sum of weighted edges in the backbone subtree \eqn{E_B} is minimal.

#' Finding such a tree can be easily shown to be NP-Complete (by reduction to the Vertex Cover problem), but we developed a fast heuristic relying on Minimum Spanning Tree to produce a reasonable approximation.
#'
#' The resulting quasi-optimal backbone tree (simply referred to as `the' backbone tree thereafter) gives a clear hierarchical representation of the cell relationship: the objective function puts pressure on finding a (small) group of prominent cells (the backbone) that are good representatives of major steps in the cell evolution (in time or space), while remaining cells are similar enough to their closest representative for their difference to be ignored. Such a tree provides a very clear visualisation of overall cell differentiation paths (including potential differentiation into sub-types).
#'
#' @example examples/compute.backbone.tree.example.R
#' @example examples/ct.plot.grouping.example.R
#'
#' @import igraph
#' @export

compute.backbone.tree <- function(lda.results, grouping = NULL, start.group.label = NULL, absolute.width = 0, width.scale.factor = 1.2, outlier.tolerance.factor = 0.1, rooting.method = NULL, only.mst = FALSE, grouping.colors = NULL, merge.sequential.backbone = FALSE) {
	.library.or.stop('igraph')
	
	### Sanity checks:
	if(is.list(grouping)) grouping = base::unlist(grouping)
	
	if(class(lda.results) != 'topics' & !is.list(lda.results)) {
		.library.or.stop('topicmodels', "to use this type of lda.results object")

		topic.distribs = lda.results@gamma
		cell.names = lda.results@documents
	} else {
		.library.or.stop('maptpx', "to use this type of lda.results object")

		topic.distribs = lda.results$omega
		cell.names = rownames(lda.results$omega)
	}
	
	if(! is.null(grouping) & (length(grouping) != length(cell.names))) stop("grouping argument must be NULL or same length as the number of cells in lda.results")

	### Give some room to recursive functions:
	save.option.expressions = getOption("expressions")
	options(expressions=length(cell.names)*5)
	
	grouping = .format.grouping(grouping, grouping.colors)
      
	start.group = .find.start.group(grouping, lda.results, start.group.label)	
   rooting.method = .pick.rooting.method(rooting.method, start.group)
   
   # Build rooted Minimum Spanning Tree
   
	if(rooting.method == "average.start.group") {
   	# add artificial node (mean of first group)
      topic.distribs = rbind(apply(topic.distribs[grouping$norm == start.group,], 2, mean), topic.distribs)
		root.v = 1
	}
	
	# can't use get.cell.dists since we may have added an artificial cell
	dists = apply(topic.distribs, 1, function(x.1) apply(topic.distribs, 1, function(x.2) .chi.square.dist(x.1, x.2)))

	g = graph.adjacency(dists, mode = "undirected", weighted = TRUE)
	min.span.tree = minimum.spanning.tree(g)
	
	if(rooting.method == "longest.path") {
      if(is.null(start.group)) {
         # doing our best to guess where the 'beginning' vertices are...
         shortest = shortest.paths(min.span.tree)
         longest.endpoints = c(which.max(apply(shortest, 1, max)), which.max(apply(shortest, 2, max)))
         # Have both ends of longest path. Picking the one closest to its neighbourhood:

         mean.dist.to.neighbours = sapply(longest.endpoints, function(v) { mean(dists[v, neighborhood(min.span.tree, 2, v)[[1]]]) })
         root.v = longest.endpoints[which.min(mean.dist.to.neighbours)]
      } else {
         shortest = shortest.paths(min.span.tree, which(grouping$norm == start.group))
         root.v = which(grouping$norm == start.group)[which.max(apply(shortest, 1, max))]
      }
  	} else if(rooting.method == "center.start.group") {
		#lowest mean squared-distance to rest of start.group
      idx = which.min(apply(dists[grouping$norm == start.group, grouping$norm == start.group]^2, 1, mean))
   	root.v = as.integer(V(min.span.tree)[grouping$norm == start.group][idx])
	}
	
	cat(paste0("Using root vertex: ", root.v, '\n'))
	
   if(only.mst) {
   	V(min.span.tree)$is.backbone = FALSE
   	V(min.span.tree)$is.root = FALSE
   	V(min.span.tree)$is.root[root.v] = TRUE
      
      cat("Returning Minimum Spanning Tree\n")
      
      return(.assign.tree.data(min.span.tree, rooting.method, grouping, topic.distribs, cell.names)
)
   }
   
   # Compute backbone tree from MST
	b.tree.dists = matrix(0, nrow(dists), nrow(dists))
	backbone = c(root.v)

	tot.branches = 0
	tot.remaining = vcount(min.span.tree)
	remaining = min.span.tree
	V(remaining)$start = FALSE
	width = NULL
   
	# Iteratively adding longest branches from MST:
	while(tot.remaining > outlier.tolerance.factor * vcount(min.span.tree)) {
		
		V(remaining)[nei(as.integer(backbone))]$start = TRUE
		V(remaining)[backbone]$start = FALSE
		remaining = remaining - E(remaining)[from(backbone)]
		
		remain.start = V(remaining)[V(remaining)$start]
		shortest = shortest.paths(remaining, remain.start)
		shortest[is.infinite(shortest)] = 0
		
		if(max(shortest) == 0) {
         cat("Remaining graph is fully disconnected: ")
         print(remaining)
         break
      }
		remaining.bottom.v = which.max(apply(shortest, 2, max))
		remaining.top.v = remain.start[which.max(shortest[,remaining.bottom.v])]
				
		branch = get.shortest.paths(remaining, remaining.top.v, remaining.bottom.v, output="both")	
		branch.v = as.integer(branch$vpath[[1]])
      if(length(branch.v) < 2)
         break # longest branch is a singleton
		
      tot.branches = tot.branches+1
      
      cat(paste0("Adding branch #", tot.branches, ":\n"))
      print(branch.v)
      
      if(is.null(width)) {
         if(absolute.width == 0) {
            # Setting width to first local minimum in multimodal density distrib of shortest dists to current path:
   			dist.density = density(dists[upper.tri(dists)])

            idx = which(diff(sign(diff(dist.density$y)))==-2)+1
            branch.dist = max(branch$epath[[1]]$weight)
            if(length(idx) > 0 & any(dist.density$x[idx] > branch.dist))
               width = dist.density$x[idx][dist.density$x[idx] > branch.dist][1]
            else
               width = branch.dist
         
            width = width.scale.factor * width

            cat(paste0("Using branch width: ", format(width, digits=3), " (width.scale.factor: ", width.scale.factor, ")\n"))
         } else
            width = absolute.width
      }
		# Add branch to the backbone tree distance matrix:
		b.tree.dists[cbind(branch.v[-length(branch.v)],branch.v[-1])] = dists[cbind(branch.v[-length(branch.v)],branch.v[-1])]
	
      # connect branch to rest of backbone:
      sprout.from = backbone[which.min(dists[remaining.top.v, backbone])]
      b.tree.dists[sprout.from, remaining.top.v] = dists[sprout.from, remaining.top.v]

		# Add branch vertices to the backbone vertex set
		backbone = c(backbone, branch.v)
				
		tot.remaining = sum(apply(dists[, backbone] > width, 1, all))
	}

   cat(paste0("Outliers: ", tot.remaining, "\nTotal number of branches: ", tot.branches, " (forks: ", (tot.branches-1), ")\n"))
   
	vertebrae.to = V(min.span.tree)[-backbone]
	# no vertebrae to root node (should check if artificial)
	if(rooting.method == "average.start.group")
		vertebrae.from = backbone[-1][apply(dists[backbone[-1],-backbone], 2, which.min)]
	else
		vertebrae.from = backbone[apply(dists[backbone,-backbone], 2, which.min)]
		
	b.tree.dists[cbind(vertebrae.from, vertebrae.to)] = dists[cbind(vertebrae.from, vertebrae.to)]
	
	b.tree = graph.adjacency(b.tree.dists, mode="directed", weighted=TRUE)
	V(b.tree)$is.backbone = FALSE
	V(b.tree)[backbone]$is.backbone = TRUE
	V(b.tree)$is.root = FALSE
	V(b.tree)$is.root[root.v] = TRUE
	V(b.tree)$name = 1:vcount(b.tree)
	
	# Recursively merge backbone:
	
	## Use mean+sd backbone-vertebrae distance
   #    if(absolute.width == 0) {
   #    backbone.vert.dists = E(b.tree)[backbone %--% -backbone]$weight
   #    width = width.scale.factor * (mean(backbone.vert.dists)+sd(backbone.vert.dists))
   # }
   #    else
   #       width = absolute.width
   
	cat(paste0("Backbone fork merge (width: ", format(width, digits=3), "): ", length(backbone)))
	
	b.tree = .recur.merge.backbone(b.tree, dists, width)
	backbone = V(b.tree)[V(b.tree)$is.backbone]
	cat(paste(" -> ", length(backbone), '\n'))	
	
	if(merge.sequential.backbone) {
		seq.merge.width = width/4
		cat(paste0("Backbone sequential vertices merge (width: ", format(seq.merge.width, digits=3), "): ", length(backbone)))
		b.tree = .recur.shorten.backbone(b.tree, dists, seq.merge.width)
		backbone = V(b.tree)[V(b.tree)$is.backbone]
		cat(paste(" -> ", length(backbone), '\n'))	
	}
	
	E(b.tree)$arrow.mode = 0
	
   b.tree = .assign.tree.data(b.tree, rooting.method, grouping, topic.distribs, cell.names)

	cat("Ranking all cells...\n")
   b.tree$ordered.branches = .recur.ordered.branches(b.tree, dists)
   
	options(expressions=save.option.expressions)
	
	return (b.tree)
}

#' Gene Expression Heatmap
#'
#' Plots a heatmap of gene expression, with cells ordered according to the structure computed by \code{\link{compute.backbone.tree}}.
#'
#' @param b.tree igraph object returned by \code{\link{compute.backbone.tree}}.
#' @inheritParams compute.lda
#' @param reorder.genes Boolean (optional). Whether the gene rows should be reordered using a dendrogram of their mean value. 
#' @return \code{data} object reordered according to the backbone tree, such as used to plot the heatmap.
#'
#' @example examples/compute.backbone.tree.example.R
#' @example examples/ct.plot.heatmap.example.R
#'
#' @import gplots
#' @export

ct.plot.heatmap <- function(data, b.tree, log.scale=TRUE, sd.filter=0.7, reorder.genes=TRUE) {
	.library.or.stop('gplots')
   
	data = .normalise.data(data, log.scale, sd.filter)
   
	mybreaks = seq(0, 2*mean(data), length.out=512)
	myheatcol = colorRampPalette(c("green3", "black", "red3"))(length(mybreaks)-1)

   cell.reordering = as.integer(base::unlist(b.tree$ordered.branches))
   
   # gene.reordering = order.genes.by.fit(data[,cell.reordering], V(b.tree)$grouping.label, deg = 4)
         
	dev.new(width=10, height=4)
	if(nrow(data) < 50) lab.row = rownames(data)[gene.reordering] else lab.row = NA
      
	gplots::heatmap.2(data[,cell.reordering], Rowv=reorder.genes, Colv=NA, dendrogram="none", labRow=lab.row, labCol=NA, trace="none", col=myheatcol, breaks=mybreaks, ColSideColors=V(b.tree)$color[cell.reordering], key=FALSE, keysize=0.8, lwid=c(0.8,10), lhei=c(0.1,3.5,0))
   
   return(data[,cell.reordering])
}

#' Plot cell tree with topic distributions
#'
#' Plots a backbone tree (or MST) that was computed with \code{\link{compute.backbone.tree}}, displaying each cell's topic distribution as a pie chart.
#' @param tree An \link[igraph]{igraph} tree, as returned by \code{\link{compute.backbone.tree}} 
#' @param file.output String (optional). Path of a file where the plot should be saved in PDF format (rendered to screen if omitted).
#' @param show.labels Boolean (optional). Whether to write each cell's row number next to its vertex.
#' @param force.recompute.layout Boolean (optional). If set to \code{TRUE}, recomputes the graph's layout coordinates even when present.
#' @param height,width Numeric (optional). Height and width (in inches) of the plot.
#' @param vertebrae.distance Numeric (optional). If non-zero: forces a specific plotting distance (in pixels) between backbone cells and related peripheral cells (`vertebrae').
#' @param backbone.vertex.size,vert.vertex.size Numeric (optional). Diameter (in pixels) of backbone and vertebrae cell vertices.
#' @return An updated \link[igraph]{igraph} object with \code{x} and \code{y} vertex coordinate attributes.
#'
#' @example examples/compute.backbone.tree.example.R
#' @example examples/ct.plot.grouping.example.R
#'
#' @import igraph
#' @export

ct.plot.topics <- function(tree, file.output = NULL, show.labels = FALSE, force.recompute.layout = FALSE, height = 20, width = 10, vertebrae.distance = 0, backbone.vertex.size = 0, vert.vertex.size = 0) {
   return(.plot.b.tree(tree, plot.topics = TRUE, plot.grouping = FALSE, file.output = file.output, show.labels = show.labels, force.recompute.layout = force.recompute.layout, height = height, width = width, vertebrae.distance = vertebrae.distance, backbone.vertex.size = backbone.vertex.size, vert.vertex.size = vert.vertex.size))
}

#' Plot cell tree with grouping information
#'
#' Plots a backbone tree (or MST) that was computed with \code{\link{compute.backbone.tree}}, displaying each cell's grouping.
#' @inheritParams ct.plot.topics
#' @return An updated \link[igraph]{igraph} object with \code{x} and \code{y} vertex coordinate attributes.
#'
#' @example examples/compute.backbone.tree.example.R
#' @example examples/ct.plot.grouping.example.R
#'
#' @import igraph
#' @export

ct.plot.grouping <- function(tree, file.output = NULL, show.labels = FALSE, force.recompute.layout = FALSE, height = 20, width = 10, vertebrae.distance = 0, backbone.vertex.size = 0, vert.vertex.size = 0) {
   return(.plot.b.tree(tree, plot.topics = FALSE, plot.grouping = TRUE, file.output = file.output, show.labels = show.labels, force.recompute.layout = force.recompute.layout, height = height, width = width, vertebrae.distance = vertebrae.distance, backbone.vertex.size = backbone.vertex.size, vert.vertex.size = vert.vertex.size))
}

#' Gene Ontology enrichment analysis
#'
#' Computes enrichment scores for Gene Ontology terms associated with genes in each topic. 
#' @param lda.results A fitted LDA model, as returned by \code{\link{compute.lda}}
#' @param go.db String. Genome-wide annotation with GO mapping for the appropriate organism (e.g. \pkg{org.Mm.eg.db} or \pkg{org.Hs.eg.db}).
#' @param ontology.type (optional). ``BP'' for Biological Process, ``MF'' for Molecular Function, and ``CC'' for Cellular Component.
#' @param reformat.gene.names Boolean. If set to \code{TRUE}, converts all gene names to capitalised lowercase.
#' @param bonferroni.correct Boolean. Unless set to \code{FALSE}, adjust statistical testing p-value threshold for multiple testing.
#' @param p.val.threshold Numeric (optional). P-value significance threshold.
#' @param go.score.class String (optional). Name of the scoring method to use for the Kolmogorov-Smirnov test (e.g. ``weigth01Score'' or ``elimScore''). See \pkg{topGO} documentation for a complete list of scoring methods.
#' @param dag.file.prefix String or \code{FALSE}. If not set to \code{FALSE}, plots individual subgraphs of significant terms for each topic using the string as filename prefix.
#' @return Returns a named list object with ranked tables of significantly enriched GO terms for each topic (`all'), terms that only appear in each topic (`unique') and terms that appear in less than half of the other topics (`rare'). In addition the list object contains an \link[igraph]{igraph} object with the full GO DAG, annotated with each term's p-value and the significance threshold adjusted for multiple testing (Bonferroni method).
#'
#' @example examples/compute.go.enrichment.example.R
#'
#' @import igraph
#' @export

compute.go.enrichment <- function(lda.results, go.db, ontology.type = "BP", reformat.gene.names=FALSE, bonferroni.correct=TRUE, p.val.threshold=if(bonferroni.correct) 0.05 else 0.01,  go.score.class="weight01Score", dag.file.prefix=FALSE) {
   .library.or.stop('topGO', bioconductor = TRUE)
	.library.or.stop('igraph')

	if(! ontology.type %in% c("BP", "MF", "CC"))
		stop("Only 'BP', 'CC' and 'MF' are supported as ontology types")
	
	if(! is.character(go.db))
		go.db = go.db$packageName
	
	if(class(lda.results) != 'topics' & !is.list(lda.results)) {
		.library.or.stop('topicmodels', "to use this type of lda.results object")
		gene.names = lda.results@terms
		k = lda.results@k
		theta = exp(lda.results@beta)
	} else {
		.library.or.stop('maptpx', "to use this type of lda.results object")
		gene.names = rownames(lda.results$theta)
		k = lda.results$K
		theta = t(lda.results$theta)
	}
	
	if(reformat.gene.names) {
		cat("Reformatting gene names (to lower case names starting with capital)")
		.simpleCap <- function(x) {
		    s <- strsplit(x, " ")[[1]]
		    paste(toupper(substring(s, 1, 1)), substring(s, 2),
		          sep = "", collapse = " ")
		}
		
		gene.names = sapply(tolower(gene.names), .simpleCap)
	} 
	
   ks.test = new(go.score.class, testStatistic = topGO::GOKSTest, name = "KS", scoreOrder="decreasing")
	
   results = list()
   tables = list()
   
	for(topic in 1:k) {
		cat(paste("Computing GO enrichment for topic:", topic, '\n'))
		
	   all.genes = theta[topic,]
	   names(all.genes) = gene.names

	   go.data <- new("topGOdata", ontology= ontology.type, allGenes = all.genes, geneSel = function(x) { length(x) }, nodeSize=8, annot=topGO::annFUN.org, mapping=go.db, ID = "symbol" )

   	if(bonferroni.correct)
         cor.p.val.threshold = p.val.threshold/length(go.data@graph@nodes)
      else
         cor.p.val.threshold = p.val.threshold
     
	   ks.results = topGO::getSigGroups(go.data, ks.test)
		
		if(is.character(dag.file.prefix)) {
         # .library.or.stop('Rgraphviz') # Seems topGO should handle that dependency on its own
			topGO::printGraph(go.data, ks.results, firstSigNodes = sum(score(ks.results) < cor.p.val.threshold), fn.prefix = paste(dag.file.prefix, "_T", topic, "_", ontology.type, "_DAG", sep=''), useInfo = "def", pdfSW = TRUE)
		}
      
      results = append(results, list(score(ks.results)))
     
      # showSigOfNodes(go.data, score(ks.results), firstSigNodes = sum(score(ks.results) < cor.p.val.threshold), useInfo = "def", swPlot=FALSE)
      
	   t = topGO::GenTable(go.data, weightKS = ks.results, numChar = 512, topNodes = sum(score(ks.results) < cor.p.val.threshold))
						
      names(t)[3] = 'Total'
      names(t)[6] = 'p-Value'
      tables = append(tables, list(t[c(1,2,3,6)]))
	}

   complete.go.dag = .graph.reverse(igraph.from.graphNEL(topGO::graph(go.data)))
   localGetTermsDefinition = getFromNamespace(".getTermsDefinition", "topGO")
   V(complete.go.dag)$def = localGetTermsDefinition(V(complete.go.dag)$name, ontology(go.data), numChar = 256)
   
   tables.unique = lapply(1:length(tables), function(i) tables[[i]][! tables[[i]]$GO.ID %in%  unique(base::unlist(lapply(tables[-i], function(l) l$GO.ID))),])
   
   tables.rare = lapply(1:length(tables), function(i) {
      t = table(base::unlist(lapply(tables[-i], function(l) l$GO.ID)))
      tables[[i]][! tables[[i]]$GO.ID %in%  names(t[t>(length(tables)/2)]),]
    })
	
	return(list(all=tables, unique=tables.unique, rare=tables.rare, results=results, complete.go.dag=complete.go.dag, adjusted.p.threshold = cor.p.val.threshold))
}

#' Ranking of cells according to backbone tree structure
#'
#' Produces a table of input cells ranked by their position in the backbone tree.
#'
#' @param b.tree An \link[igraph]{igraph} backbone tree, as returned by \code{\link{compute.backbone.tree}}.
#' @param write.to.tex.file Boolean (optional). If not \code{NULL}, outputs LaTeX version of table to file \code{write.to.tex.file}.
#' @return List of all cells, ranked by position in backbone tree, along with topic information.
#'
#' @example examples/compute.backbone.tree.example.R
#' @example examples/cell.ordering.table.example.R
#'
#' @import igraph
#' @import xtable
#' @export

cell.ordering.table <- function(b.tree, write.to.tex.file = NULL) { 
	.library.or.stop('igraph')
   
	if(is.null(b.tree$ordered.branches)) stop("The input tree does not appear to be a backbone tree computed by 'compute.backbone.tree'")
		
   branch.table = .recur.ordered.branch.table(b.tree$ordered.branches);

   branch.table = cbind(branch.table[,1:4], sapply(branch.table[,5], which.max), branch.table[,5])
   colnames(branch.table) = c("branch", "node.label", "cell.name", "cell.group", "main.topic", "topics")
   
   
   if(! is.null(write.to.tex.file)) {
      .library.or.stop('xtable', "to use the 'write.to.tex.file' option")
		
      colors = rainbow(length(V(b.tree)$pie[[1]]))
      
      latex.table = branch.table
      latex.table[,3] = gsub("_", "\\_", branch.table[,3], fixed = TRUE)
      latex.table[,4] = gsub("_", "\\_", as.character(branch.table[,4]), fixed = TRUE)
      latex.table[,6] = sapply(branch.table[,6], function(x) paste0("\\smallbars{", paste(formatC(cumsum(x), digits=2, format='f'), collapse="}{"), "}"))
      
      tex.file = paste0(write.to.tex.file, '.tex')
      sink(tex.file)
      cat('\\documentclass[9pt]{extreport}
      \\usepackage[portrait]{geometry}
      \\usepackage{longtable}
      \\usepackage[usenames,dvipsnames,svgnames,table]{xcolor}
      \\usepackage{tikz}
      \\usepackage{float}
      \\renewcommand\\topfraction{0.85}
      \\renewcommand\\bottomfraction{0.85}
      \\renewcommand\\floatpagefraction{0.85}\n')
      for(i in 1:length(colors))
         cat(paste0('\\definecolor{col', i, '}{HTML}{', substring(colors[i], 2, 7), '}\n'))
      cat('\\newcommand\\smallbars[', length(colors), ']{
      \\begin{tikzpicture}\n')
      
      cat(paste0('\\fill[col1] (0,0) rectangle (#1,0.2);\n \\draw (0,0) rectangle (#1,0.2);\n'))
      
      for(i in 2:length(colors))
         cat(paste0('\\fill[col', i, '] (#', (i-1), ',0) rectangle (#', i, ',0.2);\n \\draw (#', (i-1), ',0) rectangle (#', i, ',0.2);\n'))
      
      cat('\\end{tikzpicture}
      } 
      \\date{}
      \\begin{document}\n')
      
      cat("\\begin{center}\\Large\\textbf{Ordered cells by branch}\\end{center}\\vspace{2em}")
      cat("Legend: ")
      for(i in 1:length(colors))
         cat(paste0('\\colorbox{col', i, '}{Topic \\#', i, '}\n'))
     
      
      tables.by.branch = split(as.data.frame(latex.table[,2:6]), as.character(latex.table[,1]))
      
      for(i in 1:length(tables.by.branch)) {
         t = xtable::xtable(tables.by.branch[[i]], caption=paste0('Branch ', names(tables.by.branch)[i]), row.names=FALSE, display=c('d', 'd', 's', 's', 'd', 's'))
         xtable::align(t) = 'llllll'
         cat(toLatex(t, tabular.environment="longtable", sanitize.text.function=identity, table.placement='H', floating=FALSE, include.rownames=FALSE, caption.placement = 'top'), append=TRUE, sep="\n")
      }
      
      cat('\\end{document}', append=TRUE)
      sink()
   }
   
   return(branch.table)
}

#' LaTeX output for Gene Ontology Enrichment results
#' 
#' Outputs or writes result tables of Gene Ontology enrichment testing in LaTeX format.
#' @param go.results GO Enrichment result list object, such as returned by \code{\link{compute.go.enrichment}}.
#' @param tex.file.name String (optional). If not \code{NULL}: name of the file to save to.
#' @param topic.colors RGB colour vector (optional). Colors to use for each topic.
#' @return GO enrichment results in LaTeX format
#'
#' @example examples/go.results.to.latex.example.R
#'
#' @import xtable
#' @export

go.results.to.latex <- function(go.results, tex.file.name = 'go_terms.tex', topic.colors = rainbow(length(go.results$all))) {
   .library.or.stop('xtable')

	use.temp.file = is.null(tex.file.name)
	
	if(use.temp.file)
		tex.file = tempfile()
   else
      tex.file = tex.file.name
   
	sink(tex.file)
	
   cat('\\documentclass[9pt]{extreport}
   \\usepackage[portrait]{geometry}
   \\usepackage{tabularx}
   \\usepackage[usenames,dvipsnames,svgnames,table]{xcolor}
   \\usepackage{float}
   \\renewcommand\\topfraction{0.85}
   \\renewcommand\\bottomfraction{0.85}
   \\renewcommand\\floatpagefraction{0.85}
   \\date{}
   \\begin{document}\n')
   for(i in 1:length(topic.colors))
      cat(paste0('\\definecolor{col', i, '}{HTML}{', substring(topic.colors[i], 2, 7), '}\n'))
   for(i in 1:length(go.results$all)) {
      t = xtable::xtable(go.results$all[[i]], caption=paste0("\\textcolor{col", i, "}{Topic ", i, "} (All terms)"))
      xtable::align(t) = 'rlXrl'
      cat(toLatex(t, tabular.environment="tabularx", width="\\textwidth", table.placement='H'), append=TRUE, sep="\n")
   }
   cat('\\clearpage\n', append=TRUE)
   
   for(i in 1:length(go.results$rare)) {
      t = xtable::xtable(go.results$rare[[i]], caption=paste0("\\textcolor{col", i, "}{Topic ", i, "} (Terms that appear in less than half of other topics)"))
      xtable::align(t) = 'rlXrl'
      cat(toLatex(t, tabular.environment="tabularx", width="\\textwidth", table.placement='H'), append=TRUE, sep="\n")
   }
   cat('\\clearpage\n', append=TRUE)
   
   for(i in 1:length(go.results$unique)) {
      t = xtable::xtable(go.results$unique[[i]], caption=paste0("\\textcolor{col", i, "}{Topic ", i, "} (Terms that only appear in this topic)"))
      xtable::align(t) = 'rlXrl'
      cat(toLatex(t, tabular.environment="tabularx", width="\\textwidth", table.placement='H'), append=TRUE, sep="\n")
   }
   cat('\\end{document}', append=TRUE)

	sink()
	
	output = readChar(tex.file, file.info(tex.file)$size)
	
   if(use.temp.file) unlink(tex.file)

	return(output)
}


#' Gene Ontology enrichment sets plotting
#'
#' Plots DAG of significantly enriched terms for all topics, along with ancestor nodes.
#'
#' @param go.results GO Enrichment result list object, such as returned by \code{\link{compute.go.enrichment}}.
#' @param up.generations Integer (optional). Number of generations above significant nodes to include in the subgraph.
#' @param only.topics Integer vector (optional). If not \code{NULL}, vector of topics that should be included in the plot (otherwise all topic enrichment sets are used).
#' @param file.output String (optional). If not \code{NULL}, pathname of file to write the plot to.
#' @param p.val.threshold Numeric (optional). P-value treshold to use to select which terms should be plotted.
#' @param only.unique Only display terms that are only significant for one of the topics.
#' @param topic.colors RGB colour vector (optional). Colors to use for each topic.
#' @return An \link[igraph]{igraph} object with the annotated GO DAG.
#'
#' @example examples/ct.plot.go.dag.example.R
#'
#' @import igraph
#' @export

ct.plot.go.dag <- function(go.results, up.generations = 2, only.topics = NULL, file.output=NULL, p.val.threshold = go.results$adjusted.p.threshold, only.unique = FALSE, topic.colors = rainbow(length(go.results$results))) {
   .library.or.stop('igraph')

   if(! is.null(file.output))	pdf(paste(file.output, "pdf", sep='.'), height=10, width=10)

	topic.nums = 1:length(go.results$results)
	
	if(! is.null(only.topics)) {
		topic.nums = topic.nums[only.topics]
		topic.colors = topic.colors[only.topics]
		used.go.results = go.results$results[only.topics]
	} else used.go.results = go.results$results
	
	term.colors = lapply(1:length(used.go.results), function(i) {
	   scores = used.go.results[[i]]
	   scores = scores[scores <= p.val.threshold]
	   rgb = as.list(col2rgb(topic.colors[[i]], alpha=TRUE)[,1]/255)
	   sapply(scores, function(p.val) {
	      score.color = rgb
	      score.color$alpha = 0.3 + min(0.7, 0.1* log(p.val.threshold/p.val))
	      rgb(score.color$red, score.color$green, score.color$blue, score.color$alpha)
	   })
	})

	all.term.colors = base::unlist(term.colors)
	unique.terms = unique(names(all.term.colors))
	sig.term.colors = sapply(unique.terms, function(term) {
	  this.term.colors = all.term.colors[names(all.term.colors) == term]
	  .mixrgb(this.term.colors)
	})

	g = go.results$complete.go.dag
	E(g)$weight = 1
	V(g)[names(used.go.results[[1]])]$p.val = apply(do.call('cbind', used.go.results), 1, min)

   if(only.unique) {
      V(g)$unique = FALSE
      V(g)[names(used.go.results[[1]])]$unique = apply(do.call('cbind', used.go.results), 1, function(p.vals) { (sum(p.vals <= p.val.threshold) == 1) })
      
   	show.vertices = V(g)[V(g)$unique]
   }
   else {
   	show.vertices = V(g)[V(g)$p.val < p.val.threshold]
   }
   
	if(length(show.vertices) == 0) stop(paste0("No significant GO term for p value threshold: ", p.val.threshold, ". Try raising the threshold."))

	if(up.generations > 0) {
		cur.vertices = show.vertices
		V(g)$pred = 0

		for(i in 1:up.generations) {
			for(v in cur.vertices) {
				parents = V(g)[nei(v, mode="in")]
				for(p in parents) {
					if(p %in% show.vertices) next
			
					p = V(g)[p]
					if(p$pred > 0) {
						n = p
						while(n$pred > 0) {
							show.vertices = union(show.vertices, n)
				
							idx = as.integer(V(g)[n]$pred)
							V(g)[n]$pred = 0
							n = V(g)[idx]
						}
			
						n = V(g)[v]
						while(n$pred > 0) {
							show.vertices = union(show.vertices, V(g)[n])
				
							idx = as.integer(V(g)[n]$pred)
							V(g)[n]$pred = 0
							n = V(g)[idx]
						}
					}
					V(g)[p]$pred = v
				}
			}

			cur.vertices = V(g)[V(g)$pred > 0]
		}
	}
	sub.tree = induced.subgraph(go.results$complete.go.dag, show.vertices)
	V(sub.tree)$color = "white"
   
   sig.term.colors = sig.term.colors[intersect(show.vertices$name, names(sig.term.colors))]
	V(sub.tree)[names(sig.term.colors)]$color = sig.term.colors
	term.defs = sapply(V(sub.tree)$def, function(s) paste(strwrap(s, 20), collapse="\n"))
	V(sub.tree)$name = term.defs

	saved.par.mar = par()$mar
	saved.par.lwd = par()$lwd
	
	par(mar=c(0,0,0,0))
	par(lwd = 0.1)

   V(sub.tree)$shape = "circle"
   V(sub.tree)$size = 12
   V(sub.tree)$label.cex = 0.3
   V(sub.tree)$label.color = "black"
   V(sub.tree)$frame.color = NA
   E(sub.tree)$arrow.size = 0.4
	layout = layout.fruchterman.reingold(sub.tree, niter=10000)
   V(sub.tree)$x = layout[,1]
   V(sub.tree)$y = layout[,2]
   
	plot(sub.tree)
	
   legend("topright", legend=paste('Topic:', topic.nums), col=topic.colors, pch=19, cex=0.7)
	
	par(mar= saved.par.mar)
	par(lwd= saved.par.lwd)
   
   if(! is.null(file.output))	dev.off()
	
   return(sub.tree)
}

### Experimental functions: gene ordering using polynomial fit (not quite yet ready):

.polynom.vals <- function(x.vals, max.deg, coeffs) {
	poly.vals = base::unlist(sapply(1:(max.deg+1), function(deg) x.vals^(deg-1) * coeffs[deg]))
	if(is.null(dim(poly.vals)))
		return(sum(poly.vals))
	else
		return(apply(poly.vals, 1, sum))
}

.polynom.deriv.vals <- function(x.vals, max.deg, coeffs) {
	poly.vals = sapply(2:(max.deg+1), function(deg) (deg-1) * x.vals^(deg-2) * coeffs[deg])
	if(is.null(dim(poly.vals)))
		return(sum(poly.vals))
	else
		return(apply(poly.vals, 1, sum))
}

order.genes.by.fit <- function(y, time.grouping, deg = 4, trigger.threshold = 0.6) {
	# y = exprs / apply(exprs, 1, function(x) max(x))
	
	x = sort(sapply(time.grouping, function(x) x+rnorm(1, sd=0.05*(x + 2))))

	polys = apply(y, 1, function(ycol) lm(ycol ~ poly(x, deg, raw=TRUE)))
	min.x = min(x)-2
	max.x = max(x)+5
	genes.threshold.points = sapply(names(polys), function(name) {
		poly = polys[[name]]
		y.min.threshold = (max(poly$fitted)-max(0, min(poly$fitted)))*0.3 +max(0, min(poly$fitted))
		sols = polyroot(poly$coefficients[-1] * 1:deg)
		extrema = Re(sols[abs(Im(sols)) < 1e-10])
		extrema = extrema[extrema > min.x & extrema < max.x]
		maxima = NULL
		maximum = NULL
		next.minimum = NULL
		
		if(length(extrema) != 0) {
			directions = .polynom.deriv.vals(c(extrema+0.1), deg, poly$coefficients)
			maxima = extrema[directions < 0]
			y.vals = .polynom.vals(maxima, deg, poly$coefficients)
			maxima = maxima[y.vals > y.min.threshold]
			y.vals = y.vals[y.vals > y.min.threshold]
			
			if(length(maxima) > 0) {
				minimum = extrema[extrema < min(maxima)]
				minimum = minimum[length(minimum)]
				y.minimum = .polynom.vals(c(minimum), deg, poly$coefficients)
				if(length(extrema[extrema > min(maxima)]))
					next.minimum = extrema[extrema > min(maxima)][1]
			   
            maximum = maxima[1]
         }
		}
      
		if(is.null(maximum)) {
		   if(poly$fitted[1] < min(y.min.threshold, poly$fitted[length(poly$fitted)])) {
   			maximum = max(x)
   			next.minimum = +Inf
         } else {
            maximum = -Inf
            next.minimum = max(x)
         }
		}
      
      
		if(is.null(next.minimum)) {
			if(length(poly$fitted[x > maximum]) == 0)
				next.minimum = +Inf
			else {
				next.minimum = x[x > maximum][which.min(poly$fitted[x > maximum])]
            # if(next.minimum >= poly$fitted[x > maximum][1])
            #    next.minimum = +Inf
			}
		}
      
		c(maximum, next.minimum)
	})

	return(order(genes.threshold.points[1,], genes.threshold.points[2,]))
}

save.per.topic.gene.distribution <- function(lda.results, file.output, top.n = 100) {
   if(class(lda.results) != 'topics' & !is.list(lda.results)) {
   	.library.or.stop('topicmodels', "to use this type of lda.results object")
   	gene.names = lda.results@terms
   	k = lda.results@k
   	theta = exp(lda.results@beta)
   } else {
   	.library.or.stop('maptpx', "to use this type of lda.results object")
   	gene.names = rownames(lda.results$theta)
   	k = lda.results$K
   	theta = t(lda.results$theta)
   }

   
   lapply(1:k, function(this.k) {
      top.n.genes = order(theta[this.k,], decreasing=T)[1:top.n]
      top.n.beta = theta[this.k, top.n.genes]
      names(top.n.beta) = gene.names[top.n.genes]

      pdf(paste0(file.output, "genes_for_topic_", this.k, ".pdf"), width=5, height=12)
      barplot(rev(top.n.beta[1:top.n]), horiz=T, las=2, cex.names=0.5, col=c("blue"), border=NA, main=paste("Probability distribution for Topic", this.k))
      dev.off()
      write.table(top.n.beta, file=paste0(file.output, "genes_for_topic_", this.k, ".csv"), sep="\t", col.names=F)
      top.n.beta
   })
}
