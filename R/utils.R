# Loading required library with an error message if failure
.library.or.stop <- function(lib, msg = "to use this function", bioconductor=FALSE) {
   if(! requireNamespace(lib)) {
		if(bioconductor)
			stop(paste0("Bioconductor package '", lib, "' must be installed ", msg, ". Please install it manually and try again."))
		else
			stop(paste0("Package '", lib, "' must be installed ", msg, ". Please install it manually and try again."))
	}
}

# Using additive colour theory to produce combination of RGB colours:
.mixrgb <- function(colors) {
   rgb.colors = lapply(colors, function(col) as.list(col2rgb(col, alpha=TRUE)[,1]/255))
   sum.t = sum(sapply(rgb.colors, function(col) col$alpha))
   t = sapply(rgb.colors, function(col) col$alpha/sum.t)
   
   alpha = max(sapply(rgb.colors, function(col) col$alpha))
   
   mix.col = lapply(1:3, function(col.idx) {
      sqrt(sum(sapply(1:length(t), function(i) t[i]*rgb.colors[[i]][[col.idx]]^2)))
   })
   
   return(rgb(mix.col[1], mix.col[2], mix.col[3], alpha))
}

# computes chi-square distance
.chi.square.dist <- function(x, y) { sqrt(sum((x-y)^2/(x+y), na.rm=TRUE)) }

#DEBUG:
.sink.reset <- function(){
  for(i in seq_len(sink.number())){
	  cat("resetting sink\n")
      sink(NULL)
  }
}

# reverse direction of all edges on a directed graph 
.graph.reverse <- function (graph) {
  if (!is.directed(graph))
    return(graph)
  e <- get.data.frame(graph, what="edges")
  ## swap "from" & "to"
  neworder <- 1:length(e)
  neworder[1:2] <- c(2,1)
  e <- e[neworder]
  names(e) <- names(e)[neworder]
  graph.data.frame(e, vertices = get.data.frame(graph, what="vertices"))
}
