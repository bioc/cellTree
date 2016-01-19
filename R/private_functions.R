# log-scale and/or filter-out data with low standard-deviation
.normalise.data <- function(data, log.scale=TRUE, sd.filter=0.5) {
	if(log.scale) {
		cat('Converting to log values...\n')
		data = log(data+1)
	}
	if(sd.filter != FALSE & !is.null(sd.filter)) {
		keep.rows = (apply(data, 1, sd) > sd.filter)
		cat(paste0('Filtering out rows with standard deviation < ', sd.filter, ' (', nrow(data), ' -> ', sum(keep.rows), ')...\n'))
		data = data[keep.rows,]
	}
	return(as.matrix(data))
}

# recursively computes coordinates for a backbone-tree layout
.recur.tree.layout <- function(t, cur.coord = c(which(V(t)$is.root), 0, 0), done.vertices = numeric(), dir = 0) {
	if(max(igraph::degree(t)) > 5)
		stop("recursive layout function can only handle up to degree 5");
	
	done.vertices = c(done.vertices, cur.coord[1])
	children = difference(V(t)[nei(cur.coord[1])], done.vertices)
	
	n.children = length(children)
	if(n.children == 0) {
		all.coords = cur.coord
	} else	if(n.children == 1) {
		all.coords = rbind(cur.coord, .recur.tree.layout(t, c(as.integer(children[1]), cur.coord[2]+dir, cur.coord[3]+1+0.5*(1-abs(dir))), done.vertices, dir))
	} else if(n.children == 2) { # TODO: make a cleaner function that can handle arbitrary number of children
		all.coords = rbind(cur.coord,
			.recur.tree.layout(t, c(as.integer(children[1]), cur.coord[2]-1, cur.coord[3]+1), done.vertices, -1),
			.recur.tree.layout(t, c(as.integer(children[2]), cur.coord[2]+1, cur.coord[3]+1), done.vertices, 1))
	} else if(n.children == 3) {
		all.coords = rbind(cur.coord,
			.recur.tree.layout(t, c(as.integer(children[1]), cur.coord[2]-1, cur.coord[3]+1), done.vertices, -1),
			.recur.tree.layout(t, c(as.integer(children[2]), cur.coord[2], cur.coord[3]+1.5), done.vertices, 0),
			.recur.tree.layout(t, c(as.integer(children[3]), cur.coord[2]+1, cur.coord[3]+1), done.vertices, 1))
	} else if(n.children == 4) {
		all.coords = rbind(cur.coord,
			.recur.tree.layout(t, c(as.integer(children[1]), cur.coord[2]-1, cur.coord[3]+1), done.vertices, -1),
			.recur.tree.layout(t, c(as.integer(children[2]), cur.coord[2]-0.33, cur.coord[3]+1.25), done.vertices, -0.33),
			.recur.tree.layout(t, c(as.integer(children[3]), cur.coord[2]+0.33, cur.coord[3]+1.25), done.vertices, 0.33),
			.recur.tree.layout(t, c(as.integer(children[4]), cur.coord[2]+1, cur.coord[3]+1), done.vertices, 1))
	}
	return(all.coords)
}

# recursively merge backbone vertices that are close-enough together
.recur.merge.backbone <- function(b.tree, dists, width, cur = V(b.tree)[V(b.tree)$is.root], grand.parent = 0) {

	children = V(b.tree)[nei(as.character(cur), mode='out')]
	# children = children[children != grand.parent]

   if(length(children) > 0)
      children = children[children$is.backbone]
	
	if(length(children) > 1) {
		
		grand.children = V(b.tree)[nei(as.character(children), mode='out')]
      
		max.dists = apply(dists[children, grand.children, drop=FALSE], 1, max)
		best.child = which.min(max.dists)
		
		max.dist.to.children = max(dists[cur, children])
		
		if(max.dists[best.child] <= max(max.dist.to.children, width)) {

			V(b.tree)[as.character(children[-best.child])]$is.backbone = FALSE

			new.grand.children = V(b.tree)[nei(as.character(children[-best.child]), mode='out')] 
			
			b.tree = b.tree - E(b.tree)[inc(as.character(children[-best.child]))]
			
			b.tree = add.edges(b.tree, c(rbind(rep(as.character(children[best.child]), length(new.grand.children)), as.character(new.grand.children))), attr=list(weight=dists[children[best.child], new.grand.children]))
			b.tree = add.edges(b.tree, c(rbind(rep(as.character(children[best.child]), length(children)-1), as.character(children[-best.child]))), attr=list(weight=dists[children[best.child], children[-best.child]]))
			
			children = children[best.child]
		}

		for(child in children) {
			b.tree = .recur.merge.backbone(b.tree, dists, width, child)
		}
	} else if(length(children) == 1) {
		b.tree = .recur.merge.backbone(b.tree, dists, width, children[1])
	}
	
	return(b.tree)
}

# applies some transformations to make the grouping vector usable
.format.grouping <- function(grouping, colors = NULL) {
	if(is.null(grouping))
		return(NULL)
	
   if(is.numeric(grouping)) {
		grouping.labels = sort(unique(grouping))
		grouping.norm = sapply(grouping, function(g) which(grouping.labels == g))
		grouping.colours = colorRampPalette(c("yellow", "red"))(max(grouping.norm))
   } else {
      if(is.character(grouping)) {
         grouping.labels = unique(grouping)
			grouping.norm = match(grouping, grouping.labels)
		} else {
         grouping.factors = droplevels(as.factor(grouping))
   	   grouping.norm = as.integer(grouping.factors)
         grouping.labels = levels(grouping.factors)
      }
      if(is.null(colors))
         grouping.colours = rainbow(n=max(grouping.norm), start=0.51, end=0.40)
      else
         grouping.colours = colors[1:max(grouping.norm)]
         
	}
	
	return(list(labels = grouping.labels, norm = grouping.norm, colours = grouping.colours))
}

# heuristics + user preferences to find the start group
.find.start.group <- function(grouping, lda.results, start.group.label = NULL) {
	if(is.null(grouping))
		return(NULL)
	
	if(! is.null(start.group.label))
		start.group = which(grouping$labels == start.group.label)
   else {
		if(is.numeric(grouping$labels))
			start.group = which(grouping$labels == min(grouping$labels))
		else {
         dists = get.cell.dists(lda.results)

         # assuming that the start group is the one with the lowest intra-group mean square distance...
         min.means = sapply(1:max(grouping$norm), function(group) {
            if(sum(grouping$norm == group) <= 1)
               1000
            else
               min(apply(dists[grouping$norm == group, grouping$norm == group]^2, 1, mean))
         }) 
         start.group = which.min(min.means)
		}
	}
	
   cat(paste0("Using start group: ", grouping$labels[start.group], " (", start.group, ")\n"))
   
	return(start.group)
}

# heuristics + user preferences to pick the rooting method
.pick.rooting.method <- function(rooting.method, start.group) {
	if(is.null(rooting.method)) {
		if(! is.null(start.group))
         rooting.method = "center.start.group"
      else
         rooting.method = "longest.path" 
	} else {
		rooting.methods = c("longest.path", "center.start.group", "average.start.group")
		idx = pmatch(rooting.method, rooting.methods)
		if(is.na(idx)) stop(paste0("Unknown rooting method. Must be one of: ", paste0(rooting.methods, collapse=", ")))
		rooting.method = rooting.methods[idx]
	}
   
   if((rooting.method == 'average.start.group' | rooting.method == 'center.start.group') & is.null(start.group))
      stop(paste0("Cannot use the '", rooting.method, "' rooting method without providing a grouping. Please use 'longest.path' instead or provide a grouping for the cells."))
   
   cat(paste0("Using rooting method: ", rooting.method, "\n"))
   
   return(rooting.method)
}

# adds metadata as igraph tree attributes
.assign.tree.data <- function(tree, rooting.method, grouping, topic.distribs, cell.names) {
   
	# Following is only used if graphing without pies
   if(! is.null(grouping)) {   	
      if(rooting.method == "average.start.group") {
   		grouping$norm = c(max(grouping$norm)+1, grouping$norm)
   		grouping$colours = c(grouping$colours, 'black')
         grouping$labels = c(grouping$labels, 'root')
   	}

   	V(tree)$color = grouping$colours[grouping$norm]
      V(tree)$group.idx = grouping$norm
      V(tree)$grouping.label = grouping$labels[grouping$norm]
      tree$grouping.colors = grouping$colours
      tree$grouping.labels = grouping$labels
	}
   
	shapes = rep("pie", vcount(tree))
	if(rooting.method == "average.start.group")
		shapes[V(tree)$is.root] = "circle"
	
	V(tree)$pie = as.list(data.frame(t(topic.distribs)))
	V(tree)$name = 1:vcount(tree)
	# V(tree)$name[! V(tree)$is.backbone] = ""
	V(tree)$shape = shapes
   if(length(cell.names) == vcount(tree))
      V(tree)$cell.name = cell.names
   else 
      V(tree)$cell.name = c("Root", cell.names)
 
   return(tree)
}

# recursively lists the cells and backbone tree branches from the root
.recur.ordered.branches <- function(b.tree, dists, cur = V(b.tree)[V(b.tree)$is.root]) {   
   cur = V(b.tree)[cur] # making sure we are using tree vertices
   v = V(b.tree)[nei(cur, mode='out')]

   next.v = v[v$is.backbone]
   vert.v = v[! v$is.backbone]
   
   cur.plus.vert = cur
   
   if(length(vert.v) > 0) {
      parent = V(b.tree)[nei(cur, mode='in')]
      
		if(length(parent) == 0) {
	      if(length(next.v) == 0)
				vert.dists = dists[cur, vert.v]
			else
	         vert.dists = - apply(dists[next.v, vert.v, drop=FALSE]/(dists[next.v, vert.v]+dists[cur, vert.v]), 2, min)
		}
		else {
	      before.dists = dists[vert.v, parent]/(dists[vert.v, parent]+dists[vert.v, cur])
	      if(length(next.v) > 0) {
	         after.dists = - apply(dists[next.v, vert.v, drop=FALSE]/(dists[next.v, vert.v]+dists[cur, vert.v]), 2, min)
            				
	         vert.dists = apply(rbind(before.dists, after.dists), 2, function(x) x[which.min(abs(x))])
	      } else vert.dists = before.dists
		}
      cur.plus.vert = c(cur, vert.v)[order(c(0, vert.dists), decreasing=TRUE)]
   } else cur.plus.vert = cur
      
   if(length(next.v) > 1) {
      ret = sapply(next.v, function(x) .recur.ordered.branches(b.tree, dists, x))
      ret = list(cur.plus.vert, ret)
   } else if(length(next.v) == 1) {
      ret = .recur.ordered.branches(b.tree, dists, next.v[[1]])
      ret[[1]] = c(cur.plus.vert, ret[[1]])
   }
   else ret = list(cur.plus.vert)
   
   return(ret)
}

.compute.tree.layout <- function(tree, ratio, vertebrae.distance = 0, backbone.vertex.size = 0, vert.vertex.size = 0) {
	cat("Computing tree layout...\n")
	
	if(length(V(tree)$is.backbone) == 0 | is.null(tree$ordered.branches)) {
		# Not a b-tree. Use standard tree layout:
		coords = layout_as_tree(tree, root=V(tree)[V(tree)$is.root])
		coords = layout.norm(coords, -1, 1, -ratio, ratio)
		
	   if(backbone.vertex.size == 0)   backbone.vertex.size = 0.01
		vertices.sizes = rep(backbone.vertex.size, vcount(tree))
	}
	else {
		### Give some room to recursive functions:
		save.option.expressions = getOption("expressions")
		options(expressions=max(1000,vcount(tree)*10))
	
		### Compute btree layout coordinates:
		coords = matrix(0, vcount(tree), 2)
	
		backbone.t = induced.subgraph(tree, V(tree)$is.backbone)

		if(max(igraph::degree(backbone.t)) >= 5) {
	      cat("Degree of the graph is too high: layout out as regular tree.")

			v.coords = layout_as_tree(backbone.t, root=V(backbone.t)[V(backbone.t)$name == V(tree)[V(tree)$is.root]$name])
	      
         coords[V(tree)$is.backbone,] = v.coords
		}
		else {
		   # special case: path graph
		   if(length(tree$ordered.branches) <= 1) {
		      v.coords = .recur.tree.layout(backbone.t, dir=1)
		   } else {
				v.coords = .recur.tree.layout(backbone.t)
		   }
		
		   coords[V(tree)$is.backbone,] = v.coords[order(v.coords[,1]),c(2,3)]
		   coords[,2] = -coords[,2]
		}
	
	   coords = layout.norm(coords)
		coords = layout.norm(coords, -1, 1, -ratio, ratio)
   
	   vertebrae.dist.factor = 0.1
   
	   backbone = V(tree)[V(tree)$is.backbone]
   
	   min.backbone.dist = 5
		for(v in backbone) {
			neighbors = V(tree)[nei(v, mode='out')]
	      neighbors = neighbors[neighbors$is.backbone]
      
	      if(length(neighbors) > 0) {
	   		this.min = rowSums((coords[v,] - coords[neighbors,,drop=FALSE])^2)
	         min.backbone.dist = min(min.backbone.dist, this.min)
	      }
	   }
   
	   min.backbone.dist = sqrt(min.backbone.dist)
   
	   if(backbone.vertex.size == 0)   backbone.vertex.size = min.backbone.dist * 0.2
	   if(vert.vertex.size == 0)   vert.vertex.size = min.backbone.dist * 0.08
	   if(vertebrae.distance == 0)   vertebrae.distance = min.backbone.dist * 0.45

		for(v in backbone) {
			nei = difference(V(tree)[nei(v)], backbone)
			if(length(nei) > 0) {
				coords[nei,] = vertebrae.distance * layout_as_star(induced.subgraph(tree, union(v, nei)))[-1,]

				coords[nei,1] = coords[nei,1] + coords[v,1]
				coords[nei,2] = coords[nei,2] + coords[v,2]
			}
		}
		
		options(expressions=save.option.expressions)
		
		vertices.sizes = rep(vert.vertex.size, vcount(tree))
		vertices.sizes[V(tree)$is.backbone] = backbone.vertex.size
	}
	
   V(tree)$x = coords[,1]
   V(tree)$y = coords[,2]
   tree$ratio = ratio
   
	V(tree)$size = vertices.sizes
   
   
   return(tree)
}

# generic function to plot backbone tree with options (topics or grouping) called by public functions ct.plot.grouping and ct.plot.topics
.plot.b.tree <- function(tree, plot.topics, plot.grouping, file.output=NULL, show.labels=FALSE, force.recompute.layout=FALSE, height=20, width=10, vertebrae.distance = 0, backbone.vertex.size = 0, vert.vertex.size = 0) {
	.library.or.stop('igraph')

   if(! is.null(file.output))	pdf(paste(file.output, "pdf", sep='.'), height=height, width=width)
		
	#Sanity check:
	if(class(tree) != 'igraph') stop('input tree needs to be an igraph object')
		
   save.par.lwd = par('lwd')
   save.par.mar = par('mar')
	par(lwd = 0.1)
	par(mar=c(0.1,0.1,0.1,0.1))
	
	# Reorder so that backbone vertices are drawn on top:
	if(! is.null(V(tree)$is.backbone))
		tree = permute.vertices(tree, order(order(V(tree)$is.backbone)))
	
	if(is.null(V(tree)$is.root)) {
		cat("Tree is unrooted: using vertex 1 as root.\n")
		V(tree)$is.root = FALSE
		V(tree)[1]$is.root = TRUE
	}
   # root.v = which(V(tree)$is.root)

	ratio = par("pin")[2]/par("pin")[1]
   if(is.null(tree$ratio)) tree$ratio = 0
      
	if(is.null(V(tree)$x) | tree$ratio != ratio |  force.recompute.layout) {
	   tree = .compute.tree.layout(tree, ratio, vertebrae.distance = vertebrae.distance, backbone.vertex.size = backbone.vertex.size, vert.vertex.size = vert.vertex.size)
	}
	
	if(show.labels & !is.null(V(tree)$name)) labels = V(tree)$name else labels = NA
	
   	# layout=layout.reingold.tilford(tree, root=root.v),
   if(!is.null(V(tree)$x) & V(tree)$x[which.max(V(tree)$y)] <= 0) legend.position = "topright" else legend.position = "topleft"
   
   if(plot.grouping == TRUE) {
   	## Grouping:
      
      shapes = rep("circle", vcount(tree))
   	plot.igraph(tree, edge.arrow.mode = 0, vertex.label.cex = 0.4, edge.width = 0.5, rescale = FALSE, margin = 0,  vertex.shape=shapes, vertex.label=labels, vertex.shape=shapes, vertex.size= V(tree)$size * dev.size('px')[1] / 2 )
      if(! is.null(tree$grouping.labels))
         legend(legend.position, legend=paste('Group:', tree$grouping.labels), col=tree$grouping.colors, pch=19, cex=0.7)
      
   } else {
		plot.igraph(tree, edge.arrow.mode = 0, vertex.label.cex = 0.4, edge.width = 0.5, rescale = FALSE,  vertex.pie.color=list(rainbow(length(V(tree)$pie[[1]]))), vertex.label=labels, vertex.size= V(tree)$size * dev.size('px')[1] / 2)
   
      legend(legend.position, legend=paste('Topic', 1:length(V(tree)$pie[[1]])), col=rainbow(length(V(tree)$pie[[1]])), pch=19, cex=0.7)
      
   }


	if(! is.null(file.output))	dev.off()
	
	par(lwd = save.par.lwd)
	par(mar = save.par.mar)

	return(tree)
}

.recur.ordered.branch.table <- function(branch, cur.branch.name=c(1)) {
	if(!is.list(branch)) {
		return(cbind(paste(cur.branch.name, collapse='.'), branch$name, branch$cell.name, branch$grouping.label, branch$pie))
	}
	
	ret = NULL
	if(!is.list(branch[[1]])) {
		ret = rbind(ret, .recur.ordered.branch.table(branch[[1]], cur.branch.name=cur.branch.name))
		branch = branch[-1]
	}
   
	if(length(branch) > 0) {
		for(i in seq_along(branch[[1]])) {
			sub.branch = branch[[1]][[i]]
			ret = rbind(ret, .recur.ordered.branch.table(sub.branch, cur.branch.name=c(cur.branch.name, i)))
		}
	}
	
	return(ret)
}


# merge two backbone nodes together (and re-connects vertebrae appropriately)
.merge.backbone.node.to <- function(b.tree, dists, stay.backbone, become.vertebra, stay.backbone.children = NULL, become.vertebra.children = NULL) {	
	
	if(is.null(stay.backbone.children)) {
		stay.backbone.children = V(b.tree)[nei(as.character(stay.backbone))]
		stay.backbone.children = stay.backbone.children[! stay.backbone.children$is.backbone]
	}
	if(is.null(become.vertebra.children)) {
		become.vertebra.children = V(b.tree)[nei(as.character(become.vertebra))]
		become.vertebra.children = become.vertebra.children[! become.vertebra.children$is.backbone]
	}
			
	V(b.tree)[as.character(become.vertebra)]$is.backbone = FALSE 
	
	if(length(become.vertebra.children) > 0) {
		b.tree = b.tree - E(b.tree)[inc(as.character(become.vertebra.children))]
		b.tree = add.edges(b.tree, c(rbind(rep(as.character(stay.backbone), length(become.vertebra.children)), as.character(become.vertebra.children))), attr=list(weight=dists[stay.backbone$name, become.vertebra.children$name]))
	}
		
	to.become.vert = V(b.tree)[nei(as.character(become.vertebra), mode='in')]
	to.become.vert = to.become.vert[to.become.vert != stay.backbone]
			
	if(length(to.become.vert) > 0)
		b.tree = add.edges(b.tree,  c(as.character(to.become.vert), rbind(rep(as.character(stay.backbone), length(to.become.vert)))), attr=list(weight=dists[stay.backbone$name, to.become.vert$name]))

	from.become.vert = V(b.tree)[nei(as.character(become.vertebra), mode='out')]
	from.become.vert = from.become.vert[from.become.vert != stay.backbone]
	if(length(from.become.vert) > 0)
		b.tree = add.edges(b.tree,  c(rbind(rep(as.character(stay.backbone), length(from.become.vert)), as.character(from.become.vert))), attr=list(weight=dists[stay.backbone$name, from.become.vert$name]))
	
	b.tree = b.tree - E(b.tree)[inc(as.character(become.vertebra))]
		
	b.tree = add.edges(b.tree, c(as.character(stay.backbone), as.character(become.vertebra)), attr=list(weight=dists[stay.backbone$name, become.vertebra$name]))
			
	if(! is.null(E(b.tree)$arrow.mode))
		E(b.tree)$arrow.mode = 0
	
	return(b.tree)
}

# recursively merge sequential backbone nodes that are close enough
.recur.shorten.backbone <- function(b.tree, dists, tree.width, cur = V(b.tree)[V(b.tree)$is.root]) {
	backbone = V(b.tree)[V(b.tree)$is.backbone]
	vertebrae = V(b.tree)[!V(b.tree)$is.backbone]
	
	cur.vertebrae = V(b.tree)[nei(as.character(cur))]
	cur.vertebrae = cur.vertebrae[cur.vertebrae %in% vertebrae]
	
	edges = E(b.tree)[as.character(cur) %->% as.character(backbone)]
		
	if(length(edges) == 0) return(b.tree)

	stay.backbone = NULL

	for(i in 1:length(edges)) {
		
		edge = edges[i]
		if(edge$weight > tree.width) next
		
		succ = V(b.tree)[to(edge)]
		succ.vertebrae = V(b.tree)[nei(as.character(succ))]
		succ.vertebrae = succ.vertebrae[succ.vertebrae %in% vertebrae]
		
		if (! any(dists[cur$name, succ.vertebrae$name] > tree.width)) {
			if (! any(dists[succ$name, cur.vertebrae$name] > tree.width)) {
				if(sum(dists[cur$name, union(succ.vertebrae$name, cur.vertebrae$name)]^2) < sum(dists[succ$name, union(succ.vertebrae$name, cur.vertebrae$name)]^2)) {
					become.vertebra = succ
					stay.backbone = cur
				} else {
					become.vertebra = cur
					stay.backbone = succ
				}
			} else {
				become.vertebra = succ
				stay.backbone = cur
			}
		} else if (! any(dists[succ$name, cur.vertebrae$name] > tree.width)) {
			become.vertebra = succ
			stay.backbone = cur
		} else {
			next
		}
					
		b.tree = .merge.backbone.node.to(b.tree, dists, stay.backbone, become.vertebra)

		break
	}
	
	if(is.null(stay.backbone)) {
		successors = V(b.tree)[nei(as.character(cur), mode='out')]
		successors = successors[successors %in% backbone]
					
		for(i in 1:length(successors))
			b.tree = .recur.shorten.backbone(b.tree, dists, tree.width,  V(b.tree)[as.character(successors[i])])
	} else {
			
		b.tree = .recur.shorten.backbone(b.tree, dists, tree.width, V(b.tree)[as.character(stay.backbone)])
	}
	
	return(b.tree)
}
