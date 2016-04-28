# Plot heatmap:
data(HSMM_expr_matrix)
ct.plot.heatmap(HSMM_expr_matrix[1:2000,], b.tree, reorder.genes=FALSE)
