# Load skeletal myoblast RNA-Seq data from HSMMSingleCell package:
library(HSMMSingleCell)
data(HSMM_expr_matrix)

# Run LDA inference using 'maptpx' method for k = 4:
 lda.results = compute.lda(HSMM_expr_matrix, k.topics=4, method="maptpx")

\donttest{
# Run LDA inference using 'maptpx' method for number of topics k = 3 to 6:
 lda.results = compute.lda(HSMM_expr_matrix, k.topics=3:6, method="maptpx")

# Run LDA inference using 'Gibbs' [collapsed sampling] method for number of k = 4 topics:
 lda.results = compute.lda(HSMM_expr_matrix, k.topics=4, method="Gibbs")
}
