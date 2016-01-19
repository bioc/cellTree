# Load pre-computed LDA model for skeletal myoblast RNA-Seq data from HSMMSingleCell package:
data(HSMM_lda_model)

# Compute cell pairwise distances:
b.tree = get.cell.dists(HSMM_lda_model)
