# Load pre-computed LDA model for skeletal myoblast RNA-Seq data from HSMMSingleCell package:
data(HSMM_lda_model)

# Load GO mapping database for 'homo sapiens':
library(org.Hs.eg.db)

\donttest{
# Compute GO enrichment sets for each topic:
go.results = compute.go.enrichment(HSMM_lda_model, org.Hs.eg.db, bonferroni.correct=TRUE)

go.dag.subtree = ct.plot.go.dag(go.results, up.generations = 2)

}
