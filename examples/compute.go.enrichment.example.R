# Load pre-computed LDA model for skeletal myoblast RNA-Seq data from HSMMSingleCell package:
data(HSMM_lda_model)

\donttest{
# Load GO mapping database for 'homo sapiens':
library(org.Hs.eg.db)
# Compute Cellular Component GO enrichment sets for each topic:
go.results = compute.go.enrichment(HSMM_lda_model, org.Hs.eg.db, ontology.type="CC", bonferroni.correct=TRUE, p.val.threshold=0.01)

# Print table of terms that are only significantly enriched in each topic: 
print(go.results$unique)
}
