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
