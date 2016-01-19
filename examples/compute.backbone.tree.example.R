# Load pre-computed LDA model for skeletal myoblast RNA-Seq data from HSMMSingleCell package:
data(HSMM_lda_model)

# Recover sampling time (in days) for each cell:
library(HSMMSingleCell)
data(HSMM_sample_sheet)
days.factor = HSMM_sample_sheet$Hours
days = as.numeric(levels(days.factor))[days.factor]

# Compute near-optimal backbone tree:
b.tree = compute.backbone.tree(HSMM_lda_model, days)
