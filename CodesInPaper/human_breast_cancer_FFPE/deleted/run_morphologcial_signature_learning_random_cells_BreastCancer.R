deconvolution = FALSE
clustering = FALSE
signature_learning = TRUE

known_signature = FALSE
known_cell_types = TRUE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

############################################################
## deconvolution
############################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("BreastCancer_spot_BIC_new_get_cells_on_spot.RData")

STIE_result = result

cells_on_spot = STIE_result[['2.8']]$cells_on_spot
cell_types = STIE_result[['2.8']]$cell_types
lambda=0
steps=20
known_signature=FALSE
known_cell_types=TRUE
Signature=NULL
min_cells=2
morphology_steps=10
Signature = NULL

############################################################
## clustering
############################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("HumanBreastCancer_clustering.RData")

STIE_result = result

cells_on_spot = STIE_result[[6]]$cells_on_spot
cell_types = STIE_result[[6]]$cell_types
lambda=0
steps=30
known_signature=FALSE
known_cell_types=TRUE
Signature=NULL
min_cells=2
morphology_steps=10
Signature = NULL

################################################################################
####### randomly assign cell types
################################################################################
num = 100

set.seed(1234)
selected_cell_ids = do.call(c, lapply( unique(cell_types), function(ct) {
    cell_id_ct = unique(names(cell_types)[cell_types == ct])
    n2 = min(num, length(cell_id_ct))
    sample(cell_id_ct, n2)
} ))
cells_on_spot$cell_types = NA
index = which( as.character(cells_on_spot$cell_id) %in% selected_cell_ids )
cells_on_spot$cell_types[ index ] = cell_types[ index ]
################################################################################

result = STIE(ST_expr, Signature=Signature, cells_on_spot, features, lambda=lambda, steps=steps, 
              morphology_steps=morphology_steps,
              known_signature=FALSE, known_cell_types=TRUE)

table(STIE_result[['2.8']]$cell_types, result$cell_types)
table(STIE_result[[6]]$cell_types, result$cell_types)


