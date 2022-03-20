deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

############################################################
# lambda
############################################################
las = sort( unique(c( 0, 10^c(-3:6) )) )

result = list()
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)

for(i in 1:length(las))
{
    cat(las[i],'\n')
    result[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=las[i], steps=30, 
                       known_signature=TRUE, known_cell_types=FALSE)
}

names(result) = las

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
save(result, file="BreastCancer_lambda_BIC2.RData")


