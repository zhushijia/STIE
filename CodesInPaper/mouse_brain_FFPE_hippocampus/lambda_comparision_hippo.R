deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo.R")

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

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
save(result, file="MouseBrainHippo_lambda.RData")

score = lapply( result, function(x) calculate_BIC(x,ST_expr) )
rmse = sapply( score, function(x) mean(sqrt(x$mse)) )
plot( rmse )



