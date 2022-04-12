deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo3.R")

############################################################
# lambda
############################################################
las = sort( unique(c( 0, 10^c(-3:6) )) )

results = list()
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)

for(i in 1:length(las))
{
    cat(las[i],'\n')
    results[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=las[i], steps=30, 
                       known_signature=TRUE, known_cell_types=FALSE)
}

score = lapply( results, function(x) calculate_BIC(x, ST_expr) )
names(results) = names(score) = las


outputDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE"
dir.create(outputDir, recursive=T)
setwd(outputDir)
save(results, score, file="MouseBrainHippo_lambda_2.5Xspot.RData")

########################################################################################################################
####### Visualization
########################################################################################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")
load("MouseBrainHippo_lambda_2.5Xspot.RData")
PME_diff = get_PME_diff(results) 

################################################################################################
############ draw barplot of RMSE
################################################################################################
pdf("hippo_lambda_comparision_2.5Xspot.pdf")
par(mfrow=c(2,2))
errplot( PME_diff )
errplot( lapply(score, function(x) sqrt(x$mse) ) )
dev.off()



