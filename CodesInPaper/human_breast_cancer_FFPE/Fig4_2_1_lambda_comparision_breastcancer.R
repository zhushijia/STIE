deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")

############################################################
# lambda
############################################################
las = sort( unique(c( 0, 10^c(-3:6) )) )
las =  c(2500, 5000, 7500, 6000, 7000 )

results = list()
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)

for(i in 4:length(las))
{
    cat(las[i],'\n')
    results[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=las[i], steps=30, 
                       known_signature=TRUE, known_cell_types=FALSE, min_cells=-1)
}

score = lapply( results, function(x) calculate_BIC(x, ST_expr) )
names(results) = names(score) = las

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
save(results, score, file="BreastCancer_lambda_comparison_2.5xSpot_fullSignature_2500.RData")

########################################################################################################################
####### Visualization
########################################################################################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load( "BreastCancer_lambda_comparison_2.5xSpot_fullSignature.RData")
ind = as.character( sort( unique(c( 0, 10^c(-1:6) )) ) )
results = results[ ind ]
score = score[ind]
PME_diff = get_PME_diff(results) 

################################################################################################
############ draw barplot of RMSE
################################################################################################
pdf("BreastCancer_lambda_comparision_2.5Xspot.pdf")
par(mfrow=c(2,2))
errplot( PME_diff, "PME_diff" )
errplot( lapply(score, function(x) sqrt(x$mse) ), "RMSE" )
errplot( lapply(score, function(x) x$bic ), "BIC" )
dev.off()


