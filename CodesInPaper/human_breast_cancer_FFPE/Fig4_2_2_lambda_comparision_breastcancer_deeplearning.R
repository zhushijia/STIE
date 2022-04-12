source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/DWLS_on_Wu_etal_2021_BRCA_Signature")
x = load("Visium_FFPE_Human_Breast_Cancer_DWLS_on_Wu_etal_2021_BRCA_Signature.RData" )

############################################################
# load cell deep learning cell segmentation
############################################################
thres = 0
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/")
cells = data.frame(fread("Breast_Cancer_FFPE_information_model40_image20.csv"))
morphology_fts = with(cells, data.frame(
    cell_id = paste0("cell_",rownames(cells)),
    cell_types=nucleus_class,
    probability=probability,
    X=centroid_x0,
    Y=centroid_y0,
    Area=area,
    Major=major_axis_length,
    Minor=minor_axis_length,
    Eccentricity=eccentricity,
    Orientation=orientation,
    Solidity=solidity) )
#morphology_fts = subset(morphology_fts, ! cell_types %in% c("blood") )
morphology_fts = morphology_fts[order(morphology_fts$probability,decreasing=T),]
morphology_fts = subset(morphology_fts,probability>thres)
morphology_fts$pixel_x = morphology_fts$X * args$x_scale
morphology_fts$pixel_y = morphology_fts$Y * args$x_scale
features = c("Area", "Eccentricity")

############################################################
# lambda
############################################################
las = sort( unique(c( 0, 10^c(-1:6) )) )

results = list()
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)

for(i in 1:length(las))
{
    cat(las[i],'\n')
    results[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=las[i], steps=30, 
                        known_signature=TRUE, known_cell_types=FALSE, min_cells=-1)
}

score = lapply( results, function(x) calculate_BIC(x, ST_expr) )
names(results) = names(score) = las

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
save(results, score, file=paste0("BreastCancer_lambda_comparison_2.5xSpot_fullSignature_DL_thres",thres,".RData") )


########################################################################################################################
####### Visualization
########################################################################################################################
thres = 0
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load( paste0("BreastCancer_lambda_comparison_2.5xSpot_fullSignature_DL_thres",thres,".RData") )
PME_diff = get_PME_diff(results) 

################################################################################################
############ draw barplot of RMSE
################################################################################################
pdf( paste0("BreastCancer_lambda_comparison_2.5xSpot_fullSignature_DL_thres",thres,".pdf") )
par(mfrow=c(2,2))
errplot( PME_diff, "PME_diff" )
errplot( lapply(score, function(x) sqrt(x$mse) ), "RMSE" )
errplot( lapply(score, function(x) x$bic ), "BIC" )
dev.off()


