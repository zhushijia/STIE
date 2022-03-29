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
# run on different spot size
############################################################
#ratio = sort( c(0.5,seq(1,3,0.2),seq(1.5,8,1),seq(4,8,1)) )
ratio = seq(0.5,8,0.5)
result = list()
for(i in 1:length(ratio))
{
    cat(ratio[i],'\n')
    cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, ratio[i]*spot_radius)
    result[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                     known_signature=TRUE, known_cell_types=FALSE)
}

score = lapply( result, function(x) calculate_BIC(x, ST_expr) )
names(result) = names(score) = ratio

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
save(result, score, file=paste0("BreastCancer_spot_DL_CellSegmentation_thres",thres,".RData") )

