source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/DWLS_on_Wu_etal_2021_BRCA_Signature")
x = load("Visium_FFPE_Human_Breast_Cancer_DWLS_on_Wu_etal_2021_BRCA_Signature.RData" )

############################################################
# load cell deep learning cell segmentation
############################################################
thres = 0.5
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
ratio = seq(0.5,8,0.5)
results = list()
for(i in 1:length(ratio))
{
    cat(ratio[i],'\n')
    cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, ratio[i]*spot_radius)
    results[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                     known_signature=TRUE, known_cell_types=FALSE)
}

score = lapply( results, function(x) calculate_BIC(x, ST_expr) )
names(results) = names(score) = ratio

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
save(results, score, file=paste0("BreastCancer_spot_DL_CellSegmentation_thres",thres,".RData") )


################################################################################################
########## visualization
################################################################################################
thres = 0
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load( paste0("BreastCancer_spot_DL_CellSegmentation_thres",thres,".RData") )

ind = as.character(seq(0.5,5.5,0.5))
results = results[ind]
score = score[ind]
PME_diff = get_PME_diff(results) 


################################################################################################
############ draw barplot of RMSE
################################################################################################
pdf( paste0("BreastCancer_spot_DL_CellSegmentation_thres",thres,"_spot_comparision.pdf") )
colors = c( "steelblue", "darkred", "black", "cyan", "yellow", "darkorange", "#4DAF4A",'red','green')
celltypes = c("Bcells", "Tcells", "Plasmablasts", "Myeloid",         
              "Endothelial", "CAFs", "PVL", 
              "NormalEpithelial", "CancerEpithelial" )
par(mfrow=c(2,2))
errplot( PME_diff, "PME_diff" )
errplot( lapply(score, function(x) sqrt(x$mse) ), "RMSE")
errplot( lapply(score, function(x) x$bic ), "BIC")
barplot1(score, colors, celltypes)
dev.off()



################################################################################################
########## leaf plot
################################################################################################
prop_on_spot = lapply(score, function(x) {
    y = x$celltypes_on_spot
    t(apply(y,1,function(z) z/sum(z)))
} )

rmse = sapply(score,function(x) mean(sqrt(x$mse)))
ri = 5
ref = prop_on_spot[[ri]]

dist = lapply(prop_on_spot, function(x) {
    spots = intersect( rownames(x), rownames(ref) )
    x2 = x[match(spots,rownames(x)), ]
    ref2 = ref[match(spots, rownames(ref)), ]
    sapply( 1:length(spots), function(i) dist( rbind(x2[i,], ref2[i,]) ) )
    #sapply( 1:length(spots), function(i) cor(x2[i,], ref2[i,])^2 )
}  )

#pdf('BreastCancer_spot_new_get_cells_on_spot_full_signature_leaf_plot.pdf',w=6,h=6)
title = paste0(names(dist),"x")
leaf.plot(dist, upper=1:(ri-1), lower=(ri+1):length(dist), axis=TRUE, title=title, ref_title=title[ri])
leaf.plot(dist, upper=1:(ri-1), lower=(ri+1):length(dist), axis=TRUE)
leaf.plot(dist, upper=1:(ri-1), lower=(ri+1):length(dist), axis=FALSE)
#dev.off()

