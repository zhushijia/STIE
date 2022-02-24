library(readbitmap)

library(grid)
library(RColorBrewer)
library(cowplot)

library(Seurat)
library(hdf5r)
library(Matrix)

library(data.table)



library("argparse")
library("data.table")
library("magick")
library("magrittr")
library("EBImage")
library("ggplot2")
library("dplyr")
library("quadprog")

library("foreach") 

library("parallel") 
library("doParallel")

library("STIE")

args = list()
args$image = '/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/Visium_FFPE_Mouse_Kidney_image.tif'
args$feature_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/Visium_FFPE_Mouse_Kidney_split_images"

args$spotfile = '/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/count/Visium_FFPE_Mouse_Kidney/outs/spatial/tissue_positions_list.csv'

args$spaceranger_count_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/count/Visium_FFPE_Mouse_Kidney"
scalefactor_paths <- paste( args$spaceranger_count_dir, "outs/spatial/scalefactors_json.json", sep="/")
matrix_paths <- paste( args$spaceranger_count_dir, "outs/filtered_feature_bc_matrix.h5", sep="/")

scales <- rjson::fromJSON(file = scalefactor_paths)
args$x_scale = scales$tissue_lowres_scalef

args$croparea = NULL
args$array_type = "1k"
fct <- ifelse(args$array_type == "1k", 0.25, 50/(sqrt(2)*100))

features <- c("Area", "Solidity", "Circ.", "Skew", "Kurt", "Round", "Eccentricity")
features = c("Area", "Skew","Round","Eccentricity",'IntDen')
features = c("Area", "Round")

############################################################
########## read image
############################################################
im <- image_read(args$image)

############################################################
########## split image
############################################################
if(0) {
    split_image(image=args$image, split_dir=args$feature_dir, 
                w=3000, h=3000, margin=100, x_scale=1 )
    run_imageJ_plugin(
        imageJ = "/archive/SCCC/Hoshida_lab/s184554/Code/github/ImageJ/Fiji.app/ImageJ-linux64",
        plugin_macro = "/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/Code/STIE/v1.5/data/DeepImageJ_plugin_multi_organ_3000_3000.fiji.ijm",
        split_image_dir = args$feature_dir, 
        pattern="tif$" )
}

############################################################
########## load features
############################################################
cat("Loading morpholgy features...\n")

cell_info_file = paste0(args$feature_dir,"/cell_info.RData")
if( file.exists( cell_info_file ) ) {
    c = load(cell_info_file)
} else {
    cell_info <- merge_feature( args$feature_dir )
}
morphology_fts = cell_info$cell_feature
contour = cell_info$cell_contour
morphology_fts$pixel_x = morphology_fts$X * args$x_scale
morphology_fts$pixel_y = morphology_fts$Y * args$x_scale
contour = lapply(contour,function(x) data.frame(pixel_x=x[,1]*args$x_scale,pixel_y=x[,2]*args$x_scale) )

###########################################################
##########  Spot coordinates
############################################################
cat("Loading spot coordinates...\n")
spot_coordinates <- read.csv(file = args$spotfile, header = FALSE,
                             col.names=c("barcode","tissue","row","col","imagerow","imagecol"))
spot_coordinates <- data.frame( spot_coordinates, 
                                pixel_x=spot_coordinates$imagecol*args$x_scale,
                                pixel_y=spot_coordinates$imagerow*args$x_scale )

spot_radius <- calculate_spot_radius(spot_coordinates, fct)


############################################################
########## load 10X transcriptome
############################################################
cat("load 10X transcriptome...\n")

count <- as.data.frame(t(Read10X_h5(matrix_paths)))

spot_coordinates <- subset(spot_coordinates, barcode%in%rownames(count) )
ST_expr = count[match( as.character(spot_coordinates$barcode), rownames(count) ),]

############################################################
######## modes
############################################################

############################################################
## known signature, unkonwn cell types, konwn cell segmentation
############################################################
if( deconvolution )
{
    known_signature = TRUE
    known_cell_types = FALSE
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/Signature/NatComms_2021")
    y = load("NatComms2021_AdultMouseKidney_scRNASeq_DWLS_Signature.RData")
    Signature = Signature[rownames(Signature)%in%colnames(ST_expr),]
    group = data.frame( category = c("Endo", "Podo", 
                                     "LOH", 
                                     "PC", "IC",  
                                     "DCT", "PCT", 
                                     "PST", "EarlyPT", 
                                     "Macro", "Neutro", "Stroma2", "Proliferating"),
                        region = c("Endothelia", "Podocyte", 
                                   "EpitheliaInRenalTubule", 
                                   "EpitheliaInRenalTubule", "EpitheliaInRenalTubule",  
                                   "EpitheliaInRenalTubule", "EpitheliaInRenalTubule", 
                                   "EpitheliaInRenalTubule", "EpitheliaInRenalTubule", 
                                   "Immune", "Immune","Stroma", "Stroma"),
                        subregion = c("Endothelia", "Podocyte", 
                                      "LoopOfHenle", 
                                      "LoopOfHenle", "LoopOfHenle",  
                                      "ConvolutedTubular", "ConvolutedTubular", 
                                      "ProximalStraightTubules", "ProximalStraightTubules", 
                                      "Immune", "Immune","Stroma", "Stroma"))
    
    Signature = Signature[ ,colnames(Signature) %in%c("Endo", "Podo", 
                                                      "PC", "IC",  
                                                      "DCT", "PCT", 
                                                      "PST") ]
    
    NatCommes_fullSignature_res_s2_0 = STIE(ST_expr, Signature, cell_coordinates=morphology_fts, features, spot_coordinates, spot_radius=2*spot_radius, lambda=0,steps=15)
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/count/results/parameter_testing")
    save(NatCommes_fullSignature_res_s2_0, file="NatCommes_fullSignature_res_s2_0.RData")
    
    Signature2 = t(apply(Signature,1,function(x)tapply(x,group$region,mean)))
    NatCommes_Signature_res_s2_0 = STIE(ST_expr, Signature2, cell_coordinates=morphology_fts, features, spot_coordinates, spot_radius=2*spot_radius, lambda=0,steps=15)
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/count/results/parameter_testing")
    save(NatCommes_Signature_res_s2_0, file="NatCommes_Signature_res_s2_0.RData")
    
    Signature3 = t(apply(Signature,1,function(x)tapply(x,group$subregion,mean)))
    NatCommes_subSignature_res_s2_0 = STIE(ST_expr, Signature=Signature3, cell_coordinates=morphology_fts, features, spot_coordinates, spot_radius=2*spot_radius, lambda=0,steps=15)
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/count/results/parameter_testing")
    save(NatCommes_subSignature_res_s2_0, file="NatCommes_subSignature_res_s2_0.RData")
}

############################################################
## UNKNOWN signature, UNKNOWN cell types, known cell segmentation
############################################################
if( clustering )
{
    known_signature = FALSE
    known_cell_types = FALSE
    cluster = read.csv( paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/graphclust/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    tmp_expr = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    Signature = t(apply(tmp_expr, 2, function(x) tapply(x,cluster[,2],mean) ))
}


