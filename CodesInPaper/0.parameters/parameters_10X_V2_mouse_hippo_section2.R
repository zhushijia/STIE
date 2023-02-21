library(Matrix)
library(Seurat)
library(hdf5r)
library(data.table)
library(magick)
library(EBImage)
library(quadprog)
library(STIE)

image_path = '/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/V2_MouseBrainCoronalSection2_FFPE/CytAssist_FFPE_Mouse_Brain_Rep2_tissue_image.tif'
morph_feature_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/V2_MouseBrainCoronalSection2_FFPE/CytAssist_FFPE_Mouse_Brain_Rep2_tissue_image_split_images"
selected_spots_path = '/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/V2_MouseBrainCoronalSection2_FFPE/Section2_Hippo2.csv'

spaceranger_outs_dir <- "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/V2_MouseBrainCoronalSection2_FFPE/CytAssist_FFPE_Mouse_Brain_Rep2/outs"
spotfile_path <- paste0(spaceranger_outs_dir, "/spatial/tissue_positions.csv")
matrix_path <- paste0(spaceranger_outs_dir, "/filtered_feature_bc_matrix")
scalefactor_path <- paste0( spaceranger_outs_dir, "/spatial/scalefactors_json.json")
############################################################
########## read image
############################################################
scales <- rjson::fromJSON(file = scalefactor_path)
x_scale = scales$tissue_lowres_scalef
fct <- (55/2)/100 
############################################################
########## read image
############################################################
im <- image_read(image_path)
im <- im[1]
############################################################
########## split image
############################################################
if(0) {
    split_image(image=image_path, split_dir=morph_feature_dir, 
                w=3000, h=3000, margin=100, x_scale=1 )
    run_imageJ_plugin(
        imageJ = "/archive/SCCC/Hoshida_lab/s184554/Code/github/ImageJ/Fiji.app/ImageJ-linux64",
        plugin_macro = NULL,
        split_image_dir = morph_feature_dir, 
        pattern="tif$" )
}
############################################################
########## load features
############################################################
cat("Loading morpholgy features...\n")
cell_info_file = paste0(morph_feature_dir,"/cell_info.RData")
if( file.exists( cell_info_file ) ) {
    x = load(cell_info_file)
} else {
    cell_info <- merge_feature( morph_feature_dir )
}
morphology_fts = cell_info$cell_feature
contour = cell_info$cell_contour
morphology_fts$X = sapply( cell_info$cell_contour, function(c) mean(c$x) )
morphology_fts$Y = sapply( cell_info$cell_contour, function(c) mean(c$y) )
morphology_fts$pixel_x = morphology_fts$X * x_scale
morphology_fts$pixel_y = morphology_fts$Y * x_scale
contour = lapply(contour,function(x) data.frame(pixel_x=x[,1]*x_scale,pixel_y=x[,2]*x_scale) )
###########################################################
features_list <- list(size = c("Area", "Major", "Minor", "Width", "Height", "Feret", "Perim."), 
                      shape = c("Round", "Circ."),
                      angle = c('FeretAngle','Angle'))
PCs <- do.call(cbind, lapply(features_list, function(f) {
  X = morphology_fts[, f]
  X2 = scale(X)
  prcomp(X2)$x[, 1]
}))
morphology_fts <- cbind(morphology_fts, PCs)
###########################################################
##########  Spot coordinates
############################################################
cat("Loading spot coordinates...\n")
spot_coordinates <- read.csv(file = spotfile_path, header = TRUE,
                             col.names=c("barcode","tissue","row","col","imagerow","imagecol"))
spot_coordinates <- data.frame( spot_coordinates, 
                                pixel_x=spot_coordinates$imagecol*x_scale,
                                pixel_y=spot_coordinates$imagerow*x_scale )
spot_radius <- calculate_spot_radius(spot_coordinates, fct)
selected_spots <- read.table(selected_spots_path, sep=",", header=T)
spot_coordinates = subset(spot_coordinates, barcode %in% as.character(selected_spots$Barcode) )
############################################################
########## load 10X transcriptome
############################################################
cat("load 10X transcriptome...\n")
count <- as.data.frame(t(Read10X(matrix_path)))
spot_coordinates <- subset(spot_coordinates, barcode%in%rownames(count) )
ST_expr = count[match( as.character(spot_coordinates$barcode), rownames(count) ),]
############################################################
## known signature, unkonwn cell types, konwn cell segmentation
############################################################
if( deconvolution )
{
    known_signature = TRUE
    known_cell_types = FALSE
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
    Signature = read.delim("Major_cell_types_marker_genes.txt", header=T, row.names=1)
    Signature = Signature[,!colnames(Signature)%in%c("Ependymal.cells")]
    Signature = as.matrix(Signature)[rownames(Signature)%in%colnames(ST_expr),]
    colnames(Signature) = gsub("Pyramidal.neurons.|.interneurons|.like.cells|Granule.cells.", "", colnames(Signature))
    Signature = Signature[,order(colnames(Signature))]
}
