#### for manipulating images
library(magick)
library(EBImage)

#### for manipulating ST gene expression 
library(Seurat)

#### for the quadratic programming 
library("quadprog")

#### for spatially resolved cell-cell interaction
library(CellChat)
library(NMF)
library(ggalluvial)

library(STIE)


############################################################
image_path="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Visium_FFPE_Mouse_Brain_image.jpg"
im <- image_read(image_path)
############################################################

matrix_paths="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/Visium_FFPE_Mouse_Brain/outs/filtered_feature_bc_matrix.h5"
spotfile <- "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/Visium_FFPE_Mouse_Brain/outs/spatial/tissue_positions_list.csv"

spot_coordinates <- read.csv(file = spotfile, header = FALSE,
                             col.names=c("barcode","tissue","row","col","imagerow","imagecol"))
spot_coordinates <- data.frame(spot_coordinates, pixel_x=spot_coordinates$imagecol, pixel_y=spot_coordinates$imagerow)

ST_expr <- as.data.frame(t(Read10X_h5(matrix_paths)))

spot_coordinates <- subset(spot_coordinates, barcode%in%rownames(ST_expr) )
ST_expr = ST_expr[match( as.character(spot_coordinates$barcode), rownames(ST_expr) ),]

############################################################
feature_dir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/split_images2"
cell_info_file = paste0(feature_dir,"/cell_info.RData")
x = load(cell_info_file)
cell_contour = cell_info$cell_contour
morphology_fts = cell_info$cell_feature
morphology_fts$pixel_x <- morphology_fts$X
morphology_fts$pixel_y <- morphology_fts$Y

############################################################
fct <- 55/2/100 # spot_radius / spot_center_distance
spot_radius <- calculate_spot_radius(spot_coordinates, fct)
features = c("Area", "Round")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.3*spot_radius)
cell_contour = cell_info$cell_contour[ names(cell_info$cell_contour) %in% as.character(cells_on_spot$cell_id) ]
morphology_fts = subset( morphology_fts, cell_id%in%as.character(cells_on_spot$cell_id) )

############################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
Signature = read.delim("Major_cell_types_marker_genes.txt", header=T, row.names=1)
Signature = Signature[,!colnames(Signature)%in%c("Ependymal.cells")]
Signature = as.matrix(Signature)[rownames(Signature)%in%colnames(ST_expr),]
colnames(Signature) = gsub("Pyramidal.neurons.|.interneurons|.like.cells|Granule.cells.", "", colnames(Signature))
Signature = Signature[,order(colnames(Signature))]

############################################################
spaceranger_count_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/Visium_FFPE_Mouse_Brain"
Kmeans = list()
for(i in 2:10)
{
    cat(i,'\n')
    
    cluster = read.csv(paste0(spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    Kmeans[[i]] = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
}

names(Kmeans)[2:10] = 2:10


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")
save( ST_expr, spot_coordinates, cell_contour, morphology_fts, Signature, Kmeans, 
file="MouseBrainHippocampus_10xVisiumFFPE.RData" )


