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
args$image = '/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Visium_FFPE_Mouse_Brain_image.jpg'
args$feature_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/split_images2"

args$spotfile = '/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/Visium_FFPE_Mouse_Brain/outs/spatial/tissue_positions_list.csv'

args$spaceranger_cout_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/Visium_FFPE_Mouse_Brain"
scalefactor_paths <- paste( args$spaceranger_cout_dir, "outs/spatial/scalefactors_json.json", sep="/")

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
im <- image_read(args$image, args$croparea, args$x_scale)

############################################################
########## split image
############################################################
if(0) {
    split_image(image=args$image, split_dir=args$feature_dir, 
                w=3000, h=3000, margin=100, x_scale=1 )
    run_imageJ_plugin(
        imageJ = "/archive/SCCC/Hoshida_lab/s184554/Code/github/ImageJ/Fiji.app/ImageJ-linux64",
        plugin_macro = NULL,
        split_image_dir = args$feature_dir, 
        pattern="jpg$" )
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

STdata = load_STdata( "Visium_FFPE_Mouse_Brain", args$spaceranger_cout_dir, umi_cutoff=0, is_normalized=FALSE ) 
bcs <- STdata$bcs_merge
count = STdata$matrix[[1]]

spot_coordinates <- subset(spot_coordinates, barcode%in%as.character(bcs$barcode) )
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
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
    Signature = read.delim("Major_cell_types_marker_genes.txt", header=T, row.names=1)
    Signature = Signature[,!colnames(Signature)%in%c("Ependymal.cells")]
    Signature = as.matrix(Signature)[rownames(Signature)%in%colnames(ST_expr),]
    colnames(Signature) = gsub("Pyramidal.neurons.|.interneurons|.like.cells|Granule.cells.", "", colnames(Signature))
    Signature = Signature[,order(colnames(Signature))]
}

############################################################
## UNKNOWN signature, UNKNOWN cell types, known cell segmentation
############################################################
if( clustering )
{
    known_signature = FALSE
    known_cell_types = FALSE
    cluster = read.csv( paste0(args$spaceranger_cout_dir, "/outs/analysis/clustering/graphclust/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    tmp_expr = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    Signature = t(apply(tmp_expr, 2, function(x) tapply(x,cluster[,2],mean) ))
}



############################################################
## run STIE
############################################################

cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)

result <- STIE(ST_expr, Signature, cells_on_spot, features, 
               lambda=0, steps=30, morphology_steps=ceiling(steps/3),
               known_signature=known_signature, known_cell_types=known_cell_types)

cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]

#### selecte region
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
myCol = c( "red", "blue", "green", "black", "cyan", "yellow", "purple", "steelblue",
           "darkorange", "darkred", "grey", "blue4", "chartreuse4", "burlywood1", "darkgoldenrod4" )
colors = myCol[1:ncol(Signature)]

colors = get_my_colors(ncol(Signature),mode=2)

pdf("hippo_STIE_xoff8000_yoff10500.pdf")
plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=1, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

plot_sub_image(im=im, w=1050, h=1500, xoff=11100, yoff=11500, 
               x_scale=1, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=1, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=T, plot_cell=F  )
dev.off()

#### whole image
pdf("hippo_STIE_xoff8000_yoff10500.pdf")
plot_sub_image(im=im, 
               x_scale=args$x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T, 
               axis_tick=2000, axis_col='grey'  )

dev.off()




##########################################################################################
######## morphology features
##########################################################################################
library(vioplot)

myCol = c( "red", "blue", "green", "black", "cyan", "yellow", "purple", "steelblue",
           "darkorange", "darkred", "grey", "blue4", "chartreuse4", "burlywood1", "darkgoldenrod4" )
colors = myCol[1:ncol(Signature)]
colors = get_my_colors(ncol(Signature),mode=2)
fs = c("Area", "Solidity", "Circ.", "Skew", "Kurt", "Round", "Eccentricity")
par( mfrow=c(3,3) )
lapply( fs, function(f) {
    x = split( cells_on_spot[,f], result$cell_types )
    vioplot(cells_on_spot[,f]~result$cell_types, col=colors)
})
vioplot( FC_C_C9 , FC_C_C9[range1] , FC_C_C9[range2] ,FC_C_C9[range3] , col='gold')		

