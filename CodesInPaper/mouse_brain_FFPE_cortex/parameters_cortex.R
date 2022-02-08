library("argparse")
library("data.table")
library("magick")
library("magrittr")
library("EBImage")
library("ggplot2")
library("dplyr")

library("foreach") 

library("parallel") 
library("doParallel")

library("STIE")

args = list()
args$image = '/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Visium_FFPE_Mouse_Brain_image.jpg'
args$feature_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/split_images2"

args$spotfile = '/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/Visium_FFPE_Mouse_Brain/outs/spatial/tissue_positions_list.csv'

args$spaceranger_cout_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/Visium_FFPE_Mouse_Brain"
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
        plugin_macro = "/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/Code/STIE/v1.5/data/DeepImageJ_plugin_multi_organ_3000_3000.fiji.ijm",
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
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/allen_cortex")
    x = load("AllenCortex_MouseBrain_scRNASeq_DWLS_Signature.RData")
    Signature = Signature[rownames(Signature)%in%colnames(ST_expr),]
    Signature = Signature/100
    ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id, match(rownames(Signature),colnames(ST_expr)) ] )
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
# sub_level_signature
############################################################
if( sub_level_signature )
{
    known_signature = TRUE
    known_cell_types = FALSE
    
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/allen_cortex")
    x = load("AllenCortex_MouseBrain_scRNASeq_DWLS_Signature.RData")
    Signature = Signature[rownames(Signature)%in%colnames(ST_expr),]
    Signature = Signature/100
    
    level1 = c("Endothelial","Endothelial","Endothelial",
               "GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic","GABAergic",
               "Glutamatergic","Glutamatergic","Glutamatergic","Glutamatergic","Glutamatergic","Glutamatergic","Glutamatergic","Glutamatergic","Glutamatergic",
               "Non-Neuronal","Non-Neuronal","Non-Neuronal","Non-Neuronal")
    level2 = c("Endo","Peri","SMC",
               "Lamp5","Meis2","Pvalb","Serpinf1","Sncg","Sst","Vip",
               "CR","L2or3","L4","L5","L5","L6","L6","L6","NP",
               "Astro","Macrophage","Oligo","VLMC")
    level3 = c("Endo","Peri","SMC",
               "Lamp5","Meis2","Pvalb","Serpinf1","Sncg","Sst","Vip",
               "CR","L2or3IT","L4","L5IT","L5PT","L6CT","L6IT","L6b","NP",
               "Astro","Macrophage","Oligo","VLMC")
    Signature_levels = data.frame(level1=level1,level2=level2,level3=level3)
    Signature_levels = Signature_levels[ match( colnames(Signature), as.character(Signature_levels$level3)) , ]
}





cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)




