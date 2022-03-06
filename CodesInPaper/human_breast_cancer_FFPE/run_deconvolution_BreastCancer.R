deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

STIE.dir = system.file(package = "STIE")
load( paste0(STIE.dir,"/data/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R") )

############################################################
## run deconvolution
############################################################

cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)

result <- STIE(ST_expr, Signature, cells_on_spot, features, 
               lambda=0, steps=30, morphology_steps=ceiling(steps/3),
               known_signature=known_signature, known_cell_types=known_cell_types)

cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]


##########################################################################################
######## morphology features
##########################################################################################
library(vioplot)
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")

pdf("BreastCancer_morphology_features.pdf")

levels = c("Bcells", "Tcells", "Plasmablasts","CancerEpithelial", "NormalEpithelial", "CAFs", "Endothelial" )
levels_col = c( "steelblue", "#4DAF4A", "darkorange", "black", "yellow", "darkred", "cyan")
fs = c("Area", "Round", "Solidity","IntDen","Feret", "Circ.", "Skew", "Kurt", "Eccentricity")
par( mfrow=c(3,3) )
lapply( fs, function(f) {
    x = split( cells_on_spot[,f], result$cell_types )
    group = factor( result$cell_types, levels=levels)
    vioplot(cells_on_spot[,f]~group, col=levels_col, las=3,xlab=NULL,ylab=f)
})

dev.off()

#### selectec region

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")

pdf("BreastCancer_result.pdf")
colors = c( "steelblue", "darkred", "black", "cyan", "yellow", "darkorange", "#4DAF4A")
plot_sub_image(im=im, w=2000, h=2000, xoff=10000, yoff=10500, 
               x_scale=1, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

plot_sub_image(im=im, w=3000, h=3000, xoff=10000, yoff=10000, 
               x_scale=1, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

#### whole image
plot_sub_image(im=im, 
               x_scale=args$x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, 
               plot_spot=F, plot_cell=T, 
               axis_tick=0, axis_col='grey'  )

dev.off()



