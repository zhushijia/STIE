deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

############################################################
## run deconvolution
############################################################

cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)

result <- STIE(ST_expr, Signature, cells_on_spot, features, 
               lambda=1e3, steps=30, 
               known_signature=known_signature, known_cell_types=known_cell_types, min_cells=-1)

if(0)
{
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
    load( "BreastCancer_lambda_comparison_2.5xSpot_fullSignature.RData")
    result = results[['1000']]
    cells_on_spot = result$cells_on_spot
}
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]


##########################################################################################
##########################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")

library(vioplot)

cc = data.frame(celltypes = c("Plasmablasts","Bcells", "Tcells", "Myeloid", 
                              "Endothelial", "CAFs", "PVL", 
                              "CancerEpithelial"),
                colors = c( "#4DAF4A", "darkorange", "steelblue", "yellow", 
                            "cyan", "darkred", "purple",
                            "black" ) )
cc = subset(cc, celltypes%in%result$cell_types)
colors = as.character(cc$colors)[order(as.character(cc$celltypes))]

pdf("BreastCancer_deconvolution_2.5xSpot_lambda_1e3_fullSignature.pdf")

plot_sub_image(im=im, w=2000, h=2000, xoff=10000, yoff=10500, 
               x_scale=0.5, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

plot_sub_image(im=im, w=3000, h=3000, xoff=8000, yoff=10000, 
               x_scale=0.5, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

plot_sub_image(im=im, w=3000, h=3000, xoff=15000, yoff=18000, 
               x_scale=0.5, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

#### whole image
plot_sub_image(im=im, 
               x_scale=args$x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, 
               plot_spot=F, plot_cell=T, 
               axis_tick=0, axis_col='grey'  )


######## morphology features
levels = as.character(cc$celltypes)
levels_col = as.character(cc$colors)
fs = c("Area", "Round", "Solidity","IntDen","Feret", "Circ.", "Skew", "Eccentricity")
par( mfrow=c(4,2) )
par( mar=c(3,3,3,3) )
lapply( fs, function(f) {
    cat(f,"\n")
    s5 = function(x) substr(x,1,5)
    group = factor( s5(result$cell_types), levels=s5(levels) )
    vioplot(cells_on_spot[,f]~group, col=levels_col, las=3,xlab=NA, main=f, rectCol="grey")
})

dev.off()


