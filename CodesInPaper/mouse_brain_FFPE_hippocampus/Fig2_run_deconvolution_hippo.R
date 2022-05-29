deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo3.R")

############################################################
## run deconvolution
############################################################

cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)

result <- STIE(ST_expr, Signature, cells_on_spot, features, 
                 lambda=0, steps=30, 
                 known_signature=known_signature, known_cell_types=known_cell_types)

cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
colors = c("magenta", "blue", "green", "black", "orange", "cyan")



##########################################################################################
#### whole image
##########################################################################################
pdf("hippo_STIE_whole.pdf")
plot_sub_image(im=im, 
               x_scale=args$x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T, 
               axis_tick=2000, axis_col='grey'  )

dev.off()


##########################################################################################
#### selecte region
##########################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")

pdf("hippo_STIE_xoff8000_yoff10500.pdf")
#png("tmp.png", w=1000,h=1000,res=300)
plot_sub_image(im=im, w=5400, h=4100, xoff=8200, yoff=11200, 
               x_scale=0.2, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=0.2, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=0.2, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=T, plot_cell=F  )

plot_sub_image(im=im, w=1050, h=1500, xoff=11100, yoff=11500, 
               x_scale=0.2, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

plot_sub_image(im=im, image_transparency=1, w=1050, h=1500, xoff=11100, yoff=11500, 
               x_scale=0.2, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

dev.off()


##########################################################################################
#### in and out spot
##########################################################################################
cells_on_spot_in = get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, spot_radius)
cells_on_spot_out = cells_on_spot[ ! as.character(cells_on_spot$cell_id) %in% as.character(cells_on_spot_in$cell_id) , ]
cell_types_in = cell_types[ match( as.character(cells_on_spot_in$cell_id), names(cell_types) ) ]
cell_types_out = cell_types[ match( as.character(cells_on_spot_out$cell_id), names(cell_types) ) ]
contour_in = cell_info$cell_contour[ match(names(cell_types_in), names(cell_info$cell_contour)) ]
contour_out = cell_info$cell_contour[ match(names(cell_types_out), names(cell_info$cell_contour)) ]

a = length(unique(names(cell_types_in)))
b = length(unique(names(cell_types_out)))
c = a+b
a/c; b/c

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")
colors = c("magenta", "blue", "green", "black", "orange", "cyan")

pdf("hippo_STIE_xoff8000_yoff10500_spot_in_out.pdf")

plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=0.2, spot_coordinates=spot_coordinates, 
               contour=contour_in, cell_types=cell_types_in, color_use=colors, plot_spot=T, plot_cell=T  )

plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=0.2, spot_coordinates=spot_coordinates, 
               contour=contour_out, cell_types=cell_types_out, color_use=colors, plot_spot=T, plot_cell=T  )

dev.off()
##########################################################################################
######## morphology features
##########################################################################################
library(vioplot)

pdf("hippo_STIE_morphology_features.pdf", w=6, h=4)

colors = c("magenta", "blue", "green", "black", "orange", "cyan")
fs = c("Area", "Solidity", "Circ.", "Skew", "Kurt", "Round", "Eccentricity","IntDen","Feret")
par( mfrow=c(3,3) )
par(mar = c(2, 2, 2, 2))
lapply( fs, function(f) {
    x = split( cells_on_spot[,f], result$cell_types )
    vioplot(cells_on_spot[,f]~result$cell_types, col=colors, main=f, xlab=f, rectCol="darkgrey")
})

dev.off()

##########################################################################################
######## clustering
##########################################################################################



##########################################################################################
######## spatial CCI
##########################################################################################



