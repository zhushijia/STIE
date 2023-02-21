deconvolution = TRUE
source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/parameters/parameters_10X_Visium_mouse_brain_hippo3.R")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)
############################################################
## run deconvolution
############################################################
result = STIE(ST_expr, Signature, cells_on_spot, features=c('shape','size'),
              steps=30, known_signature=TRUE, known_cell_types=FALSE, 
              lambda=0 )
############################################################
## plot 
############################################################
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
colors = c("magenta", "blue", "green", "black", "orange", "cyan")
plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=0.2, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

############################################################
library(vioplot)
colors = c("magenta", "blue", "green", "black", "orange", "cyan")
features_list <- list( size = c("Area", 'Major', 'Minor', 'Width', 'Height', 'Feret','Perim.'),
                       shape = c("Round", 'Circ.'), 
                       angle = c('FeretAngle','Angle') )
fs = do.call(c,features_list[1:2])
par( mfrow=c(3,3) )
par(mar = c(2, 2, 2, 2))
lapply( fs, function(f) {
    x = split( cells_on_spot[,f], result$cell_types )
    vioplot(cells_on_spot[,f]~result$cell_types, col=colors, main=f, xlab=f, rectCol="darkgrey")
})
