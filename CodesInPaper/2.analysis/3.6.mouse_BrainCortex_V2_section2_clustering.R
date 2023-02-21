source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/parameters/parameters_10X_V2_mouse_cortex_section2.R")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)
############################################################
## run deconvolution
############################################################
k = 6
pc = prcomp(ST_expr)$x[,1:10]
cluster = kmeans(pc,k)$cluster
Signature_ = t(apply(ST_expr, 2, function(x) tapply(x,cluster,mean) ))

result = STIE(ST_expr, Signature_, cells_on_spot, features=c('shape'),
              steps=30, known_signature=FALSE, known_cell_types=FALSE, 
              min_cells=5, lambda=0 )
############################################################
## plot 
############################################################
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
colors = c("cyan", "blue", "green", "magenta", "black", "orange")
plot_sub_image(im=im, w=13500, h=10200, xoff=6700, yoff=3500, 
               x_scale=x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T)
spot_col = colors[ cluster ]
plot_sub_image(im=im, w=13500, h=10200, xoff=6700, yoff=3500, 
               x_scale=x_scale, spot_coordinates=spot_coordinates, 
               spot_col=spot_col, fill_spot=T,
               plot_spot=T, plot_cell=F  )

