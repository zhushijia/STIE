deconvolution = TRUE
source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/parameters/parameters_10X_V2_mouse_hippo_section2.R")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)

############################################################
## run deconvolution
############################################################
result = STIE(ST_expr, Signature, cells_on_spot, features=c('shape'),
              steps=30, known_signature=TRUE, known_cell_types=FALSE, min_cells=5,
              lambda=0 )

############################################################
## plot 
############################################################
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
colors = c("magenta", "blue", "green", "black", "orange", "cyan")
plot_sub_image(im=im, w=8200, h=4500, xoff=11000, yoff=7300,
               x_scale=0.2, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T)

