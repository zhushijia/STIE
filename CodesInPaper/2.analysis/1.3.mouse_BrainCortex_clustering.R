source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/parameters/parameters_10X_Visium_mouse_brain_cortex.R")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)
############################################################
## run clustering
############################################################
k = 6
km = read.csv(paste0(spaceranger_outs_dir, "/analysis/clustering/kmeans_",k,"_clusters/clusters.csv") )
km = km[ match(rownames(ST_expr),as.character(km[,1])), ]
cluster = km[,2]; names(cluster)=as.character(km[,1])
Signature_ = t(apply(ST_expr, 2, function(x) tapply(x,cluster,mean) ))
result = STIE(ST_expr, Signature_, cells_on_spot, features=c('shape'),
              steps=30, known_signature=FALSE, known_cell_types=FALSE, 
              lambda=0 )
############################################################
## STIE plot 
############################################################
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
colors = c( "red", "steelblue", "blue", "black", "yellow", "cyan")
plot_sub_image(im=im, image_transparency=0, w=14000, h=7000, xoff=6000, yoff=8000, 
               x_scale=x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, 
               plot_spot=F, plot_cell=T  )

