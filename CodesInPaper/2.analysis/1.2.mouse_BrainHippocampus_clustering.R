source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/parameters/parameters_10X_Visium_mouse_brain_hippo3.R")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)
############################################################
## run clustering
############################################################
k=5
km = read.csv(paste0(spaceranger_outs_dir, "/analysis/clustering/kmeans_",k,"_clusters/clusters.csv") )
km = km[ match(rownames(ST_expr),as.character(km[,1])), ]
cluster = km[,2]; names(cluster)=as.character(km[,1])
Signature_ = t(apply(ST_expr, 2, function(x) tapply(x,cluster,mean) ))
result = STIE(ST_expr, Signature_, cells_on_spot, features=c('shape'),
              steps=30, known_signature=FALSE, known_cell_types=FALSE, 
              lambda=0 )
############################################################
## plot 
############################################################
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
colors = c("orange", "green", "black", "magenta", "blue")
plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=0.2, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )



