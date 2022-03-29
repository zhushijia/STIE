deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo.R")

############################################################
## run deconvolution
############################################################

cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.8*spot_radius)
res_real <- STIE(ST_expr, Signature, cells_on_spot, features, 
               lambda=0, steps=30, 
               known_signature=known_signature, known_cell_types=known_cell_types)

STdata_sim <- simulate_STdata(cells_coordinates=res_real$cells_on_spot, 
                              cell_types=res_real$cell_types, 
                              Signature=res_real$Signature,
                              spot_coordinates_ref=spot_coordinates, 
                              spot_diameter_pixel_ref=spot_radius*2, 
                              spot_size_ref=55, 
                              spot_size_sim=5, 
                              x_scale=args$x_scale )

res_sim <- STIE(ST_expr=STdata_sim$ST_expr, Signature, 
                 cells_on_spot=STdata_sim$cells_on_spot, 
                 features, lambda=0, steps=30, 
                 known_signature=known_signature, known_cell_types=known_cell_types)


cells = unique(names(res_real$cell_types))
cell_types_real = res_real$cell_types[ match(cells, names(res_real$cell_types)) ]
cell_types_sim = res_sim$cell_types[ match(cells, names(res_sim$cell_types)) ]
table(cell_types_real,cell_types_sim)




cell_types = res_sim$cell_types
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

