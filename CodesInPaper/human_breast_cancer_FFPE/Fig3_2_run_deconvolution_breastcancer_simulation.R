deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

############################################################
## run deconvolution
############################################################
if(0) {
    cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)
    
    res_real <- STIE(ST_expr, Signature, cells_on_spot, features, 
                     lambda=1e3, steps=30, 
                     known_signature=known_signature, known_cell_types=known_cell_types, min_cells=-1)
} else {
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
    load( "BreastCancer_lambda_comparison_2.5xSpot_fullSignature.RData")
    res_real = results[['1000']]
    cells_on_spot = res_real$cells_on_spot
}

STdata_sim <- simulate_STdata(cells_coordinates=res_real$cells_on_spot, 
                              cell_types=res_real$cell_types, 
                              Signature=res_real$Signature,
                              spot_coordinates_ref=spot_coordinates, 
                              spot_diameter_pixel_ref=spot_radius*2, 
                              spot_size_ref=55, 
                              spot_size_sim=5, 
                              x_scale=args$x_scale )

res_sim <- STIE(ST_expr, Signature, cells_on_spot, features, 
                lambda=1e3, steps=30, 
                known_signature=known_signature, known_cell_types=known_cell_types, min_cells=-1)


cells = unique(names(res_real$cell_types))
cell_types_real = res_real$cell_types[ match(cells, names(res_real$cell_types)) ]
cell_types_sim = res_sim$cell_types[ match(cells, names(res_sim$cell_types)) ]
table(cell_types_real,cell_types_sim)


############################################################
## visualization
############################################################

cell_types = res_sim$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
spot_coordinates_sim = STdata_sim$spot_coordinates

cc = data.frame(celltypes = c("Plasmablasts","Bcells", "Tcells", "Myeloid", 
                              "Endothelial", "CAFs", "PVL", 
                              "CancerEpithelial"),
                colors = c( "#4DAF4A", "darkorange", "steelblue", "yellow", 
                            "cyan", "darkred", "purple",
                            "black" ) )
cc = subset(cc, celltypes%in%res_sim$cell_types)
colors = as.character(cc$colors)[order(as.character(cc$celltypes))]

#### select region
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")

pdf("BreastCancer_deconvolution_2.5xSpot_lambda_1e3_fullSignature_simulation.pdf")

plot_sub_image(im=im, w=2000, h=2000, xoff=10000, yoff=10500, 
               x_scale=0.5, spot_coordinates=spot_coordinates_sim, 
               plot_spot=T, plot_cell=F, spot_radius=spot_radius/args$x_scale/11,
               spot_cols='darkred')

plot_sub_image(im=im, w=3000, h=3000, xoff=8000, yoff=10000, 
               x_scale=0.5, spot_coordinates=spot_coordinates_sim, 
               plot_spot=T, plot_cell=F, spot_radius=spot_radius/args$x_scale/11,
               spot_cols='darkred')


plot_sub_image(im=im, w=3000, h=3000, xoff=15000, yoff=18000,  
               x_scale=0.5, spot_coordinates=spot_coordinates_sim, 
               plot_spot=T, plot_cell=F, spot_radius=spot_radius/args$x_scale/11,
               spot_cols='darkred')

#### whole image
plot_sub_image(im=im, 
               x_scale=args$x_scale, spot_coordinates=spot_coordinates_sim, 
               contour=contour2, cell_types=cell_types, color_use=colors, 
               plot_spot=F, plot_cell=T, 
               axis_tick=0, axis_col='grey'  )

dev.off()



