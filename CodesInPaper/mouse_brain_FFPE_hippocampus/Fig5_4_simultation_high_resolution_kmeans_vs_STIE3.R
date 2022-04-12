deconvolution = TRUE
source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo3.R")

############################################################
## deconvolution on multiorgan
############################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")

if( deconvolution ) {
    
    cat("simulation from deconvolution\n")
    x = load("MouseBrainHippo_spots.RData")
    result = results[["2.5"]]
    
    cells_on_spot = result$cells_on_spot
    cells_on_spot$cell_types = result$cell_types
    spot_id = as.character(cells_on_spot$spot)
    cell_id = as.character(cells_on_spot$cell_id)
    ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id,  ] )
    
    coefs = table( as.character(cells_on_spot$spot), as.character(cells_on_spot$cell_types) )
    coefs = coefs[match( colnames(ST_expr2), rownames(coefs) ), ]
    result$Signature = t( apply( ST_expr2, 1, function(x) solveNNLS( coefs, as.matrix(x), scaled=F ) ) )

} else {
    
    cat("simulation from clustering\n")
    y = load("MouseBrainHippocampus_clustering_2.5Xspot.RData")
    result = results[[5]]
}

STdata_sim <- simulate_STdata(cells_coordinates=result$cells_on_spot, 
                              cell_types=result$cell_types, 
                              Signature=result$Signature,
                              spot_coordinates_ref=spot_coordinates, 
                              spot_diameter_pixel_ref=spot_radius*2, 
                              spot_size_ref=55, 
                              spot_size_sim=5, 
                              x_scale=args$x_scale )

##########################################################################################
Signature_sim = result$Signature
ST_expr_sim = STdata_sim$ST_expr
cells_on_spot_sim = STdata_sim$cells_on_spot
spot_coordinates_sim = STdata_sim$spot_coordinates
spot_coordinates_sim = spot_coordinates_sim[ match( rownames(ST_expr_sim), as.character(spot_coordinates_sim$barcode) ) , ]
##########################################################################################
k=5
set.seed(123)
Km = kmeans(ST_expr_sim, k)
Signature_kmeans = t( apply(ST_expr_sim, 2, function(x) tapply(x,Km$cluster,mean) ))

stie = STIE(ST_expr_sim, Signature_kmeans, cells_on_spot_sim, features, 
           lambda=0, steps=30, 
           known_signature=FALSE, known_cell_types=FALSE)
Signature_stie = stie$Signature

##########################################################################################
kmeans_nnls = sapply( 1:ncol(Signature_kmeans), function(i) solveNNLS( Signature_sim, Signature_kmeans[,i] ) )
stie_nnls = sapply( 1:ncol(Signature_stie), function(i) solveNNLS( Signature_sim, Signature_stie[,i] ) )
kmeans_dwls = dwls(S=Signature_sim, B=Signature_kmeans)
stie_dwls = dwls(S=Signature_sim, B=Signature_stie)

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")
save( ST_expr_sim, cells_on_spot_sim, spot_coordinates_sim,
      Signature_sim, Signature_kmeans, Signature_stie, 
      stie, Km,
      kmeans_nnls, stie_nnls, kmeans_dwls, stie_dwls,
      file="hippo_simulation_highresolution_5nm_from_STIEdeconvolution2.5xspot.RData" )
##########################################################################################


lapply(1:5, function(i) sort( kmeans_nnls[,i], decreasing=T) )
lapply(1:5, function(i) sort( stie_nnls[,i], decreasing=T) )

real_cluster = with( cells_on_spot_sim, tapply(cell_types, spot, function(x) paste(sort(unique(x)),collapse="_") ))
sim_cluster = Km$cluster
sim_cluster = sim_cluster[ match( names(real_cluster), names(sim_cluster) ) ]
table(real_cluster, sim_cluster)


##########################################################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE/clustering_visualization")
pdf( "hippo_simulation_highresolution_5nm_from_STIEdeconvolution2.5xspot.pdf", w=10, h=8 )

myCol= c('magenta','blue','green','black','orange')
myCol2= c('purple','aquamarine','green4','steelblue','yellow')
index = apply(prop_stie,2,which.max)

layout(matrix(1))
par(mar=c(2,2,2,2))
par(mfrow=c(2,3))
barplot(kmeans_nnls[, index], col=myCol, ylim=c(0,1))
barplot(stie_nnls[, index], col=myCol, ylim=c(0,1))
plot.new()
barplot(kmeans_dwls[, index], col=myCol, ylim=c(0,1))
barplot(stie_dwls[, index], col=myCol, ylim=c(0,1))

cell_types = stie$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]

colors = myCol[index]
plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=0.1, spot_coordinates=spot_coordinates_sim, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

colors = myCol2[index]
plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=0.1, spot_coordinates=spot_coordinates_sim, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

spot_cols = myCol[ match( Km$cluster, index ) ]
plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=0.1, spot_coordinates=spot_coordinates_sim, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=T, plot_cell=F,
               spot_cols=spot_cols, fill_spot=T)

spot_cols = myCol2[ match( Km$cluster, index ) ]
plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=0.1, spot_coordinates=spot_coordinates_sim, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=T, plot_cell=F,
               spot_cols=spot_cols, fill_spot=T)


dev.off()







