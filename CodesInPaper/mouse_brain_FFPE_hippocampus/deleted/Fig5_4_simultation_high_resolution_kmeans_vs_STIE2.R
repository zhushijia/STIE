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

##########################################################################################
k=5
Km = kmeans(ST_expr_sim, k, algorithm="Forgy")
Signature_kmeans = t( apply(ST_expr_sim, 2, function(x) tapply(x,Km$cluster,mean) ))

prop_kmeans = lapply( 1:ncol(Signature_kmeans), function(i) solveNNLS( Signature_sim, Signature_kmeans[,i] ) )
data.frame(sapply(prop_kmeans,which.max),sapply(prop_kmeans,max))
lapply(prop_kmeans,sort,decreasing=T)

stie = STIE(ST_expr_sim, Signature_kmeans, cells_on_spot_sim, features, 
           lambda=0, steps=30, 
           known_signature=FALSE, known_cell_types=FALSE)
Signature_stie = stie$Signature

prop_stie = lapply( 1:ncol(Signature_stie), function(i) solveNNLS( Signature_sim, Signature_stie[,i] ) )
data.frame(sapply(prop_stie,which.max),sapply(prop_stie,max))
lapply(prop_stie,sort,decreasing=T)

##########################################################################################

prop_kmeans = lapply( 1:ncol(Signature_kmeans), function(i) solveNNLS( Signature_sim, Signature_kmeans[,i] ) )
data.frame(sapply(prop_kmeans,which.max),sapply(prop_kmeans,max))
lapply(prop_kmeans,sort,decreasing=T)

prop_stie = lapply( 1:ncol(Signature_stie), function(i) solveNNLS( Signature_sim, Signature_stie[,i] ) )
data.frame(sapply(prop_stie,which.max),sapply(prop_stie,max))
lapply(prop_stie,sort,decreasing=T)


##########################################################################################
real_cluster = with( cells_on_spot_sim, tapply(cell_types, spot, function(x) paste(sort(unique(x)),collapse="_") ))
sim_cluster = Km$cluster
sim_cluster = sim_cluster[ match( names(real_cluster), names(sim_cluster) ) ]
table(real_cluster, sim_cluster)












