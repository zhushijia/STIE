deconvolution = FALSE
clustering = TRUE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_cortex/parameters_cortex.R")

############################################################
## run clustering
############################################################

cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)

results = list()
for(i in 2:10)
{
    cat(i,'\n')
    
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    Signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    
    results[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                        known_signature=FALSE, known_cell_types=FALSE)
}

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE")
save(results, file="MouseBrainCortex_clustering_2.5xSpot_lambda0.RData")


############################################################
# visualization
############################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE")
#load("MouseBrainCortex_clustering.RData")
load("MouseBrainCortex_clustering_2.5xSpot_lambda0.RData")

x_scale = 0.05

for(i in 2:length(results))
{
    
    cell_types = results[[i]]$cell_types
    contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
    
    myCol = c( "red", "blue", "green", "black", "cyan", "yellow", "steelblue",
               "darkorange", "darkred", "grey", "blue4", "purple", "chartreuse4", "burlywood1", "darkgoldenrod4" )
    
    
    colors = myCol[ 1:ncol(results[[i]]$Signature) ]
    colors = colors[ length(colors)+1-rank(table(cell_types)) ]
    
    #colors = get_my_colors(ncol(result[[i]]$Signature),mode=2)
    
    if(i==6) colors = c( "red", "steelblue", "blue", "black", "yellow", "cyan")
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE/clustering_visualization")
    pdf( paste0("MouseBrainCortex_clustering",i,"_scale", x_scale, ".pdf") )
    plot_sub_image(im=im, image_transparency=0, w=14000, h=7000, xoff=6000, yoff=8000, 
                   x_scale=x_scale, spot_coordinates=spot_coordinates, 
                   contour=contour2, cell_types=cell_types, color_use=colors, 
                   plot_spot=F, plot_cell=T  )
    dev.off()
    
}


