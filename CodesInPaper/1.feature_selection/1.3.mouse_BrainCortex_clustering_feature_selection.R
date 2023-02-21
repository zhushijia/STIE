deconvolution = FALSE
clustering = TRUE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/parameters/parameters_10X_Visium_mouse_brain_cortex.R")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)
outputDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE_search"

criterion = "rmse"
outputDir = paste0(outputDir, "_", criterion)
lambdas=c(0,1e1,1e2,1e3,1e4,1e5,1e6)

############################################################
## run clustering
############################################################
dir.create(outputDir)
setwd(outputDir)

for(i in 3:10)
{
  cat(i,'\n')
  # i=6
  cluster = read.csv(paste0(spaceranger_outs_dir, "/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
  cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
  all( as.character(cluster[,1]) %in% rownames(ST_expr) )
  ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
  Signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
  
  pdf( paste0("MouseBrainCortex_STIE_searchPaths_clustering_",i ,".pdf") )
  paths = STIE_search(ST_expr, Signature, cells_on_spot, 
                      steps=30, known_signature=FALSE, known_cell_types=FALSE, 
                      lambdas=lambdas,
                      criterion = criterion )
  dev.off()
  save(paths, file=paste0("MouseBrainCortex_STIE_searchPaths_clustering_",i ,".RData") )
}

############################################################
## plot 
############################################################
result = paths[[1]]$best_result[[2]]
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
colors = c( "red", "steelblue", "blue", "black", "yellow", "cyan")
plot_sub_image(im=im, image_transparency=0, w=14000, h=7000, xoff=6000, yoff=8000, 
               x_scale=x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, 
               plot_spot=F, plot_cell=T  )
