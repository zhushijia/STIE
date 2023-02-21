deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/parameters/parameters_10X_Visium_humanBreastCancer.R")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)
outputDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE_search"

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
  # i = 6
  cluster = read.csv(paste0(spaceranger_outs_dir, "/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
  cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
  all( as.character(cluster[,1]) %in% rownames(ST_expr) )
  ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
  Signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
  
  pdf( paste0("HumanBreastCancer_STIE_searchPaths_clustering_",i ,".pdf") )
  paths = STIE_search(ST_expr, Signature, cells_on_spot, 
                      steps=30, known_signature=FALSE, known_cell_types=FALSE, 
                      lambdas=lambdas,
                      criterion = criterion )
  dev.off()
  save(paths, file=paste0("HumanBreastCancer_STIE_searchPaths_clustering_",i ,".RData") )
}

############################################################
## plot 
############################################################
setwd(outputDir)
i=6
load(paste0("HumanBreastCancer_STIE_searchPaths_clustering_",i ,".RData"))
result = paths[[3]]$best_result[[2]]
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]

myCol = c( "red", "blue", "green", "black", "cyan", "yellow", "purple", "steelblue",
           "darkorange", "darkred", "grey", "blue4", "chartreuse4", "burlywood1", "darkgoldenrod4" )
myCol = c( "black", "steelblue", "darkorange", "purple", "yellow", "darkred", "purple", "steelblue",
           "darkorange", "darkred", "grey", "blue4", "chartreuse4", "burlywood1", "darkgoldenrod4" )

colors = myCol[1:ncol(result$Signature)]

plot_sub_image(im=im, 
               x_scale=0.01, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

plot_sub_image(im=im, w=2000, h=2000, xoff=10000, yoff=10500, 
               x_scale=0.5, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )
