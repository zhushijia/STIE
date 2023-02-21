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
## run deconvolution
############################################################
dir.create(outputDir)
setwd(outputDir)

pdf("HumanBreastCancer_STIE_searchPaths_deconvolution.pdf")
paths = STIE_search(ST_expr, Signature, cells_on_spot, 
                    steps=30, known_signature=TRUE, known_cell_types=FALSE, 
                    lambdas=lambdas, min_cells=-1,
                    criterion = criterion )
dev.off()

save(paths, file="HumanBreastCancer_STIE_searchPaths_deconvolution.RData")

############################################################
## plot 
############################################################
setwd(outputDir)
load("HumanBreastCancer_STIE_searchPaths_deconvolution.RData")
result = paths[[2]]$best_result[[3]]
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]

cc = data.frame(celltypes = c("Plasmablasts","Bcells", "Tcells", "Myeloid", 
                              "Endothelial", "CAFs", "PVL", 
                              "CancerEpithelial"),
                colors = c( "#4DAF4A", "darkorange", "steelblue", "yellow", 
                            "cyan", "darkred", "purple",
                            "black" ) )
cc = subset(cc, celltypes%in%result$cell_types)
colors = as.character(cc$colors)[order(as.character(cc$celltypes))]

plot_sub_image(im=im, 
               x_scale=0.05, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

plot_sub_image(im=im, w=2000, h=2000, xoff=10000, yoff=10500, 
               x_scale=0.5, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

############################################################
## plot 
############################################################
scores = lapply(paths, function(path) lapply(path$best_result,function(result)get_summary(result,ST_expr))  )
par(mfrow=c(2,2))
plot_pathScore( lapply(scores, function(y) sapply(y,function(x) mean(sqrt(x$mse))) ), "RMSE")
plot_pathScore( lapply(scores, function(y) sapply(y,function(x) mean(x$logLik)) ), "logLik")

plot_pathScore( lapply(scores, function(y) sapply(y,function(x) mean(x$logLik_Expr)) ), "logLik_Expr" )
plot_pathScore( lapply(scores, function(y) sapply(y,function(x) mean(x$logLik_Morp)) ),"logLik_Morp" )

par(mfrow=c(1,1))
plot_pathScore( lapply(scores, function(y) sapply(y,function(x) mean(x$logLik)) ), "logLik")
plot_pathScore( lapply(scores, function(y) sapply(y,function(x) mean(sqrt(x$mse))) ), "RMSE")
