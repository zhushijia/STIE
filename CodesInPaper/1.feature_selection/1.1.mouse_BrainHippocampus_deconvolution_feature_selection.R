deconvolution = TRUE
source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/parameters/parameters_10X_Visium_mouse_brain_hippo3.R")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)
outputDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE_search"

criterion = "rmse"
outputDir = paste0(outputDir, "_", criterion)
lambdas=c(0,1e1,1e2,1e3,1e4,1e5,1e6)

############################################################
## run deconvolution
############################################################
dir.create(outputDir)
setwd(outputDir)

pdf("MouseBrainHippocampus_STIE_searchPaths_deconvolution.pdf")
paths = STIE_search(ST_expr, Signature, cells_on_spot, 
                    steps=30, known_signature=TRUE, known_cell_types=FALSE, 
                    lambdas=lambdas,
                    criterion = criterion )
dev.off()

save(paths, file="MouseBrainHippocampus_STIE_searchPaths_deconvolution.RData")



############################################################
## run deconvolution uneqal_prior
############################################################
deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE
source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo3.R")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)
outputDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE_search2_unequal_prior"
dir.create(outputDir)
setwd(outputDir)
pdf("MouseBrainHippocampus_STIE_searchPaths_deconvolution.pdf")
paths = STIE_search(ST_expr, Signature, cells_on_spot, 
                    steps=30, known_signature=TRUE, known_cell_types=FALSE, 
                    lambdas=c(0,1e1,1e2,1e3,1e4,1e5,1e6),
                    criterion = "L2sum" )
dev.off()
save(paths, file="MouseBrainHippocampus_STIE_searchPaths_deconvolution.RData")

############################################################
## plot 
############################################################
setwd(outputDir)
load("MouseBrainHippocampus_STIE_searchPaths_deconvolution.RData")
result = paths[[3]]$best_result[[1]]
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
colors = c("magenta", "blue", "green", "black", "orange", "cyan")
plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=0.2, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

############################################################
## score 
############################################################
setwd(outputDir)
load("MouseBrainHippocampus_STIE_searchPaths_deconvolution.RData")
scores = lapply(paths, function(path) lapply(path$best_result,function(result)get_summary(result,ST_expr))  )
par(mfrow=c(3,3))
lapply(scores, function(s) plot(sapply(s, function(x) mean(sqrt(x$rmse)) )) )
lapply(scores, function(s) plot(sapply(s, function(x) x$logLik )) )
lapply(scores, function(s) boxplot(sapply(s, function(x) x$logLik_Expr ),outline=F,las=3) )
lapply(scores, function(s) boxplot(sapply(s, function(x) x$logLik_Morp ),outline=F,las=3) )

par(mfrow=c(2,2))
plot_pathScore( lapply(scores, function(y) sapply(y,function(x) mean(x$logLik_Expr)) ), "logLik_Expr" )
plot_pathScore( lapply(scores, function(y) sapply(y,function(x) mean(x$logLik_Morp)) ),"logLik_Morp" )
plot_pathScore( lapply(scores, function(y) sapply(y,function(x) mean(x$logLik)) ), "logLik")
plot_pathScore( lapply(scores, function(y) sapply(y,function(x) mean(sqrt(x$mse))) ), "RMSE")




