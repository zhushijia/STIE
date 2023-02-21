source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/parameters/parameters_10X_V2_mouse_hippo_section1.R")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)
outputDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/V2_MouseBrainCoronalSection1_FFPE/CytAssist_FFPE_Mouse_Brain_Rep1/results/STIE_search"

criterion = "rmse"
outputDir = paste0(outputDir, "_", criterion)
lambdas=c(0,1e1,1e2,1e3,1e4,1e5,1e6)

############################################################
## run deconvolution
############################################################
dir.create(outputDir)
setwd(outputDir)

for(i in 3:10)
{
    cat(i,'\n')
    #i = 5
    pdf( paste0("V2_MouseBrainHippocampus_STIE_searchPaths_clustering_",i ,".pdf") )
    pc = prcomp(ST_expr)$x[,1:10]
    cluster = kmeans(pc,i)$cluster
    Signature = t(apply(ST_expr, 2, function(x) tapply(x,cluster,mean) ))
    
    paths = STIE_search(ST_expr, Signature, cells_on_spot, 
                        steps=30, known_signature=FALSE, known_cell_types=FALSE, 
                        lambdas=lambdas, min_cells=5,
                        criterion = criterion )
    dev.off()
    
    save(paths, file=paste0("V2_MouseBrainHippocampus_STIE_searchPaths_clustering_",i ,".RData") )
}


############################################################
## plot 
############################################################
setwd(outputDir)
i = 5
load(paste0("V2_MouseBrainHippocampus_STIE_searchPaths_clustering_",i ,".RData"))
result = paths[[2]]$best_result[[1]]
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
colors = c("magenta", "blue", "green", "black", "orange", "cyan")
plot_sub_image(im=im, w=9000, h=5000, xoff=15500, yoff=11500,
               x_scale=x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=T, plot_cell=T)



############################################################
########## test
############################################################
set.seed(1234)
i = 5
pc = prcomp(ST_expr)$x[,1:10]
cluster = kmeans(pc,i)$cluster
Signature = t(apply(ST_expr, 2, function(x) tapply(x,cluster,mean) ))

colors = c("magenta", "blue", "green", "black", "orange", "cyan")
colors = c("green", "magenta", "black", "orange", "blue", "green")
spot_cols = colors[cluster]
plot_sub_image(im=im, w=9000, h=5000, xoff=15500, yoff=11500,
               x_scale=x_scale, spot_coordinates=spot_coordinates, 
               spot_cols = spot_cols, plot_spot=T, plot_cell=F, fill_spot=T)


features_list <- list( size = c("Area", 'Major', 'Minor', 'Width', 'Height', 'Feret','Perim.'),
                       shape = c("Round", "Eccentricity", 'Circ.'), 
                       angle = c('FeretAngle','Angle'),
                       solidity = c("Solidity") )
PCs <- do.call(cbind, lapply( features_list, function(f) {
    X = cells_on_spot[,f]
    X2 = scale(X)
    prcomp(X2)$x[,1]
}))
cells_on_spot2 <- cbind(cells_on_spot,PCs)
result <- STIE(ST_expr, Signature, cells_on_spot2, features=c('angle','shape'), 
               lambda=0, steps=30, known_signature=FALSE, known_cell_types=FALSE)

cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
colors = c("green", "magenta", "orange", "blue", "black", "cyan")
colors = c("green", "magenta", "black", "orange", "blue", "cyan")
plot_sub_image(im=im, w=9000, h=5000, xoff=15500, yoff=11500,
               x_scale=x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=T, plot_cell=T)




