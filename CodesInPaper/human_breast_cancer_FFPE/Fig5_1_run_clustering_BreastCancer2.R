deconvolution = FALSE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

############################################################
## run clustering
############################################################
ri = 2.5
la = 1e3
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, ri*spot_radius)

results = list()
for(i in 2:10)
{
    cat(i,'\n')
    
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    Signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    
    results[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=la, steps=30, 
                     known_signature=FALSE, known_cell_types=FALSE)
}

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
names(results)[2:10] = 2:10
#save(results, file=paste0("HumanBreastCancer_clustering_",ri,"xSpot.RData") )
save(results, file=paste0("HumanBreastCancer_clustering_",ri,"xSpot_lambda",la,".RData") )



############################################################
# visualization
############################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("HumanBreastCancer_clustering.RData")

pdf("HumanBreastCancer_clustering.pdf")

for(i in 2:10)
{
    cell_types = result[[i]]$cell_types
    contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
    
    myCol = c( "red", "blue", "green", "black", "cyan", "yellow", "purple", "steelblue",
               "darkorange", "darkred", "grey", "blue4", "chartreuse4", "burlywood1", "darkgoldenrod4" )
    colors = myCol[1:ncol(result[[i]]$Signature)]
    
    #colors = get_my_colors(ncol(Signature),mode=2)
    
    plot_sub_image(im=im, 
                   x_scale=args$x_scale, spot_coordinates=spot_coordinates, 
                   contour=contour2, cell_types=cell_types, color_use=colors, 
                   plot_spot=F, plot_cell=T, 
                   axis_tick=0, axis_col='grey'  )
}

dev.off()





plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
               x_scale=1, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )

