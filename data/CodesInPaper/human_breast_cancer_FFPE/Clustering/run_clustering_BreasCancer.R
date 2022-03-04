deconvolution = FALSE
clustering = FALSE
signature_learning = FALSE

library(STIE)
STIE.dir = system.file(package = "STIE")
load( paste0(STIE.dir,"/data/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R") )

############################################################
## run clustering
############################################################

cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)

result = list()
for(i in 2:10)
{
    cat(i,'\n')
    
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    Signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    
    workdir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/Kmeans/Seurat"
    fileName = paste0(workdir,"/cluster_",i,"/cluster_",i,".Seurat.markers.txt") 
    marker = read.delim( fileName, sep="\t", header=T)
    Signature = Signature[rownames(Signature)%in%marker$gene, ]
    
    result[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=1e6, steps=30, 
                     morphology_steps=ceiling(steps/3), 
                     known_signature=FALSE, known_cell_types=FALSE)
}

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
save(result, file="HumanBreastCancer_clustering.RData")

score1 = lapply( result[2:10], function(x) calculate_BIC(x, ST_expr) )
score2 = sapply(score1,function(x)mean(x$bic))
names(score2) = names(score2) = rr
plot(score2)
rr[which.min(score2)]


i = 10
markers = find_sig_markers(STIE_result=result[[i]], ST_expr, transform="log", DEG_pthres=1)
features = unique(markers$features)

markers2 = do.call(cbind, lapply( 1:i, function(j) {
    a = markers$pct.1[ markers$clusters == j ]
    names(a) = markers$features[ markers$clusters == j ]
    a = a[order(names(a))]
    a
} ) )

Signature = Signature[rownames(Signature) %in% rownames(markers2), ]
markers2 = markers2[ match( rownames(Signature), rownames(markers2) ), ]




############################################################
# run clustering on known signature genes
############################################################
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
SS = as.matrix(Signature)[rownames(Signature)%in%colnames(ST_expr),]

res4 = list()
for(i in 2:10)
{
    cat(i,'\n')
    
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    Signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    Signature = Signature[rownames(Signature)%in%rownames(SS), ]
    
    res4[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                     morphology_steps=ceiling(steps/3), 
                     known_signature=FALSE, known_cell_types=FALSE)
}

score1 = lapply( res4[2:10], function(x) BIC(x) )
score2 = sapply(score1,mean)
names(score2) = names(score2) = rr
plot(score2)


############################################################
# visualization
############################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")

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

