deconvolution = FALSE
clustering = TRUE
signature_learning = FALSE

STIE.dir = system.file(package = "STIE")
load( paste0(STIE.dir,"/data/CodesInPaper/XXX.R") )

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
    
    result[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                     morphology_steps=ceiling(steps/3), 
                     known_signature=FALSE, known_cell_types=FALSE)
}

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/count/results/STIE")
save(result, file="MouseKidney_clustering.RData")

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

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
Signature = read.delim("Major_cell_types_marker_genes.txt", header=T, row.names=1)
SS = Signature[,!colnames(Signature)%in%c("Ependymal.cells")]
SS = as.matrix(SS)[rownames(SS)%in%colnames(ST_expr),]

res4 = list()
for(i in 2:10)
{
    cat(i,'\n')
    
    cluster = read.csv(paste0("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/Visium_FFPE_Mouse_Brain/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
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

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseKidney_FFPE/count/results/STIE")
load("MouseKidney_clustering.RData")

x_scale = args$x_scale
pdf( paste0("MouseKidney_clustering_scale", x_scale, ".pdf") )

for(i in 2:length(result))
{
    cell_types = result[[i]]$cell_types
    contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
    
    myCol = c( "red", "blue", "green", "black", "cyan", "yellow", "steelblue",
               "darkorange", "darkred", "grey", "blue4", "purple", "chartreuse4", "burlywood1", "darkgoldenrod4" )
    colors = myCol[ 1:ncol(result[[i]]$Signature) ]
    
    colors = get_my_colors(ncol(result[[i]]$Signature),mode=2)
    
    plot_sub_image(im=im, 
                   x_scale=x_scale, spot_coordinates=spot_coordinates, 
                   contour=contour2, cell_types=cell_types, color_use=colors, 
                   plot_spot=F, plot_cell=T, 
                   axis_tick=0, axis_col='grey'  )
}

dev.off()





