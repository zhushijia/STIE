deconvolution = FALSE
clustering = TRUE
signature_learning = FALSE

STIE.dir = system.file(package = "STIE")
source( paste0(STIE.dir,"/data/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo.R") )

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
                     known_signature=FALSE, known_cell_types=FALSE)
}

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
save(result, file="MouseBrainHippocampus_clustering.RData")

score1 = lapply( result[2:9], function(x) calculate_BIC(x, ST_expr) )
score2 = sapply(score1,function(x) mean(x$bic) )
plot(score2)
rr[which.min(score2)]



> mean(b5[[1]])
[1] 63203.42
> mean(b6[[1]])
[1] 65479.65
> mean(b4[[1]])
[1] 63405.09

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
# clustering deconvolution
############################################################


dwls = function(S, B)
{
    library(Matrix)
    source("/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/DWLS/Deconvolution_functions.R")
    
    nb = ncol(B)
    ns = ncol(S)
    prop_mat <- as.data.frame(matrix(0, nrow = nb, ncol = ns) )
    rownames(prop_mat) <- colnames(B)
    colnames(prop_mat) <- colnames(S)
    for (ib in 1:nb) {
        print( sprintf("Estimating proportion for spot : %d / %d", ib, nb) )
        b <- B[,ib]
        tr <- trimData(S, b)
        tr$sig <- tr$sig[,colSums(tr$sig) > 0]
        is_pd <- eigen(t(tr$sig)%*%tr$sig)$values
        is_pd <- all(is_pd > 10e-6)
        if (!(is_pd)) { next}
        try(solDWLS <- solveDampenedWLS(tr$sig,tr$bulk),next)
        print("Proportions >> ")
        prop_mat[ib,names(solDWLS)] <- solDWLS
    }
    
    t(prop_mat)
    
}


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
Signature = read.delim("Major_cell_types_marker_genes.txt", header=T, row.names=1)
Signature = Signature[,!colnames(Signature)%in%c("Ependymal.cells")]
Signature = as.matrix(Signature)[rownames(Signature)%in%colnames(ST_expr),]
colnames(Signature) = gsub("Pyramidal.neurons.|.interneurons|.like.cells|Granule.cells.", "", colnames(Signature))
sc_signature = Signature[,order(colnames(Signature))]

prop = list()
for(i in 2:10)
{
    cat(i,'\n')
    # i==7, has CA2 great
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    cluster_signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    stie_signature = result[[i]]$Signature
    
    genes = intersect( rownames(stie_signature) , intersect( rownames(cluster_signature), rownames(sc_signature) ) )
    sc_signature = sc_signature[match(genes,rownames(sc_signature)), ]
    cluster_signature = cluster_signature[match(genes,rownames(cluster_signature)), ]
    stie_signature = stie_signature[match(genes,rownames(stie_signature)), ]
    
    cluster_prop = apply( cluster_signature, 2, function(b) solveOLS( sc_signature, b ) )
    stie_prop = apply( stie_signature, 2, function(b) solveOLS( sc_signature, b ) )
    
    par(mfrow=c(2,2))
    myCol=c('yellow','black','red','green','blue','purple')
    barplot(cluster_prop,col=myCol)
    barplot(stie_prop,col=myCol)
    
    snames = colnames(sc_signature)
    data.frame( apply(cluster_prop, , function(x) snames[which.max(x)] ), apply(cluster_prop, 2, max) )
    data.frame( apply(stie_prop, 2, function(x) snames[which.max(x)] ), apply(stie_prop, 2, max) )
    
    cluster_max = apply( cluster_prop, 1, max )
    stie_max = apply( stie_prop, 1, max )
    boxplot(cluster_max, stie_max)
    
    
    stie_dwls = dwls(S=sc_signature, B=stie_signature)
    cluster_dwls = dwls(S=sc_signature, B=cluster_signature)
    
    cluster_max = apply( cluster_prop, 2, max )
    stie_max = apply( stie_prop, 2, max )
    boxplot(cluster_max, stie_max)
    
    layout(matrix(1))
    par(mar=c(2,2,2,2))
    par(mfrow=c(2,2))
    index = 1:ncol(sc_signature)
    barplot(stie_prop[index,],col=myCol,ylim=c(0,1))
    barplot(cluster_prop[index,],col=myCol,ylim=c(0,1))
    barplot(stie_dwls[index,],col=myCol,ylim=c(0,1))
    barplot(cluster_dwls[index,],col=myCol,ylim=c(0,1))
    
    
    barplot(apply(stie_dwls[index,],2,function(x)x/sum(x)),col=myCol,ylim=c(0,1))
    barplot(apply(cluster_dwls[index,],2,function(x)x/sum(x)),col=myCol,ylim=c(0,1))
    
    cluster_max = apply( cluster_prop[index,], 1, max )
    stie_max = apply( stie_prop[index,], 1, max )
    boxplot(cluster_max, stie_max)
    wilcox.test(cluster_max, stie_max,'less',paired=T)
    
    cluster_max = apply( cluster_prop[index,1:6], 1, max )
    stie_max = apply( stie_prop[index,1:6], 1, max )
    boxplot(cluster_max, stie_max)
    wilcox.test(cluster_max, stie_max,'less',paired=T)
    
    
    
}



############################################################
# visualization
############################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
load("MouseBrainHippocampus_clustering.RData")

x_scale = 0.1
pdf( paste0("MouseBrainHippocampus_clustering_scale", x_scale, ".pdf") )

for(i in 2:length(result))
{
    cell_types = result[[i]]$cell_types
    contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
    
    myCol = c( "red", "blue", "green", "black", "cyan", "yellow", "purple", "steelblue",
               "darkorange", "darkred", "grey", "blue4", "chartreuse4", "burlywood1", "darkgoldenrod4" )
    colors = myCol[ 1:ncol(result[[i]]$Signature) ]
    
    colors = get_my_colors(ncol(result[[i]]$Signature),mode=2)
    
    plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
                   x_scale=x_scale, spot_coordinates=spot_coordinates, 
                   contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )
}

dev.off()
