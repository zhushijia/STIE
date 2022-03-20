deconvolution = FALSE
clustering = TRUE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo.R")

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

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
load("MouseBrainHippocampus_clustering.RData")

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
Signature = read.delim("Major_cell_types_marker_genes.txt", header=T, row.names=1)
Signature = Signature[,!colnames(Signature)%in%c("Ependymal.cells")]
Signature = as.matrix(Signature)[rownames(Signature)%in%colnames(ST_expr),]
colnames(Signature) = gsub("Pyramidal.neurons.|.interneurons|.like.cells|Granule.cells.", "", colnames(Signature))
sc_signature = Signature[,order(colnames(Signature),decreasing=T)]

prop = list()
for(i in 2:10)
{
    cat(i,'\n')
    # i==7, has CA2 great
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    kmeans_signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    stie_signature = result[[i]]$Signature
    
    ###### bayes space
    b = load( paste0("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/BayesSpace/cluster",
                     "/cluster", i, "/cluster",i,".RData") )
    bs_signature = t(apply( counts, 1, function(x) {
        tapply(x, cluster, function(y) mean(exp(y),na.rm=T) )   
    } ))
    
    ################################################################################################
    genes = intersect( rownames(kmeans_signature), rownames(sc_signature) )
    genes = intersect( rownames(bs_signature) , genes )
    genes = intersect( rownames(stie_signature) , genes )
    
    sc_signature = sc_signature[match(genes,rownames(sc_signature)), ]
    kmeans_signature = kmeans_signature[match(genes,rownames(kmeans_signature)), ]
    bs_signature = bs_signature[match(genes,rownames(bs_signature)), ]
    stie_signature = stie_signature[match(genes,rownames(stie_signature)), ]
    
    stie_ols = apply( stie_signature, 2, function(b) solveOLS( sc_signature, b ) )
    bs_ols = apply( bs_signature, 2, function(b) solveOLS( sc_signature, b ) )
    kmeans_ols = apply( kmeans_signature, 2, function(b) solveOLS( sc_signature, b ) )
    
    stie_dwls = dwls(S=sc_signature, B=stie_signature)
    bs_dwls = dwls(S=sc_signature, B=bs_signature)
    kmeans_dwls = dwls(S=sc_signature, B=kmeans_signature)
    
    
    ################################################################################################
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
    pdf( paste0("STIE_cluster_",i,".pdf"), w=10, h=8 )
    
    myCol=c('yellow','black','red','green','blue','purple','cyan','white','orange')
    myCol=rev( c('red','blue','green','black','cyan','yellow') )
    
    layout(matrix(1))
    par(mar=c(2,2,2,2))
    par(mfrow=c(2,3))
    index = 1:ncol(sc_signature)
    index = c(7,5,2,3,6,1,4)
    index = c(3,6,2,4,8,1,7,5)
    barplot(stie_ols[,index],col=myCol,ylim=c(0,1))
    barplot(bs_ols,col=myCol,ylim=c(0,1))
    barplot(kmeans_ols,col=myCol,ylim=c(0,1))
    barplot(stie_dwls[,index],col=myCol,ylim=c(0,1))
    barplot(bs_dwls,col=myCol,ylim=c(0,1))
    barplot(kmeans_dwls,col=myCol,ylim=c(0,1))
    
    ############################################################
    # visualization
    ############################################################
    
    cell_types = result[[i]]$cell_types
    contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
    colors = get_my_colors(ncol(result[[i]]$Signature),mode=2)
    colors = c('yellow','green','black','cyan','blue','red','white')
    colors = c('orange','green','red','black','yellow','blue','cyan','purple')
    plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
                   x_scale=x_scale, spot_coordinates=spot_coordinates, 
                   contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )
    
    dev.off()
    
    
    
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
