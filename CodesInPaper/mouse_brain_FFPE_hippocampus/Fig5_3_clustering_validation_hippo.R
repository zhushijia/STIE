deconvolution = FALSE
clustering = TRUE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo3.R")

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

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")
load("MouseBrainHippocampus_clustering_2.5Xspot.RData")

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
Signature = read.delim("Major_cell_types_marker_genes.txt", header=T, row.names=1)
Signature = Signature[,!colnames(Signature)%in%c("Ependymal.cells")]
Signature = as.matrix(Signature)[rownames(Signature)%in%colnames(ST_expr),]
colnames(Signature) = gsub("Pyramidal.neurons.|.interneurons|.like.cells|Granule.cells.", "", colnames(Signature))
sc_signature = Signature[,order(colnames(Signature),decreasing=T)]

for(i in 2:10)
{
    cat(i,'\n')
    # i==5, has CA2 great
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    kmeans_signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    stie_signature = results[[i]]$Signature
    
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
    
    stie_ols = apply( stie_signature, 2, function(b) solveNNLS( sc_signature, b ) )
    bs_ols = apply( bs_signature, 2, function(b) solveNNLS( sc_signature, b ) )
    kmeans_ols = apply( kmeans_signature, 2, function(b) solveNNLS( sc_signature, b ) )
    
    library(nnls)
    s_nnls = function(A,b) { x=coef(nnls(A,b)); x/sum(x) }
    stie_nnls = apply( stie_signature, 2, function(b) s_nnls( sc_signature, b ) )
    bs_nnls = apply( bs_signature, 2, function(b) s_nnls( sc_signature, b ) )
    kmeans_nnls = apply( kmeans_signature, 2, function(b) s_nnls( sc_signature, b ) )
    
    stie_dwls = dwls(S=sc_signature, B=stie_signature)
    bs_dwls = dwls(S=sc_signature, B=bs_signature)
    kmeans_dwls = dwls(S=sc_signature, B=kmeans_signature)
    
    ################################################################################################
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE/clustering_visualization")
    pdf( paste0("STIE_cluster_",i,"_signature_deconvolution.pdf"), w=10, h=8 )
    
    myCol=rev( c('magenta','blue','green','black','cyan','orange') )
    
    layout(matrix(1))
    par(mar=c(2,2,2,2))
    par(mfrow=c(2,3))
    if (i==5) {
        stie_index = c(4,5,2,3,1)
        bs_index = c(4,3,2,5,1)
        kmeans_index = c(4,5,2,3,1)
    } else {
        stie_index = 1:ncol(stie_ols)
        bs_index = 1:ncol(bs_ols)
        kmeans_index = 1:ncol(kmeans_ols)
    }   
    
    barplot(kmeans_ols[,kmeans_index], col=myCol, ylim=c(0,1))
    barplot(bs_ols[,bs_index], col=myCol, ylim=c(0,1))
    barplot(stie_ols[,stie_index], col=myCol, ylim=c(0,1))
    
    if(0) {
        barplot(kmeans_nnls[,kmeans_index], col=myCol, ylim=c(0,1))
        barplot(bs_nnls[,bs_index], col=myCol, ylim=c(0,1))
        barplot(stie_nnls[,stie_index], col=myCol, ylim=c(0,1))
    }
    
    barplot(kmeans_dwls[,kmeans_index], col=myCol, ylim=c(0,1))
    barplot(bs_dwls[,bs_index], col=myCol, ylim=c(0,1))
    barplot(stie_dwls[,stie_index], col=myCol, ylim=c(0,1))
    
    ############################################################
    # visualization
    ############################################################
    
    cell_types = results[[i]]$cell_types
    contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
    
    if(i==5) {
        colors = c('magenta','blue','green','black','orange')[order(stie_index)]
    } else {
        colors = get_my_colors(i, mode=2)
    }
    
    #png("Hippo_cluster5.png", w=2000, h=2000, res=300)
    plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
                   x_scale=0.1, spot_coordinates=spot_coordinates, 
                   contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=F, plot_cell=T  )
    
    dev.off()
    
}


