deconvolution = FALSE
clustering = TRUE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_cortex/parameters_cortex.R")


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

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE")
load("MouseBrainCortex_clustering.RData")

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/allen_cortex")
x = load("AllenCortex_MouseBrain_scRNASeq_DWLS_Signature.RData")
sc_signature = Signature[rownames(Signature)%in%colnames(ST_expr),]
sc_signature = sc_signature/100
sc_signature = sc_signature[,order(colnames(sc_signature))]


for(i in 2:10)
{
    cat(i,'\n')
    # i=9
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    kmeans_signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    stie_signature = result[[i]]$Signature
    
    ###### bayes space
    b = load( paste0("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/BayesSpace/cluster",
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
    
    kmeans_ols = apply( kmeans_signature, 2, function(b) solveOLS( sc_signature, b ) )
    bs_ols = apply( bs_signature, 2, function(b) solveOLS( sc_signature, b ) )
    stie_ols = apply( stie_signature, 2, function(b) solveOLS( sc_signature, b ) )
    
    stie_dwls = dwls(S=sc_signature, B=stie_signature)
    bs_dwls = dwls(S=sc_signature, B=bs_signature)
    kmeans_dwls = dwls(S=sc_signature, B=kmeans_signature)
    
    
    
    
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE")
    
    
    myCol = c('palegreen4','green3','green4','green','royalblue','royalblue4','red4','red3','red2')
    myCol=c('yellow','black','red','green','blue','purple','cyan','white','orange')
    
    myCol=c('yellow','black','red3','green3','green4','royalblue4','steelblue','royalblue')
    
    layout(matrix(1))
    par(mar=c(2,2,2,2))
    par(mfrow=c(2,3))
    index = grep("L[0-9]|Astro",colnames(sc_signature))
    colnames(sc_signature)[index]
    barplot(stie_ols[index,],col=myCol,ylim=c(0,1))
    barplot(bs_ols[index,],col=myCol,ylim=c(0,1))
    barplot(kmeans_ols[index,],col=myCol,ylim=c(0,1))
    
    barplot(stie_dwls[index,],col=myCol,ylim=c(0,1))
    barplot(bs_dwls[index,],col=myCol,ylim=c(0,1))
    barplot(kmeans_dwls[index,],col=myCol,ylim=c(0,1))
    
    
    
    
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
    
    
    x_scale = 0.1
    cell_types = result[[i]]$cell_types
    contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
    myCol=c('yellow','black','red3','green3','green4','royalblue4','steelblue','royalblue')
    
    myCol = c( "red", "blue", "purple", "green", "black", "yellow", "white", "white" )
    colors = myCol[ 1:ncol(result[[i]]$Signature) ]
    colors = colors[ length(colors)+1-rank(table(cell_types)) ]
    plot_sub_image(im=im, image_transparency=0, w=14000, h=7000, xoff=6000, yoff=8000, 
                   x_scale=x_scale, spot_coordinates=spot_coordinates, 
                   contour=contour2, cell_types=cell_types, color_use=colors, 
                   plot_spot=F, plot_cell=T  )
    
    
}

