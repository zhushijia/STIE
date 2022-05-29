deconvolution = FALSE
clustering = TRUE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_cortex/parameters_cortex.R")

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
load("MouseBrainCortex_clustering_2.5xSpot_lambda0.RData")

cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/allen_cortex")
x = load("AllenCortex_MouseBrain_scRNASeq_DWLS_Signature.RData")
sc_signature = Signature[rownames(Signature)%in%colnames(ST_expr),]
sc_signature = sc_signature/100
sc_signature = sc_signature[,order(colnames(sc_signature))]


for(i in 2:10)
{
    cat(i,'\n')
    # i=6 is the best
    # i=6 
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    kmeans_signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    stie_signature = results[[i]]$Signature
    
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
    layout(matrix(1))
    par(mar=c(2,2,2,2))
    par(mfrow=c(2,3))
    index = grep("L[0-9]|Astro|Oligo",colnames(sc_signature))
    colnames(sc_signature)[index]
    myCol=get_my_colors(i,mode=3)
    barplot(kmeans_ols[index,],col=myCol,ylim=c(0,1))
    barplot(bs_ols[index,],col=myCol,ylim=c(0,1))
    barplot(stie_ols[index,],col=myCol,ylim=c(0,1))
    
    barplot(kmeans_dwls[index,],col=myCol,ylim=c(0,1))
    barplot(bs_dwls[index,],col=myCol,ylim=c(0,1))
    barplot(stie_dwls[index,],col=myCol,ylim=c(0,1))
    
    outputDir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE/clustering_visualization"
    dir.create(outputDir, recursive = T)
    setwd(outputDir)
    
    
    if(i==6)
    {
        pdf( paste0("STIE_cluster_",i,"_signature_deconvolution.pdf"), w=10, h=8 )
        layout(matrix(1))
        par(mar=c(2,2,2,2))
        par(mfrow=c(2,3))
        index = grep("L[0-9]|Astro|Oligo",colnames(sc_signature))
        myCol=c('grey','orange','darkred','green2','green3','dodgerblue','dodgerblue2','dodgerblue3','purple')
        
        ref_col = c( "cyan", "black", "red", "blue", "steelblue", "yellow")
        kmeans_col = c("red", "steelblue", "blue", "black", "yellow", "cyan")
        bs_col = c("yellow", "blue", "red", "steelblue", "black", "cyan")
        stie_col = c( "red", "steelblue", "blue", "black", "yellow", "cyan")
        
        kmeans_index = c()
        kmeans_index = match(ref_col, kmeans_col)
        bs_index = match(ref_col, bs_col)
        stie_index = match(ref_col, stie_col)
        
        barplot(kmeans_ols[index,][, kmeans_index],col=myCol,ylim=c(0,1), border='white')
        barplot(bs_ols[index,][, bs_index],col=myCol,ylim=c(0,1), border='white')
        barplot(stie_ols[index,][, stie_index],col=myCol,ylim=c(0,1), border='white')
        
        barplot(kmeans_dwls[index,][, kmeans_index],col=myCol,ylim=c(0,1), border='white')
        barplot(bs_dwls[index,][, bs_index],col=myCol,ylim=c(0,1), border='white')
        barplot(stie_dwls[index,][, stie_index],col=myCol,ylim=c(0,1), border='white')
        dev.off()
        
    }
    
    
}


