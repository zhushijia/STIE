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

############################################################
# clustering deconvolution
############################################################

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("HumanBreastCancer_clustering.RData")
stie_clustering = result

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("BreastCancer_spot_BIC_new_get_cells_on_spot.RData")
stie_deconvolution = result

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
sc_signature = Signature[,order(colnames(Signature))]

############################################################
# clustering deconvolution
############################################################

de_i = stie_deconvolution[['2']]
cl_i = stie_clustering[[6]]

cells = unique( intersect( names(de_i$cell_types), names(cl_i$cell_types) ) )

de_i_uni_celltypes = de_i$cell_types[ match( cells, names(de_i$cell_types) ) ]
cl_i_uni_celltypes = cl_i$cell_types[ match( cells, names(cl_i$cell_types) ) ]

table(de_i_uni_celltypes, cl_i_uni_celltypes)


############################################################
# clustering deconvolution
############################################################



for(i in 2:10)
{
    cat(i,'\n')
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
    myCol=c('yellow','black','red','green','blue','purple','cyan','white','orange')
    barplot(stie_prop,col=myCol)
    barplot(cluster_prop,col=myCol)
    
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
    barplot(stie_prop[index,],col=myCol)
    barplot(cluster_prop[index,],col=myCol)
    barplot(stie_dwls[index,],col=myCol,ylim=c(0,1))
    barplot(cluster_dwls[index,],col=myCol,ylim=c(0,1))
    
    
    
    barplot(apply(stie_dwls[index,],2,function(x)x/sum(x)),col=myCol,ylim=c(0,1))
    barplot(apply(cluster_dwls[index,],2,function(x)x/sum(x)),col=myCol,ylim=c(0,1))
    
    cluster_max = apply( cluster_prop[index,], 1, max )
    stie_max = apply( stie_prop[index,], 1, max )
    boxplot(cluster_max, stie_max)
    wilcox.test(cluster_max, stie_max,'less',paired=T)
    
    cluster_max = apply( cluster_dwls[index,], 1, max )
    stie_max = apply( stie_dwls[index,], 1, max )
    boxplot(cluster_max, stie_max)
    wilcox.test(cluster_max, stie_max,'less',paired=T)
    

    
}


