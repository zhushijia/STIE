source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

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

table2df <- function(t) {
    data.frame(ID=rownames(t), do.call(rbind, lapply(1:nrow(t),function(i) t[i,])) )
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
#sc_signature = sc_signature[, ! colnames(sc_signature) %in% c("Myeloid","PVL") ]
sc_group = data.frame( celltype = c( "Bcells", "Tcells", "Plasmablasts", "Myeloid", "CancerEpithelial", "NormalEpithelial", "CAFs", "Endothelial", "PVL"),
                      category = c( "Immune", "Immune", "Immune", "Immune", "Epithelial", "Epithelial", "Stroma", "Stroma", "Stroma"),
                       myCol = c('palegreen4','green3','green4','green','royalblue','royalblue4','red4','red3','red2'))
sc_group = sc_group[sc_group$celltype %in% colnames(sc_signature), ]
sc_signature = sc_signature[ , match( sc_group$celltype, colnames(sc_signature) ) ]

############################################################
# clustering deconvolution
############################################################

for(i in 2:10)
{
    cat(i,'\n')
    # i=6
    ##### kmeans
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    kmeans_signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    stie_signature = stie_clustering[[i]]$Signature
    
    ###### bayes space
    b = load( paste0("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/BayesSpace/cluster",
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
    
    kmeans_ols = apply( kmeans_signature, 2, function(b) solveNNLS( sc_signature, b ) )
    bs_ols = apply( bs_signature, 2, function(b) solveNNLS( sc_signature, b ) )
    stie_ols = apply( stie_signature, 2, function(b) solveNNLS( sc_signature, b ) )
    
    kmeans_dwls = dwls(S=sc_signature, B=kmeans_signature)
    bs_dwls = dwls(S=sc_signature, B=bs_signature)
    stie_dwls = dwls(S=sc_signature, B=stie_signature)
    
    myCol=c('yellow','black','red','green','blue','purple','cyan','white','orange')
    myCol=c('palegreen4','green3','green4','royalblue','royalblue4','red4','red3')
    myCol=as.character(sc_group$myCol)
    
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
    pdf( paste0("overlap_STIE_deconvolute_and_STIE_cluster_",i,".pdf"), w=10, h=8 )
    
    index1 = order( apply(stie_ols,2,max), decreasing=T)
    index2 = order( apply(bs_ols,2,max), decreasing=T)
    index1 = c(1,6,3,2,4,5)
    index2 = c(5,6,3,4,1,2)
    index3 = c(1,6,3,2,4,5)
    
    layout(matrix(1))
    par(mar=c(2,2,2,2))
    par(mfrow=c(2,3))
    barplot(stie_ols[,index1], col=myCol, ylim=c(0,1), border='white')
    barplot(bs_ols[,index2], col=myCol, ylim=c(0,1), border='white')
    barplot(kmeans_ols[,index3], col=myCol, ylim=c(0,1), border='white')
    
    barplot(stie_dwls[,index1], col=myCol, ylim=c(0,1), border='white')
    barplot(bs_dwls[,index2], col=myCol, ylim=c(0,1), border='white')
    barplot(kmeans_dwls[,index3], col=myCol, ylim=c(0,1), border='white')

    ############################################################
    
    d = 1
    kmeans_max = apply( kmeans_ols, d, max)
    bs_max = apply( bs_ols, d, max)
    stie_max = apply( stie_ols, d, max)
    
    wilcox.test( stie_max, bs_max, 'greater', paired=T )
    wilcox.test( stie_max, kmeans_max, 'greater', paired=T )
    
    boxplot( kmeans_max, bs_max, stie_max )
    n = length(kmeans_max)
    points( 1+0.5*runif(n,-0.5,0.5), kmeans_max, col=myCol, pch=15, cex=4)
    points( 2+0.5*runif(n,-0.5,0.5), bs_max, col=myCol, pch=15, cex=4)
    points( 3+0.5*runif(n,-0.5,0.5), stie_max, col=myCol, pch=15, cex=4)
    
    ############################################################
    if(1) {
        de_i = stie_deconvolution[['2']]
        cl_i = stie_clustering[[i]]
        
        cells = unique( intersect( names(de_i$cell_types), names(cl_i$cell_types) ) )
        de_i_uni_celltypes = de_i$cell_types[ match( cells, names(de_i$cell_types) ) ]
        cl_i_uni_celltypes = cl_i$cell_types[ match( cells, names(cl_i$cell_types) ) ]
        
        #t1 = table(de_i_uni_celltypes, cl_i_uni_celltypes)
        #t1 = t1[match( colnames(sc_signature), rownames(t1)), index1]
        
        de_i_uni_celltypes2 = as.character(sc_group$category)[ match( de_i_uni_celltypes, as.character(sc_group$celltype) ) ]
        t2 = table(de_i_uni_celltypes2, cl_i_uni_celltypes)
        t2 = t2[,index1]
        #barplot(t1,col=myCol,border='white')
        barplot(t2,col=c("steelblue",'darkgreen','darkred'),border='white')
        plot(NA,xlim=c(1,10),ylim=c(1,10))
        legend('topright',col=myCol, pch=15, cex=2, legend=colnames(sc_signature))
        
        #write.table(table2df(t1), paste0("overlap_STIE_deconvolute_and_STIE_cluster_",i,".txt"),
        #            sep="\t",col.names=T,row.names=F,quote=F)
        write.table(table2df(t1), paste0("overlap_large_category_STIE_deconvolute_and_STIE_cluster_",i,".txt"),
                    sep="\t",col.names=T,row.names=F,quote=F)
    }
    
    ############################################################
    cell_types = cl_i$cell_types
    contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
    
    colors = c("steelblue",'darkgreen',"sienna1",'mediumpurple','darkred','brown')
    colors = c("steelblue",'darkgreen',"darkorange",'darkred','purple','brown')
    
    plot_sub_image(im=im, image_transparency=0,
                   x_scale=args$x_scale, spot_coordinates=spot_coordinates, 
                   contour=contour2, cell_types=cell_types, color_use=colors, 
                   plot_spot=F, plot_cell=T, 
                   axis_tick=0, axis_col='grey'  )
    
    dev.off()
    
}
