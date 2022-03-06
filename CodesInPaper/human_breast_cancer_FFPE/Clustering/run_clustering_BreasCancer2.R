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
#for(i in 2:10)
#{
    i = 9
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
#}



###########################################################################
###########################################################################
###########################################################################
workdir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq"
setwd(workdir)
annot = read.csv("metadata.csv",header=T,row.names=1)
counts <- readMM("count_matrix_sparse.mtx")
genes <- read.delim("count_matrix_genes.tsv",header=F)
cell_ids <- read.delim("count_matrix_barcodes.tsv",header=F)
rownames(counts) <- as.character(genes[,1])
colnames(counts) <- as.character(cell_ids[,1])

counts = counts[,match( rownames(annot), colnames(counts) )]
all(rownames(annot)==colnames(counts))
celltypes = gsub(" |-", "", as.character(annot$celltype_major) )

########## downsampling
uni_celltypes = unique(celltypes)
set.seed(1234567)
#selected = do.call(c, lapply(uni_celltypes, function(x) {
#    index = which(celltypes==x)
#    sample( index, ceiling(length(index)*0.15) )
#} ))

#uni_celltypes = setdiff(uni_celltypes, "NormalEpithelial")

tmp = lapply(uni_celltypes, function(x) {
    index = which(celltypes==x)
    sample( index, ceiling(length(index)*0.15) )
} )
selected = c()
for(k in 1:length(tmp)) selected = c(selected, tmp[[k]])

countData = counts[,selected]
celltypes1 = celltypes[selected]

i = 9
markers = find_sig_markers(STIE_result=result[[i]], ST_expr, transform="log", DEG_pthres=1)
markers2 = lapply( split(markers, markers$clusters), function(x) {
    x = x[order(x$pvalues),]
    k = sum(x$pvalues<1e-10)
    cat(x$pvalues[k],"\n")
    x$features[1:k]
    
})
sapply(markers2,length)

#markers = find_sig_markers(STIE_result=result[[i]], ST_expr, transform="log", DEG_pthres=1)
markers2 = lapply( split(markers, markers$clusters), function(x) {
    x = x[order(x$pvalues),]
    x = subset(x, diff>0 & pvalues<1e-2 )
    intersect( as.character(x$features) , as.character(cluster0$gene) )
})
sapply(markers2,length)

min.pct = 0.25, logfc.threshold = 0.25
###########################################################################
###########################################################################
###########################################################################

workdir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/Kmeans/Seurat"
fileName = paste0(workdir,"/cluster_",i,"/cluster_",i,".Seurat.markers.txt") 
cluster0 = read.delim( fileName, sep="\t", header=T)
cluster = tapply(cluster0$gene, cluster0$cluster, as.character)


###########################################################################
###########################################################################
###########################################################################

runSCINA <- function(countData, signatures_list)
{
    
    library(SCINA)
    library(preprocessCore)
    
    exp_raw = log(as.matrix(countData)+1)
    exp = normalize.quantiles(exp_raw)
    rownames(exp) = rownames(exp_raw)
    colnames(exp) = colnames(exp_raw)
    
    #signatures_list = cluster
    #signatures_list = markers2
    signatures_list = lapply(signatures_list, function(x) {
        x = intersect(x,rownames(countData))
        x = x[ 1:min(200,length(x)) ]
        x
    } ) 
    signatures_list = signatures_list[!sapply(signatures_list,function(x)any(is.na(x)))]
    
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE/tmp")
    res <- SCINA(exp, signatures_list, max_iter = 100, convergence_n = 10, 
                 convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
                 rm_overlap=TRUE, allow_unknown=FALSE, log_file='SCINA.log')
    probs = res$probabilities
    res
}

#buildSeurat <- function( countData, celltypes, projectName, outputDir )
#{
    Seurat <- CreateSeuratObject(counts = countData, project = 'test', min.cells = 0)
    Seurat <- NormalizeData(object = Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    Seurat <- FindVariableFeatures(object = Seurat, selection.method = "vst", nfeatures = 2000)
    dim(Seurat$RNA@counts)
    dim(Seurat$RNA@data)
    
    # Run the standard workflow for visualization and clustering
    Seurat <- ScaleData(object = Seurat, verbose = FALSE)
    Seurat <- RunPCA(object = Seurat, npcs = 30, verbose = FALSE)
    Seurat <- RunUMAP(object = Seurat, reduction = "pca", dims = 1:30)
    Seurat <- FindNeighbors(object = Seurat, reduction = "pca", dims = 1:30)
    
    Idents(Seurat) = as.character(celltypes1)
    #Seurat.markers <- FindAllMarkers(object = Seurat.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    
    p1 <- DimPlot(object = Seurat, reduction = "umap", label = TRUE)
    
    res2 = runSCINA(countData, cluster)
    celltypes2 = res2$cell_labels
    Idents(Seurat) = as.character(celltypes2)
    p2 <- DimPlot(object = Seurat, reduction = "umap", label = TRUE)
    
    res3 = runSCINA(countData, markers2)
    celltypes3 = res3$cell_labels
    Idents(Seurat) = as.character(celltypes3)
    p3 <- DimPlot(object = Seurat, reduction = "umap", label = TRUE)
    
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE/tmp")
    pdf(paste0('test_cluster',i,'.pdf'),w=20,h=10)
    plot_grid(p1,p2,p3)
    dev.off()
#}

    
    
    
    

    
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")    
x = load("Wu_etal_2021_BRCA_scRNASeq_meanExpr.RData")

cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
all( as.character(cluster[,1]) %in% rownames(ST_expr) )
ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
cluster_signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))

stie_signature = result[[i]]$Signature

genes = intersect( intersect( rownames(stie_signature) , rownames(cluster_signature)  ), rownames(meanExpr) )

meanExpr = meanExpr[ match(genes, rownames(meanExpr)), ]
cluster_signature = cluster_signature[ match(genes, rownames(cluster_signature) ), ]
stie_signature = stie_signature[ match(genes, rownames(stie_signature)), ]


cor(meanExpr, cluster_signature)
cor(meanExpr, stie_signature)


a = apply( cluster_signature, 2, function(x) solveOLS(meanExpr, x, scaled = T) )
b = apply( stie_signature, 2, function(x) solveOLS(meanExpr, x, scaled = T) )

a1 = apply( a, 2, function(x) colnames(meanExpr)[which.max(x)] )
a2 = apply( a, 2, function(x) max(x) )
data.frame(a1,a2)

b1 = apply( b, 2, function(x) colnames(meanExpr)[which.max(x)] )
b2 = apply( b, 2, function(x) max(x) )
data.frame(b1,b2)

