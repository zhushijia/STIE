
buildSeurat <- function( countData, celltypes, projectName, outputDir )
{
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
    p1 <- DimPlot(object = Seurat, reduction = "umap", label = TRUE)
    
    celltypes2 = apply(probs, 2, function(x) rownames(probs)[which.max(x)] )
    celltypes2 = res$cell_labels
    Idents(Seurat) = as.character(celltypes2)
    p2 <- DimPlot(object = Seurat, reduction = "umap", label = TRUE)
    
    celltypes3 = apply(probs, 2, function(x) rownames(probs)[which.max(x)] )
    celltypes3 = res$cell_labels
    Idents(Seurat) = as.character(celltypes3)
    p3 <- DimPlot(object = Seurat, reduction = "umap", label = TRUE)
    
    plot_grid(p1,p2,p3)
    
}


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
        x = x[ 1:min(100,length(x)) ]
        x
    } ) 
    
    
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE/tmp")
    res <- SCINA(exp, signatures_list, max_iter = 100, convergence_n = 10, 
                     convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
                     rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
    probs = res$probabilities

}

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

tmp = lapply(uni_celltypes, function(x) {
    index = which(celltypes==x)
    sample( index, ceiling(length(index)*0.15) )
} )
selected = c()
for(k in 1:length(tmp)) selected = c(selected, tmp[[k]])

countData = counts[,selected]
celltypes1 = celltypes[selected]


###########################################################################
###########################################################################
###########################################################################
workdir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq"
projectName = "Wu_etal_2021_BRCA_scRNASeq"
fileName = paste0(workdir,"/Seurat/",projectName,".Seurat.markers.txt") 
truth = read.delim( fileName, sep="\t", header=T)
truth = tapply(truth$gene, truth$cluster, as.character)


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("HumanBreastCancer_clustering.RData")

i = 6
markers = find_sig_markers(STIE_result=result[[i]], ST_expr, transform="log", DEG_pthres=1)

if(0)
{
    x = result[[i]]$Signature
    y = counts
    g = intersect(rownames(x),rownames(y))
    x = x[match(g,rownames(x)),]
    y = y[match(g,rownames(y)),]
    probs = cor(as.matrix(x),as.matrix(log2(y+1)),method='pearson')
}

markers2 = subset(markers, pvalues<1e-10 & diff>0)
#markers2 = tapply(markers2$features, markers2$clusters, as.character)
markers2 = lapply( split(markers, markers$clusters), function(x) {
    x = x[order(x$pvalues),]
    cat(x$pvalues[50],"\n")
    ss = sum(x$pvalues<1)
    x$features[1:150]
})
sapply(markers2,length)


overlap = matrix( nrow=length(markers2), ncol=length(truth), data=0 )
rownames(overlap) = names(markers2)
colnames(overlap) = names(truth)

for(i in 1:length(markers2) )
{
    for(j in 1:length(truth) )
    {
        overlap[i,j] = sum( markers2[[i]]%in%truth[[j]] )
    }
}

cbind(overlap, sapply(markers2,length) )

x = do.call(rbind, lapply(1:nrow(overlap), function(i) signif( overlap[i,]/length(markers2[[i]]),2) ))
x

workdir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/Kmeans/Seurat"
fileName = paste0(workdir,"/cluster_",i,"/cluster_",i,".Seurat.markers.txt") 
cluster = read.delim( fileName, sep="\t", header=T)
cluster = tapply(cluster$gene, cluster$cluster, as.character)


overlap2 = matrix( nrow=length(cluster), ncol=length(truth), data=0 )
rownames(overlap2) = names(cluster)
colnames(overlap2) = names(truth)

for(i in 1:length(cluster) )
{
    for(j in 1:length(truth) )
    {
        overlap2[i,j] = sum( cluster[[i]]%in%truth[[j]] )
    }
}

cbind(overlap2, sapply(cluster,length) )
y = do.call(rbind, lapply(1:nrow(overlap2), function(i) signif( overlap2[i,]/length(cluster[[i]]),2) ))


################################################################################################











