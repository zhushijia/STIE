source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

runSCINA <- function(countData, signatures_list, allow_unknown=TRUE)
{
    library(SCINA)
    library(preprocessCore)
    
    exp_raw = log(as.matrix(countData)+1)
    exp = normalize.quantiles(exp_raw)
    rownames(exp) = rownames(exp_raw)
    colnames(exp) = colnames(exp_raw)
    
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE/tmp")
    res <- SCINA(exp, signatures_list, max_iter = 100, convergence_n = 10, 
                 convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
                 rm_overlap=TRUE, allow_unknown=allow_unknown, log_file='SCINA.log')
    #probs = res$probabilities
    
    res
}


###########################################################################
###########################################################################
###########################################################################
library(Matrix)
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

downsampling = FALSE
if(downsampling)
{
    countData = counts[,selected]
    celltypes1 = celltypes[selected]
    myCol2 = c('black','cyan',"darkorange",'yellow',"steelblue",'darkgreen','darkred')
    myCol3 = c('cyan',"darkorange",'yellow',"steelblue",'darkgreen','darkred')
} else {
    countData = counts
    celltypes1 = celltypes
    myCol2 = c('black','cyan',"darkorange",'yellow','darkred',"steelblue",'darkgreen')
    myCol3 = myCol2[-1]
}


###########################################################################
###########################################################################
###########################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
x = load("HumanBreastCancer_clustering_2.5xSpot_lambda1000.RData")
markers = find_sig_markers(STIE_result=results[[6]], ST_expr, transform="log", DEG_pthres=1)
markers = subset(markers, features %in%rownames(countData))
signatures_list = lapply( split(markers, markers$clusters), function(x) {
    x = x[order(x$pvalues),]
    as.character(x$features)[1:100]
})

SCINA_res = runSCINA(countData, signatures_list)

###########################################################################
###########################################################################
###########################################################################
Seurat <- CreateSeuratObject(counts = countData, project = 'STIE_cluster6', min.cells = 0)
Seurat <- NormalizeData(object = Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
Seurat <- FindVariableFeatures(object = Seurat, selection.method = "vst", nfeatures = 2000)
dim(Seurat$RNA@counts)
dim(Seurat$RNA@data)

# Run the standard workflow for visualization and clustering
Seurat <- ScaleData(object = Seurat, verbose = FALSE)
Seurat <- RunPCA(object = Seurat, npcs = 30, verbose = FALSE)
Seurat <- RunUMAP(object = Seurat, reduction = "pca", dims = 1:30)
Seurat <- FindNeighbors(object = Seurat, reduction = "pca", dims = 1:30)



setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")

pdf("STIE_cluster_match_scRNAseq_allcells.pdf", w=15, h=5)

Idents(Seurat) = as.character(celltypes1)
p1.1 <- DimPlot(object = Seurat, reduction = "umap", label = TRUE)
p1.2 <- DimPlot(object = Seurat, reduction = "umap", label = FALSE)
plot_grid(p1.1, p1.2)


celltypes2 = SCINA_res$cell_labels
Idents(Seurat) = as.character(celltypes2)
p2.1 <- DimPlot(object = Seurat, reduction = "umap", label = TRUE, cols=myCol2)
p2.2 <- DimPlot(object = Seurat, reduction = "umap", label = FALSE, cols=myCol2)
plot_grid(p2.1, p2.2)

Seurat$knownTag = sapply( celltypes2, function(x) ifelse(x!="unknown","known","unknown") )
Seurat2 = subset(Seurat, subset = knownTag == "known" )
p3.1 <- DimPlot(object = Seurat2, reduction = "umap", label = TRUE, cols=myCol3)
p3.2 <- DimPlot(object = Seurat2, reduction = "umap", label = FALSE, cols=myCol3)
plot_grid(p3.1, p3.2)

#plot_grid(p1.1, p3.1)

dev.off()






