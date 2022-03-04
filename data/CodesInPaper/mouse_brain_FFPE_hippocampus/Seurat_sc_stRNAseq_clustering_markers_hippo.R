library(dplyr)
library(Seurat)
library(cowplot)
library(data.table)
library(Matrix)

buildSeurat <- function( counts, celltypes, projectName, outputDir )
{
    Seurat <- CreateSeuratObject(counts = counts, project = projectName, min.cells = 0)
    Seurat <- NormalizeData(object = Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    Seurat <- FindVariableFeatures(object = Seurat, selection.method = "vst", nfeatures = 2000)
    dim(Seurat$RNA@counts)
    dim(Seurat$RNA@data)
    
    # Run the standard workflow for visualization and clustering
    Seurat <- ScaleData(object = Seurat, verbose = FALSE)
    Seurat <- RunPCA(object = Seurat, npcs = 30, verbose = FALSE)
    Seurat <- RunUMAP(object = Seurat, reduction = "pca", dims = 1:30)
    Seurat <- FindNeighbors(object = Seurat, reduction = "pca", dims = 1:30)
    
    Idents(Seurat) = as.character(celltypes)
    Seurat.markers <- FindAllMarkers(object = Seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    
    dir.create(outputDir)
    setwd(outputDir)
    #saveRDS(Seurat, file = paste0(projectName,".data.rds"))
    write.table(Seurat.markers, paste0(projectName,".Seurat.markers.txt"), sep="\t", quote=F, col.names=T, row.names=F)
    
    pdf( paste0(projectName,".umap.pdf"),w=20,h=10)
    p2 <- DimPlot(object = Seurat, reduction = "umap", label = TRUE)
    plot_grid(p2)
    dev.off()
    
    pdf( paste0(projectName,".markers.pdf"),w=20,h=30)
    top10 <- Seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    DoHeatmap(object = Seurat, features = top10$gene) 
    dev.off()
}

##############################################################################
########## single cell
##############################################################################
workdir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell"
setwd(workdir)

count = data.frame(fread("DATA_MATRIX_LOG_TPM.txt"))
rownames(count) = as.character(count$GENE)
colnames(count) = gsub("[.]","_",colnames(count))
count = count[,-1]
#count = exp(count) - 1
annot = read.delim("CLUSTER_AND_SUBCLUSTER_INDEX.txt",sep='\t',skip=1,header=T)
annot$TYPE = gsub("-","_",as.character(annot$TYPE))
annot$TYPE = gsub("[.]","_",as.character(annot$TYPE))
colnames(annot)[2] = "subclass"
rownames(annot) = as.character(annot$TYPE)

all( colnames(count) %in% as.character(annot$TYPE) )
annot = annot[ match(colnames(count), as.character(annot$TYPE) ) , ]
all( colnames(count) == rownames(annot) )

index = which( as.character(annot$subclass)!="Ependymal" )
dataSC = count[,index]
celltypes = as.character(annot$subclass[index])

projectName = "BroadInstitute_MouseHippo_scRNASeq"
outputDir = paste0(workdir,"/Seurat")
buildSeurat( counts=dataSC, celltypes=celltypes, projectName=projectName, outputDir=outputDir )


##############################################################################
########## spatial transcriptome clustering
##############################################################################
workdir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/Kmeans/Seurat"

spaceranger_count_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/Visium_FFPE_Mouse_Brain"
matrix_paths <- paste( spaceranger_count_dir, "outs/filtered_feature_bc_matrix.h5", sep="/")
dataST <- as.matrix(Read10X_h5(matrix_paths))

for(i in 3:10)
{
    cat(i,'\n')
    cluster = read.csv(paste0(spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    all( colnames(dataST)==as.character(cluster[,1]) )
    cluster = cluster[ as.character(cluster[,1])%in%colnames(dataST), ]
    all( as.character(cluster[,1]) %in% colnames(dataST) )
    counts = dataST[ ,match(as.character(cluster[,1]),colnames(dataST)) ]
    
    projectName = paste0("cluster_",i)
    outputDir = paste0(workdir,"/",projectName)
    celltypes = as.character(cluster[,2])
    buildSeurat( counts, celltypes, projectName, outputDir )
}





