

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("HumanBreastCancer_clustering.RData")

i = 6
markers = find_sig_markers(STIE_result=result[[i]], ST_expr, transform="inverse.transform", DEG_pthres=1)
features = unique(markers$features)




workdir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq"
projectName = "Wu_etal_2021_BRCA_scRNASeq"
fileName = paste0(workdir,"/Seurat/",projectName,".Seurat.markers.txt") 
truth = read.delim( fileName, sep="\t", header=T)
truth = tapply(truth$gene, truth$cluster, as.character)

markers2 = subset(markers, pvalues<1e-40 & diff>0)
markers2 = tapply(markers2$features, markers2$clusters, as.character)

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














