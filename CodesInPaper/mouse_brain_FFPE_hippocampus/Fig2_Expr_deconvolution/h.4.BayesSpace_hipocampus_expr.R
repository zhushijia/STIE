library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
library(hdf5r)

####################################################################################
############## BayesSpace enhanced transcriptome
####################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/BayesSpace")
x = load("ST_enhanced.RData")

colData = data.frame(colData(ST.enhanced))
counts = logcounts(ST.enhanced)
index=apply(counts, 1, function(x) any(is.na(x)) )
counts_2 = counts[!index,]
sig_genes = rownames(counts_2)
data = data.frame( colData, expr=colSums(counts_2) )

pdf("ST_enhanced_expr.pdf")
sp2 <- ggplot(data, aes(x=row, y=col, color=expr)) + geom_point(size = 3) +
       theme_bw() +
       theme(panel.grid = element_blank())
#sp2+scale_color_gradient(low="blue", high="red")
mid<-quantile(data$expr,0.5)
sp2+scale_color_gradient2(midpoint=mid, 
                          low="blue", mid="white",
                          high="red", space ="Lab" )
dev.off()

####################################################################################
##############  BayesSpace original
####################################################################################
colData = data.frame(colData(ST))
counts = logcounts(ST)
counts_2 = counts[rownames(counts)%in%sig_genes,]
data = data.frame( colData, expr=colSums(counts_2) )

pdf("ST_original_expr.pdf",w=3.8,h=4)
sp2 <- ggplot(data, aes(x=col, y=row, color=expr)) + geom_point(size = 3) +
    theme_bw() +
    theme(panel.grid = element_blank())
#sp2+scale_color_gradient(low="blue", high="red")
mid<-quantile(data$expr,0.5)
sp2+scale_color_gradient2(midpoint=mid, 
                          low="blue", mid="white",
                          high="red", space ="Lab" )
dev.off()

