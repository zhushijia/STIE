module load R/4.1.1-gccmkl 

library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)

sce <- readVisium("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/Visium_FFPE_Mouse_Brain/outs")
set.seed(102)
ST <- spatialPreprocess(sce, platform="ST", n.PCs=7, n.HVGs=2000, log.normalize=TRUE)
ST <- qTune(ST, qs=seq(2, 10), platform="ST", d=7)
#qPlot(ST)

Qn = 8
set.seed(149)
ST <- spatialCluster(ST, q=Qn, 
                   platform="ST", d=7,
                   init.method="mclust", model="t", gamma=2,
                   nrep=1000, burn.in=100,
                   save.chain=TRUE)

ST.enhanced <- spatialEnhance(ST, q=Qn, platform="ST", d=7,
                            model="t", gamma=2,
                            jitter_prior=0.3, jitter_scale=3.5,
                            nrep=1000, burn.in=100,
                            save.chain=TRUE)

################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/allen_cortex")
x = load("AllenCortex_MouseBrain_scRNASeq_DWLS_Signature.RData")
Signature = Signature[,order(colnames(Signature))]

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/SPOTlight")
load( "cortex_sc_cluster_markers_all.RData"  )

clusters = split( cluster_markers_all, cluster_markers_all$cluster)
markers = unique(do.call(c, lapply(clusters,function(x) {
    tmp = x[order(x$p_val),]
    as.character(x$gene)[1:min(100,nrow(tmp))]
} ) ))

expr = as.matrix( assays(ST)[[1]] )
markers <- intersect( markers, rownames(expr) )

ST.enhanced <- enhanceFeatures(ST.enhanced, ST,
                               feature_names=markers,
                               nrounds=0)

outputDir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/BayesSpace"
dir.create(outputDir)
setwd(outputDir)
save(ST.enhanced,file="ST_enhanced.RData")

################################################################################################
colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", 
                "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", 
                "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", 
                "#8DD3C7", "#999999")

BS_dir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/BayesSpace"
setwd(BS_dir)
x = load("ST_enhanced.RData")
bs_enhanced = ST.enhanced

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/allen_cortex")
x = load("AllenCortex_MouseBrain_scRNASeq_DWLS_Signature.RData")
Signature = Signature[,order(colnames(Signature))]

for(i in 2:10)
{
    cat(i,"\n")
    set.seed(149)
    bs_enhanced <- spatialCluster(bs_enhanced, q=i, 
                                  platform="ST", d=7,
                                  init.method="mclust", model="t", gamma=2,
                                  nrep=1000, burn.in=100,
                                  save.chain=TRUE)
    
    counts = logcounts(bs_enhanced)
    index = apply(counts, 1, function(x) any(is.na(x)) )
    counts = counts[!index,]
    counts = counts[ rownames(counts)%in%rownames(Signature), ]
    
    colData = data.frame(colData(bs_enhanced))
    spots = rownames(colData)
    cluster = colData$spatial.cluster
    names(cluster) = spots
    all( colnames(counts) == spots )
    
    bs_signature = t(apply( counts, 1, function(x) {
        tapply(x, cluster, function(y) mean(exp(y),na.rm=T) )   
    } ))
    
    outputDir = paste0(BS_dir,"/cluster/cluster",i)
    dir.create(outputDir, recursive=T)
    
    setwd(outputDir)
    pdf( paste0("cluster",i,".pdf"), 20, 20 )
        cols = colorSpace[1:i]
        clusterPlot(bs_enhanced, palette=cols, color="darkgrey")
    dev.off()
    
    if(i==6)
    {
        pdf( paste0("cluster",i,"_forFigure.pdf"), 20, 20 )
        cols = c("yellow", "blue", "red", "steelblue", "black", "cyan")
        clusterPlot(bs_enhanced, palette=cols)
        clusterPlot(bs_enhanced, color="black") +
            theme_bw() +
            xlab("Column") +
            ylab("Row") +
            labs(fill="BayesSpace\ncluster", title="Spatial clustering of ST_mel1_rep2")
        dev.off()
        
    }
    
}









