module load R/4.1.1-gccmkl 

library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)

sce <- readVisium("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/outs")
set.seed(102)
breastcancer <- spatialPreprocess(sce, platform="ST", n.PCs=7, n.HVGs=2000, log.normalize=TRUE)
breastcancer <- qTune(breastcancer, qs=seq(2, 10), platform="ST", d=7)
qPlot(breastcancer)

Qn = 8
set.seed(149)
breastcancer <- spatialCluster(breastcancer, q=Qn, 
                               platform="ST", d=7,
                               init.method="mclust", model="t", gamma=2,
                               nrep=1000, burn.in=100,
                               save.chain=TRUE)

breastcancer.enhanced <- spatialEnhance(breastcancer, q=Qn, platform="ST", d=7,
                                    model="t", gamma=2,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=1000, burn.in=100,
                                    save.chain=TRUE)

################################################################################################
####### Signature
################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
Signature = Signature[,order(colnames(Signature))]

expr = as.matrix( assays(breastcancer)[[1]] )
markers <- intersect( rownames(Signature), rownames(expr) )

if(0)
{
    breastcancer.enhanced <- enhanceFeatures(breastcancer.enhanced, breastcancer,
                                             feature_names=markers,
                                             nrounds=0)
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/BayesSpace")
    save(breastcancer.enhanced,file="breastcancer_enhanced_onlySignatureGenes.RData")
}

if(0)
{
    breastcancer.enhanced <- enhanceFeatures(breastcancer.enhanced, breastcancer,
                                             nrounds=0)
    
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/BayesSpace")
    save(breastcancer.enhanced,file="breastcancer_enhanced_allGenes.RData")
}


################################################################################################
####### output cluster
################################################################################################
colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", 
                "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", 
                "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", 
                "#8DD3C7", "#999999")

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/BayesSpace")
x = load("breastcancer_enhanced_allGenes.RData")
bs_enhanced = breastcancer.enhanced

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
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
    
    parentDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/BayesSpace/cluster"
    outputDir = paste0(parentDir,"/cluster",i)
    dir.create(outputDir)
    setwd(outputDir)
    save( counts, cluster, file=paste0("cluster",i,".RData") )
    
    pdf( paste0("cluster",i,".pdf"), 20, 20 )
    
    if(1) {
        cols = colorSpace[1:i]
        clusterPlot(bs_enhanced, palette=cols, color="darkgrey")
    } else {
        clusterPlot(bs_enhanced, color="black") +
            theme_bw() +
            xlab("Column") +
            ylab("Row") +
            labs(fill="BayesSpace\ncluster", title="Spatial clustering")
    }
    dev.off()
    
    if(i==6) {
        pdf( paste0("cluster",i,"_for_figure.pdf"), w=20, h=20 )
        #png( paste0("cluster",i,"_for_figure.png"), w=3000, h=3000, res=300 )
        cc = data.frame( cluster = c(5,4,3,1,2,6), 
                         colors = c("steelblue",'darkgreen',"darkorange",'yellow','cyan','darkred') )
        cols = as.character(cc$colors)[ order(cc$cluster) ]
        #clusterPlot(bs_enhanced, palette=cols)
        clusterPlot(bs_enhanced, palette=cols, color="darkgrey")
        #clusterPlot(bs_enhanced, palette=cols, color="black")
        dev.off()
    }
    if(i==7) {
        pdf( paste0("cluster",i,"_for_figure.pdf"), w=20, h=20 )
        #png( paste0("cluster",i,"_for_figure.png"), w=3000, h=3000, res=300 )
        cc = data.frame( cluster = c(6,4,3,1,2,5,7), 
                         colors = c("steelblue",'darkgreen',"darkorange",'yellow','cyan','darkred','black') )
        cols = as.character(cc$colors)[ order(cc$cluster) ]
        clusterPlot(bs_enhanced, palette=cols, color="darkgrey")
        dev.off()
    }
    
}









