module load R/4.1.1-gccmkl 

library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)

sce <- readVisium("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/Visium_FFPE_Mouse_Brain/outs")
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

ST.enhanced <- spatialEnhance(ST, q=Qn, platform="Visium", d=7,
                            model="t", gamma=2,
                            jitter_prior=0.3, jitter_scale=3.5,
                            nrep=1000, burn.in=100,
                            save.chain=TRUE)

################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/Signature/BroadInstitute_SingleCell")
Signature = read.delim("Major_cell_types_marker_genes.txt", header=T, row.names=1)
Signature = Signature[,!colnames(Signature)%in%c("Ependymal.cells")]

markers <- rownames(Signature)
ST.enhanced <- enhanceFeatures(ST.enhanced, ST,
                               feature_names=markers,
                               nrounds=0)

outputDir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/BayesSpace"
dir.create(outputDir)
setwd(outputDir)
save(ST.enhanced,file="ST_enhanced.RData")

