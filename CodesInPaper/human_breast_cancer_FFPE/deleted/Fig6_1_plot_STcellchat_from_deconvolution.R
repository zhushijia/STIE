source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")
library(NMF)
library(ggalluvial)


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
#x = load("HumanBreastCancer_clustering_2.5xSpot.RData")
x = load("HumanBreastCancer_clustering_2.5xSpot_lambda1000.RData")
stie_clustering = results

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
#x2 = load("BreastCancer_spot_new_get_cells_on_spot_full_signature.RData")
x2 = load("BreastCancer_lambda_comparison_2.5xSpot_fullSignature.RData")
stie_deconvolution = results

STIE_result = stie_deconvolution[['1000']]

cellchat <- get_STcellchat(STIE_result, ST_expr, database="human", db_category=NULL, max_reps=1, nboot=1, thres=2 )
cellchat <- get_STcellchat(STIE_result, ST_expr, database="human", db_category=NULL, max_reps=1000, nboot=100, thres=0.05 )
cellchat <- get_STcellchat(STIE_result, ST_expr, database="human", db_category=NULL, max_reps=NULL )
cellchat <- get_STcellchat(STIE_result, ST_expr, database="human", db_category=NULL, max_reps=1000, nboot=1, thres=2 )

######################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/CellChat")
pdf("STcellchat_from_deconvolution_selectK.pdf")
selectK(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")
dev.off()

######################################################################################################
pdf("STcellchat_from_deconvolution_summary.pdf", w=20, h=20)
par(mfrow = c(1,2), xpd=TRUE)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
#ht1 + ht2
ht1[1:20,] + ht2[1:20,]

plot.new()
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

plot.new()
plot.new()
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

dev.off()

######################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/CellChat")
save(cellchat, file="STcellchat_from_deconvolution.RData")


######################################################################################################






pattern = c("outgoing", "incoming")
patternSignaling <- methods::slot(cellchat, "netP")$pattern  [[pattern]]
data1 = patternSignaling$pattern$cell
data2 = patternSignaling$pattern$signaling
data = patternSignaling$data
data3 = merge(data1, data2, by.x = "Pattern", by.y = "Pattern")
data3$Contribution <- data3$Contribution.x * data3$Contribution.y
data3$Contribution <- data3$Contribution.x * data3$Contribution.y
#data3 <- data3[, colnames(data3) %in% c("CellGroup", "Signaling", "Contribution")]

count <- cellchat@net$count
weight <- cellchat@net$weight

pathway_names = cellchat@netP$pathways
pathway_probs = cellchat@netP$prob



