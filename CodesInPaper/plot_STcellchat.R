cellchat <- get_cellchat(STIE_result, ST_expr, database="human", db_category=NULL, max_reps=NULL )

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
ht1[1:20,] + ht2[1:20,]


library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")


selectK(cellchat, pattern = "incoming")
nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")


cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)







pattern = c("outgoing", "incoming")
patternSignaling <- methods::slot(cellchat, "netP")$pattern[[pattern]]
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



