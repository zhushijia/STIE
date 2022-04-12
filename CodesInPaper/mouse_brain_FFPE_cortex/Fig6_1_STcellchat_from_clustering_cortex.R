source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_cortex/parameters_cortex.R")
library(NMF)
library(ggalluvial)

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE")
x = load("MouseBrainCortex_clustering.RData")

STIE_result = result[[8]]

cellchat <- get_STcellchat(STIE_result, ST_expr, database="mouse", db_category=NULL, max_reps=NULL )

######################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/CellChat")
pdf("STcellchat_from_clustering_selectK.pdf")
selectK(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")
dev.off()

######################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/CellChat")
pdf("STcellchat_from_clustering_summary.pdf", w=20, h=20)
par(mfrow = c(1,2), xpd=TRUE)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
#ht1 + ht2
ht1[1:20,] + ht2[1:20,]

plot.new()
nPatterns = 8
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", 
                                          k = nPatterns, width=4, height=30, font.size=7)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

plot.new()
plot.new()
nPatterns = 8
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", 
                                          k = nPatterns, width=4, height=30, font.size=7)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

dev.off()

######################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/CellChat")
save(cellchat, file="STcellchat_from_clustering.RData")


######################################################################################################
ref_1nm = spot_radius*2/55/args$x_scale

cell_dist = calculate_cell_dist2(cells_on_spot = STIE_result$cells_on_spot, 
                                 dist_type="boundary", 
                                 dist_cutoff=10*ref_1nm,
                                 axis="Minor" )

cell_dist$d_nm = cell_dist$d/ref_1nm

plot_STcellchat(cellchat, cell_dist, 
                im, image_transparency=0, 
                w=14000, h=7000, xoff=6000, yoff=8000, x_scale=0.1,
                plot_object="PTN", direction="outgoing", 
                color_use=NULL, lwd=10,
                plot_cell=FALSE, contour=cell_info$cell_contour)
    
netVisual_STcellchat(cellchat, cell_dist, STIE_result$cell_types)


