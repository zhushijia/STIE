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
cellchat <- get_STcellchat(STIE_result, ST_expr, database="human", db_category=NULL, max_reps=NULL )

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
cc = data.frame(celltypes = c("Plasmablasts","Bcells", "Tcells", "Myeloid", 
                              "Endothelial", "CAFs", "PVL", 
                              "CancerEpithelial"),
                colors = c( "#4DAF4A", "darkorange", "steelblue", "yellow", 
                            "cyan", "darkred", "purple",
                            "black" ) )
cc = subset(cc, celltypes %in% STIE_result$cell_types )
colors = as.character(cc$colors)[order(as.character(cc$celltypes))]

ref_1nm = spot_radius*2/55/args$x_scale

cell_dist = calculate_cell_dist2(cells_on_spot = STIE_result$cells_on_spot, 
                                 dist_type="boundary", 
                                 dist_cutoff=10*ref_1nm,
                                 axis="Minor" )

cell_dist$d_nm = cell_dist$d/ref_1nm

netVisual_STcellchat(cellchat, cell_dist, STIE_result$cell_types)

plot_STcellchat(cellchat, cell_dist, 
                im, image_transparency=0, 
                image_transparency=0, x_scale=0.1,
                plot_object="Pattern 1", direction="outgoing", 
                color_use=colors, lwd=10,
                plot_cell=FALSE, contour=cell_info$cell_contour)


plot_STcellchat(cellchat, cell_dist, 
                im, image_transparency=0, 
                image_transparency=0, w=3000, h=3000, xoff=8000, yoff=10000, x_scale=0.1,
                plot_object="Pattern 1", direction="outgoing", 
                color_use=colors, lwd=10,
                plot_cell=TRUE, contour=cell_info$cell_contour)

