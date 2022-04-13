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

######################################################################################################

cc = data.frame(celltypes = c("Plasmablasts","Bcells", "Tcells", "Myeloid", 
                              "Endothelial", "CAFs", "PVL", 
                              "CancerEpithelial"),
                colors = c( "#4DAF4A", "darkorange", "steelblue", "yellow", 
                            "cyan", "darkred", "purple",
                            "black" ) )
colors = as.character(cc$colors)[order(as.character(cc$celltypes))]

ref_1nm = spot_radius*2/55/args$x_scale

cell_dist = calculate_cell_dist(cells_on_spot = STIE_result$cells_on_spot, 
                                dist_type="boundary", 
                                dist_cutoff=3*ref_1nm,
                                axis="Major" )

cell_dist$d_nm = cell_dist$d/ref_1nm

######################################################################################################

cellchat <- get_STcellchat(STIE_result, ST_expr, database="human", db_category=NULL, max_reps=NULL )

nPatterns_out = 3
nPatterns_in = 7
lwd = 10
topN = 30
outputDir = paste0( "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/CellChat/deconvolution_nPatternOut",
                    nPatterns_out,"_nPatternIn_",nPatterns_in,"_lwd_",lwd)
dir.create(outputDir)
setwd(outputDir)
#save(cellchat, file="STcellchat_from_deconvolution.RData")

######################################################################################################
if(0)
{
    pdf("STcellchat_from_deconvolution_selectK.pdf", w=8, h=4)
    selectK(cellchat, pattern = "outgoing")
    selectK(cellchat, pattern = "incoming")
    dev.off()
}


######################################################################################################
pdf("identifyCommunicationPatterns.pdf", w=10, h=20)

plot.new()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns_out, color.use=colors)
plot.new()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns_out, color.use=colors, 
                                          width=3, height=30, font.size=4)
plot.new()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns_in, color.use=colors)
plot.new()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns_in, color.use=colors, 
                                          width=3, height=30, font.size=4)
dev.off()

patternSignaling <- methods::slot(cellchat, "netP")$pattern[["outgoing"]]
data = patternSignaling$data
pathway.show = colnames(data)[1:topN]
######################################################################################################
pdf("netVisual_STcellchat.pdf", w=6, h=10)
netVisual_STcellchat(cellchat, cell_dist, STIE_result$cell_types, color.use=colors)
dev.off()
######################################################################################################
pdf("netAnalysis_signalingRole_heatmap.pdf", w=10, h=5)
par(mfrow = c(1,2), xpd=TRUE)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, color.use=colors, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, color.use=colors, pattern = "incoming")
#ht1 + ht2
ht1[1:topN,] + ht2[1:topN,]
dev.off()
######################################################################################################
pdf("netAnalysis_river.pdf", w=10, h=20)
netAnalysis_river(cellchat, color.use=colors, pattern = "outgoing")
netAnalysis_river(cellchat, color.use=colors, pattern = "incoming")
dev.off()
######################################################################################################
pdf("netAnalysis_dot.pdf", w=6, h=2.5)
netAnalysis_dot(cellchat, color.use=colors, pattern = "outgoing", pathway.show=pathway.show, shape = 21)
netAnalysis_dot(cellchat, color.use=colors, pattern = "incoming", pathway.show=pathway.show, shape = 22)
dev.off()
######################################################################################################

for(i in 1:nPatterns_out)
{
    png( paste0("plot_STcellchat_outgoing_pattern_", i ,".png"), width=2000,height=2000,res=300)
    plot_STcellchat(cellchat, cell_dist, 
                    im, image_transparency=1, x_scale=args$x_scale,
                    plot_object=paste0("Pattern ",i), direction="outgoing", 
                    color_use=colors, lwd=lwd,
                    plot_cell=FALSE, contour=cell_info$cell_contour)
    dev.off()
}


for(i in 1:nPatterns_in)
{
    png( paste0("plot_STcellchat_incoming_pattern_", i ,".png"), width=2000,height=2000,res=300)
    plot_STcellchat(cellchat, cell_dist, 
                    im, image_transparency=1, x_scale=args$x_scale,
                    plot_object=paste0("Pattern ",i), direction="incoming", 
                    color_use=colors, lwd=lwd,
                    plot_cell=FALSE, contour=cell_info$cell_contour)
    dev.off()
}
    
######################################################################################################

for(i in 1:length(pathway.show))
{
    png( paste0("plot_STcellchat_pathway_", i, "_", pathway.show[i] ,"_outgoing.png"),width=2000,height=2000,res=300)
    plot_STcellchat(cellchat, cell_dist, 
                    im, image_transparency=1, x_scale=args$x_scale,
                    plot_object=pathway.show[i], direction="outgoing", 
                    color_use=colors, lwd=lwd,
                    plot_cell=FALSE, contour=cell_info$cell_contour)
    dev.off()
    
    png( paste0("plot_STcellchat_pathway_", i, "_", pathway.show[i] ,"_incoming.png"),width=2000,height=2000,res=300)
    plot_STcellchat(cellchat, cell_dist, 
                    im, image_transparency=1, x_scale=args$x_scale,
                    plot_object=pathway.show[i], direction="incoming", 
                    color_use=colors, lwd=lwd,
                    plot_cell=FALSE, contour=cell_info$cell_contour)
    dev.off()
}




