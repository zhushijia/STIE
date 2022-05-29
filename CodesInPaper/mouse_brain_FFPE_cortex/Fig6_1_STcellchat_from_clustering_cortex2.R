source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_cortex/parameters_cortex.R")
library(NMF)
library(ggalluvial)

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE")
#x = load("MouseBrainCortex_clustering.RData")
x = load("MouseBrainCortex_clustering_2.5xSpot_lambda0.RData")

clusteri = 6
STIE_result = results[[clusteri]]

######################################################################################################

colors = NULL

ref_1nm = spot_radius*2/55/args$x_scale

cell_dist = calculate_cell_dist(cells_on_spot = STIE_result$cells_on_spot, 
                                dist_type="boundary", 
                                dist_cutoff=3*ref_1nm,
                                axis="Major" )

cell_dist$d_nm = cell_dist$d/ref_1nm

######################################################################################################

cellchat <- get_STcellchat(STIE_result, ST_expr, database="mouse", db_category=NULL, max_reps=NULL )

######################################################################################################

parentDir0="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/CellChat"
parentDir=paste0(parentDir0,"/Cluster",clusteri)
dir.create(parentDir, recursive=T)
setwd(parentDir)

if(0)
{
    pdf("STcellchat_from_deconvolution_selectK.pdf", w=8, h=4)
    selectK(cellchat, pattern = "outgoing")
    selectK(cellchat, pattern = "incoming")
    dev.off()
}

######################################################################################################

nPatterns_out = 2
nPatterns_in = 2
lwd = 10
topN = 30
prob_cutoff = 0
outputDir = paste0( parentDir, "/deconvolution_nPatternOut", nPatterns_out,"_nPatternIn_",nPatterns_in,
                    "_lwd_",lwd, "_prob_cutoff_", prob_cutoff)
dir.create(outputDir)
setwd(outputDir)
#save(cellchat, file="STcellchat_from_deconvolution.RData")

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
netVisual_STcellchat(cellchat, cell_dist, STIE_result$cell_types, color_use=colors)
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
                    plot_object=paste0("Pattern ",i), direction="outgoing", prob_cutoff=prob_cutoff,  
                    color_use=colors, lwd=lwd,
                    plot_cell=FALSE, contour=cell_info$cell_contour)
    dev.off()
}


for(i in 1:nPatterns_in)
{
    png( paste0("plot_STcellchat_incoming_pattern_", i ,".png"), width=2000,height=2000,res=300)
    plot_STcellchat(cellchat, cell_dist, 
                    im, image_transparency=1, x_scale=args$x_scale,
                    plot_object=paste0("Pattern ",i), direction="incoming", prob_cutoff=prob_cutoff,  
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
                    plot_object=pathway.show[i], direction="outgoing", prob_cutoff=prob_cutoff,  
                    color_use=colors, lwd=lwd,
                    plot_cell=FALSE, contour=cell_info$cell_contour)
    dev.off()
    
    png( paste0("plot_STcellchat_pathway_", i, "_", pathway.show[i] ,"_incoming.png"),width=2000,height=2000,res=300)
    plot_STcellchat(cellchat, cell_dist, 
                    im, image_transparency=1, x_scale=args$x_scale,
                    plot_object=pathway.show[i], direction="incoming", prob_cutoff=prob_cutoff,  
                    color_use=colors, lwd=lwd,
                    plot_cell=FALSE, contour=cell_info$cell_contour)
    dev.off()
}


######################################################################################################
pdf("COLLAGEN_example.pdf")
plot_STcellchat(cellchat, cell_dist, 
                im, image_transparency=0.8, 
                w=3000, h=3000, xoff=8000, yoff=10000, x_scale=0.2,
                plot_object="COLLAGEN", direction="outgoing", prob_cutoff=prob_cutoff,  
                color_use=colors, lwd=lwd,
                plot_cell=FALSE, contour=cell_info$cell_contour)

plot_STcellchat(cellchat, cell_dist, 
                im, image_transparency=0.8, 
                w=3000, h=3000, xoff=8000, yoff=10000, x_scale=0.2,
                plot_object="COLLAGEN", direction="incoming", prob_cutoff=prob_cutoff,  
                color_use=colors, lwd=lwd,
                plot_cell=FALSE, contour=cell_info$cell_contour)
dev.off()













