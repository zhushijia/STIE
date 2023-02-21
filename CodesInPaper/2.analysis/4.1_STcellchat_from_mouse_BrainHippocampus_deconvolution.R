library(NMF)
library(ggalluvial)

deconvolution = TRUE
source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/parameters/parameters_10X_Visium_mouse_brain_hippo3.R")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)
############################################################
## run deconvolution
############################################################
STIE_result = STIE(ST_expr, Signature, cells_on_spot, features=c('shape','size'),
              steps=30, known_signature=TRUE, known_cell_types=FALSE, 
              lambda=0 )
######################################################################################################

colors = c('magenta','blue','green','black','orange')

######################################################################################################
##### calculate the cell-cell distance
######################################################################################################

ref_1nm = spot_radius*2/55/x_scale

cell_dist = calculate_cell_dist(cells_on_spot = STIE_result$cells_on_spot, 
                                dist_type="boundary", 
                                dist_cutoff=3*ref_1nm,
                                axis="Major" )

cell_dist$d_nm = cell_dist$d/ref_1nm

######################################################################################################
##### run CellChat
######################################################################################################


cellchat <- get_STcellchat(STIE_result, ST_expr, database="mouse", db_category=NULL, max_reps=NULL )

######################################################################################################
##### select patterns
######################################################################################################

pdf("STcellchat_from_deconvolution_selectK.pdf", w=8, h=4)
selectK(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")
dev.off()

######################################################################################################
##### set pattern number
######################################################################################################

nPatterns_out = 5
nPatterns_in = 5
lwd = 1
topN = 30
prob_cutoff = 0

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
##### visualize the outgoing and incoming patterns
######################################################################################################
for(i in 1:nPatterns_out)
{
    png( paste0("plot_STcellchat_outgoing_pattern_", i ,".png"), width=2000,height=2000,res=300)
    plot_STcellchat(cellchat, cell_dist, 
                    im, image_transparency=1, x_scale=x_scale,
                    plot_object=paste0("Pattern ",i), direction="outgoing", prob_cutoff=prob_cutoff,  
                    color_use=colors, lwd=lwd,
                    plot_cell=FALSE, contour=cell_info$cell_contour)
    dev.off()
}


for(i in 1:nPatterns_in)
{
    png( paste0("plot_STcellchat_incoming_pattern_", i ,".png"), width=2000,height=2000,res=300)
    plot_STcellchat(cellchat, cell_dist, 
                    im, image_transparency=1, x_scale=x_scale,
                    plot_object=paste0("Pattern ",i), direction="incoming", prob_cutoff=prob_cutoff,  
                    color_use=colors, lwd=lwd,
                    plot_cell=FALSE, contour=cell_info$cell_contour)
    dev.off()
}

######################################################################################################
##### visualize top N outgoing and incoming pathways
######################################################################################################
for(i in 1:length(pathway.show))
{
    png( paste0("plot_STcellchat_pathway_", i, "_", pathway.show[i] ,"_outgoing.png"),width=2000,height=2000,res=300)
    plot_STcellchat(cellchat, cell_dist, 
                    im, image_transparency=1, x_scale=x_scale,
                    plot_object=pathway.show[i], direction="outgoing", prob_cutoff=prob_cutoff,  
                    color_use=colors, lwd=lwd,
                    plot_cell=FALSE, contour=cell_info$cell_contour)
    dev.off()
    
    png( paste0("plot_STcellchat_pathway_", i, "_", pathway.show[i] ,"_incoming.png"),width=2000,height=2000,res=300)
    plot_STcellchat(cellchat, cell_dist, 
                    im, image_transparency=1, x_scale=x_scale,
                    plot_object=pathway.show[i], direction="incoming", prob_cutoff=prob_cutoff,  
                    color_use=colors, lwd=lwd,
                    plot_cell=FALSE, contour=cell_info$cell_contour)
    dev.off()
}



