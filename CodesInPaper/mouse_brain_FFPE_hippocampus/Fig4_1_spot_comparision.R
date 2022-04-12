deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo3.R")

############################################################
# spot
############################################################
ratio = seq(0.5,10,0.5)
results = list()
for(i in 1:length(ratio))
{
    cat(ratio[i],'\n')
    cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, ratio[i]*spot_radius)
    results[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                       known_signature=TRUE, known_cell_types=FALSE)
}

score = lapply( results, function(x) calculate_BIC(x, ST_expr) )
names(results) = names(score) = ratio

outputDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE"
dir.create(outputDir, recursive=T)
setwd(outputDir)
save(results, score, file="MouseBrainHippo_spots.RData")

################################################################################################
########## visualization
################################################################################################
outputDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE"
setwd(outputDir)
load("MouseBrainHippo_spots.RData")

ind = as.character(seq(0.5,7,0.5))
results = results[ind]
score = score[ind]
PME_diff = get_PME_diff(results) 

################################################################################################
############ draw barplot of RMSE
################################################################################################
pdf("hippo_spot_size_comparison.pdf")
colors = c( "magenta", "blue", "green", "black", "orange")
celltypes = c("CA1", "CA2", "CA3", "DG", "Glia")
par(mfrow=c(2,2))
errplot( PME_diff, "PME_diff" )
errplot( lapply(score, function(x) sqrt(x$mse) ), "RMSE" )
errplot( lapply(score, function(x) x$bic ), "BIC" )
barplot1(score, colors, celltypes)
dev.off()


################################################################################################
########## leaf plot
################################################################################################
prop_on_spot = lapply(score, function(x) {
    y = x$celltypes_on_spot
    t(apply(y,1,function(z) z/sum(z)))
} )

rmse = sapply(score,function(x) mean(sqrt(x$mse)))
ri = 5
ref = prop_on_spot[[ri]]

dist = lapply(prop_on_spot, function(x) {
    spots = intersect( rownames(x), rownames(ref) )
    x2 = x[match(spots,rownames(x)), ]
    ref2 = ref[match(spots, rownames(ref)), ]
    sapply( 1:length(spots), function(i) dist( rbind(x2[i,], ref2[i,]) ) )
    #sapply( 1:length(spots), function(i) cor(x2[i,], ref2[i,])^2 )
}  )

pdf('hippo_leaf_plot.pdf',w=6,h=6)
title = paste0(names(dist),"x")
leaf.plot(dist, upper=1:(ri-1), lower=(ri+1):length(dist), axis=TRUE, title=title, ref_title=title[ri])
leaf.plot(dist, upper=1:(ri-1), lower=(ri+1):length(dist), axis=TRUE)
leaf.plot(dist, upper=1:(ri-1), lower=(ri+1):length(dist), axis=FALSE)
dev.off()

