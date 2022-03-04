deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

load("parameters_hippo.R")

############################################################
# spot
############################################################
#ratio = sort( unique(seq(0.5,8,0.5)) )
ratio = sort( unique( c(seq(0.5,10,0.5), seq(1,10,0.2)) ) )
result = list()
for(i in 1:length(ratio))
{
    cat(ratio[i],'\n')
    cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, ratio[i]*spot_radius)
    result[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                       morphology_steps=ceiling(steps/3), 
                       known_signature=TRUE, known_cell_types=FALSE)
}

names(result) = ratio
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE")
save(result, score, bic, cell_count, file="MouseBrainCortex_spot_BIC_new_cells_on_spot.RData")


score = lapply( result, function(x) BIC(x) )

bic = sapply(score, function(x) mean(x$bic))
cell_count = sapply(score, function(x) mean(rowSums(x[[2]])) )
names(result) = names(score) = names(bic) = names(cell_count) = ratio
plot(bic)
plot(cell_count)
ratio[which.min(bic)]

f = function(x) (x-min(x))/(max(x)-min(x))
plot( ratio, f(bic), type="l")
lines( ratio, f(cell_count) )

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE")
save(result, score, bic, cell_count, file="MouseBrainCortex_spot_BIC.RData")


############################################################
# lambda
############################################################
result = list()
las = sort( c( 10^(seq(-6,-1,1)), 0, 1:10 ) )
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.4*spot_radius)

for(i in 15:length(las))
{
    cat(las[i],'\n')
    result[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=las[i], steps=30, 
                       morphology_steps=ceiling(steps/3), 
                       known_signature=TRUE, known_cell_types=FALSE)
}

score = lapply( result, function(x) BIC(x) )

bic = sapply(score, function(x) mean(x$bic))
cell_count = sapply(score, function(x) mean(rowSums(x[[2]])) )
names(result) = names(score) = names(bic) = names(cell_count) = ratio

plot(bic)
plot(cell_count)
las[which.min(bic)]

f = function(x) (x-min(x))/(max(x)-min(x))
plot( las[1:length(bic)], f(bic), type="l")
lines( las[1:length(bic)], f(cell_count) )

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/STIE")
save(result, score, bic, cell_count, file="MouseBrainCortex_lambda_BIC.RData")

