deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

load("parameters_hippo.R")

############################################################
# spot
############################################################
ratio = sort( unique(c(0.5,seq(1,20,0.2),seq(1.5,19.5,1))) )
#ratio = c(0.5, seq(1,6,1) )
result = list()
for(i in 1:length(ratio))
{
    cat(ratio[i],'\n')
    cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, ratio[i]*spot_radius)
    result[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                       morphology_steps=ceiling(steps/3), 
                       known_signature=TRUE, known_cell_types=FALSE)
}

score = lapply( result, function(x) BIC(x) )

bic = sapply(score, function(x) mean(x$bic))
mse = sapply(score, function(x) mean(x$mse))
cell_count = sapply(score, function(x) mean(rowSums(x[[2]])) )
names(result) = names(score) = names(bic) = names(cell_count) = ratio
plot(bic)
plot(cell_count)
ratio[which.min(bic)]

f = function(x) (x-min(x))/(max(x)-min(x))
plot( ratio, f(bic), type="l")
lines( ratio, f(cell_count) )

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
save(result, score, bic, cell_count, file="MouseBrainHippo_spot_BIC.RData")


############################################################
# lambda
############################################################
las = sort( unique(c( 10^(seq(-6,-1,1)), 0, 1:6 ) ) )
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)

result = list()
for(i in 1:length(las))
{
    cat(las[i],'\n')
    result[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=las[i], steps=30, 
                       morphology_steps=ceiling(steps/3), 
                       known_signature=TRUE, known_cell_types=FALSE)
}

score = lapply( result, function(x) calculate_BIC(x) )

bic = sapply(score, function(x) mean(x$bic))
cell_count = sapply(score, function(x) mean(rowSums(x[[3]])) )
names(result) = names(score) = names(bic) = names(cell_count) = las

plot(bic)
plot(cell_count)
las[which.min(bic)]

f = function(x) (x-min(x))/(max(x)-min(x))
plot( las[1:length(bic)], f(bic), type="l")
lines( las[1:length(bic)], f(cell_count) )

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
save(result, score, bic, cell_count, file="MouseBrainHippo_lambda_BIC.RData")

############################################################
# cluster on all genes
############################################################
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)

result = list()
for(i in 2:10)
{
    cat(i,'\n')
    
    cluster = read.csv(paste0("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/Visium_FFPE_Mouse_Brain/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ as.character(cluster[,1])%in%rownames(ST_expr), ]
    all( as.character(cluster[,1]) %in% rownames(ST_expr) )
    ST_expr3 = ST_expr[ match(as.character(cluster[,1]),rownames(ST_expr)), ]
    Signature = t(apply(ST_expr3, 2, function(x) tapply(x,cluster[,2],mean) ))
    
    res3[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                     morphology_steps=ceiling(steps/3), 
                     known_signature=FALSE, known_cell_types=FALSE)
}

score1 = lapply( res3[2:9], function(x) BIC(x) )
score2 = sapply(score1,mean)
names(score2) = names(score2) = rr
plot(score2)
rr[which.min(score2)]

