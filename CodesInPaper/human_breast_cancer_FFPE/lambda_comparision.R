deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

load("parameter_BreastCancer.R")

############################################################
# spot
############################################################
ratio = sort( c(0.5,seq(1,3,0.2),seq(1.5,8,1),seq(4,8,1)) )
#ratio = seq(1,8,1)
result = list()
for(i in 1:length(ratio))
{
    cat(ratio[i],'\n')
    cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, ratio[i]*spot_radius)
    result[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=0, steps=30, 
                     morphology_steps=ceiling(steps/3), 
                     known_signature=TRUE, known_cell_types=FALSE)
}

score = lapply( result, function(x) calculate_BIC(x, ST_expr) )
mse = sapply(score, function(x) mean(x$mse))
bic = sapply(score, function(x) mean(x$bic))
cell_count = sapply(score, function(x) mean(rowSums(x[[3]])) )
names(result) = names(score) = names(bic) = names(cell_count) = ratio

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
save(result, score, ratio, mse, bic, cell_count, file="BreastCancer_spot_BIC_new_get_cells_on_spot.RData")


plot(bic)
plot(cell_count)
ratio[which.min(bic)]


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
#save(result, score, bic, cell_count, file="BreastCancer_spot_BIC.RData")

a = load("BreastCancer_spot_BIC.RData")
cell_count2 = do.call(cbind, lapply(score, function(x) colMeans(x[[3]]) ) )

f = function(x) (x-min(x))/(max(x)-min(x))
plot( ratio, f(bic), type="l")
lines( ratio, f(cell_count) )

plot(ratio, cell_count, type='l', ylim=c(0,max(cell_count)))
apply(cell_count2,1,function(x)lines(ratio,x))

############################################################
# lambda
############################################################
las = ( unique(c( 0, 1e-3, 1, 1e3, 1e6, 1e4, 1e5, 10, 100 )) )

result = list()
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)
for(i in 8:length(las))
{
    cat(las[i],'\n')
    result[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=las[i], steps=30, 
                       morphology_steps=ceiling(steps/3), 
                       known_signature=TRUE, known_cell_types=FALSE)
}

score = lapply( result, function(x) calculate_BIC(x,ST_expr) )

bic = sapply(score, function(x) mean(x$bic))
cell_count = sapply(score, function(x) mean(rowSums(x[[3]])) )
names(result) = names(score) = names(bic) = names(cell_count) = las

plot(bic)
plot(cell_count)
las[which.min(bic)]

f = function(x) (x-min(x))/(max(x)-min(x))
plot( las[1:length(bic)], f(bic), type="l")
lines( las[1:length(bic)], f(cell_count) )

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
save(result, score, bic, cell_count, file="BreastCancer_lambda_BIC.RData")


