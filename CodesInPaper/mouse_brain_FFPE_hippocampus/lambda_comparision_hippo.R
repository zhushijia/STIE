deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo.R")

############################################################
# lambda
############################################################
las = sort( unique(c( 0, 10^c(-3:6) )) )

result = list()
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2*spot_radius)

for(i in 1:length(las))
{
    cat(las[i],'\n')
    result[[i]] = STIE(ST_expr, Signature, cells_on_spot, features, lambda=las[i], steps=30, 
                       known_signature=TRUE, known_cell_types=FALSE)
}

names(result) = las

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
save(result, file="MouseBrainHippo_lambda.RData")

########################################################################################################################
########################################################################################################################
########################################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
load("MouseBrainHippo_lambda.RData")
score = lapply( result, function(x) calculate_BIC(x,ST_expr) )

PME_diff = lapply( result, function(res) {
    PE_on_spot = res$PE_on_spot
    PM_on_spot = apply(res$PM_on_cell, 2, function(x) tapply(x, as.character(res$cells_on_spot$spot), mean) )
    PM_on_spot = t( apply(PM_on_spot,1,function(x)x/sum(x)) )
    PM_on_spot = PM_on_spot[ match( rownames(PE_on_spot), rownames(PM_on_spot) ) , ]
    rowSums((PE_on_spot-PM_on_spot)^2)
} )

errplot <- function(X)
{
    m <- sapply(X,function(x) mean(x) )
    se <- sapply(X,function(x) {
        sd(x)/sqrt(length(x))
    } )
    
    barx = 1:length(m)
    plot( m, type="l", ylim=range( c( m-se, m+se) ) )
    axis(side=1,at=barx,label=names(X),las=3)
    points( 1:length(m), m, pch=16, cex=1.5 )
    arrows(barx , m+se, barx , m, angle=90, code=3, length=0.05)
    arrows(barx , m-se, barx , m, angle=90, code=3, length=0.05)
}

par(mfrow=c(2,2))
errplot( PME_diff )
errplot( lapply(score, function(x) sqrt(x$mse) ) )

rmse = sapply( score, function(x) mean(sqrt(x$mse)) )

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
pdf("MouseBrainHippo_lambda_comparision.pdf",w=4,h=4)
errplot(score)
dev.off()





