deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

load("parameters_hippo.R")

############################################################
# spot
############################################################
#ratio = sort( unique(c(0.5,seq(1,20,0.2),seq(1.5,19.5,1))) )
ratio = sort( unique( c(seq(0.5,10,0.5), seq(1,10,0.2)) ) )
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

names(result) = ratio
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
save(result, file="MouseBrainHippo_spot_BIC_new_cells_on_spot.RData")

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


################################################################################################
########## visualization
################################################################################################

leaf.plot = function(dist, xlim=NULL, ylim=NULL, upper=1:length(dist), lower=NULL, 
                     axis=TRUE, title=NULL, ref_title=NULL, main="leaf.plot")
{
    
    plot_dist <- function(d, col, size=0.5, title=NULL) {
        d = sort(d,decreasing=T)
        angles <- seq(0, pi, length.out = length(d))
        x = d*cos(angles)
        y = d*sin(angles)
        points(x, y, col=col, pch=16, cex=size)
        if( !is.null(title) )
        {
            i = which.max(abs(y))
            text(x[i], y[i], title)
        }
    }
    
    get_dist <- function(d) {
        d = sort(d,decreasing=T)
        angles <- seq(0, pi, length.out = length(d))
        x = d*cos(angles)
        y = d*sin(angles)
        data.frame(x,y)
    }
    
    all_dist = sort(do.call(c,dist))
    n = length(all_dist)
    
    hues = seq(15, 375, length = n + 1)
    colors <- grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
    all_cols = rev( colors )
    
    values = lapply(dist, get_dist)
    if(is.null(xlim)) xlim = range( do.call(c,lapply(values,function(z)z$x)) )
    if(is.null(ylim)) {
        ylim1 = range( do.call(c,lapply(values[upper],function(z) z$y )) )
        ylim2 = range( do.call(c,lapply(values[lower],function(z) -z$y )) )
        ylim = range(c(ylim1,ylim2))
    }
    #plot(NA, xlim=xlim, ylim=ylim, xlab="", ylab="")
    plot(NA, axes=FALSE, frame.plot=FALSE, xlim=xlim, ylim=ylim, xlab="",ylab="", main=main)
    tmp = sapply(1:length(dist), function(i) {
        d = sort(dist[[i]],decreasing=T)
        col = all_cols[ match(d, all_dist) ]
        if( i %in% upper ) plot_dist (d, col, s=0.5, title=title[i])
        if( i %in% lower ) plot_dist (-d, rev(col), s=0.5, title=title[i])
    } )
    
    if(axis) {
        points( all_dist, rep(0,length(all_dist)), col=all_cols, cex=1, pch=15 )
        x = seq(0, max(all_dist), length=10)
        x = x[-length(x)]
        y = rep(0,length(x))
        points( x, y, pch=3, cex=1 )
        arrows( 0, 0, max(all_dist), 0 )
    } else {
        points( all_dist, rep( min(ylim)-0.1,length(all_dist)), col=all_cols, cex=1, pch=15 )
    }
    
    points(0,0,pch=16,col=all_cols[3000],cex=3)
    if( !is.null(ref_title) ) text(0,0,ref_title)
    
}


barplot2 <- function(score)
{
    colors = c( "red", "blue", "green", "black", "cyan", "yellow")
    colors = col2rgb(colors)
    colors = apply(colors,2,function(x) rgb(x[1],x[2],x[3],max=255,alpha=0.7*255) )
    
    data <- do.call(cbind, lapply(score,function(x) colMeans(x$celltypes_on_spot)) )
    data = data[order(rownames(data)),]
    
    m <- sapply(score,function(x) mean(rowSums(x$celltypes_on_spot)) )
    se <- sapply(score,function(x) {
        y = x$celltypes_on_spot
        sd(rowSums(y))/sqrt(nrow(y))
    } )
    
    # Get the stacked barplot
    barx <- barplot(data, 
                    col=colors , 
                    border="white", 
                    space=0.04, 
                    font.axis=2, 
                    xlab="spot size",
                    las=3,
                    axes=FALSE,
                    ylim=c(-10,200))
    
    lines( barx, m )
    arrows(barx , m+se, barx , m, angle=90, code=3, length=0.05)
    
    labels = seq(0,150,50)
    axis(side=4, at=labels, labels=labels )
    
    ######################################################################
    #f = function(x) (x-min(x))/(max(x)-min(x))
    f = function(x) ( x-min(rmse)+min(rmse_se) ) / ( max(rmse)-min(rmse)+2*min(rmse_se) ) * ( max(m)-min(m) )
    
    rmse = sapply(score,function(x) mean(sqrt(x$mse)))
    rmse_se = sapply(score,function(x) {
        y = sqrt(x$mse)
        sd(y)/sqrt(length(y))
    })
    rmse2 = ( rmse-min(rmse) ) / ( max(rmse)-min(rmse) ) * ( max(m)-min(m) )
    rmse_se1 = ( rmse-rmse_se-min(rmse) ) / ( max(rmse)-min(rmse) ) * ( max(m)-min(m) )
    rmse_se2 = ( rmse+rmse_se-min(rmse) ) / ( max(rmse)-min(rmse) ) * ( max(m)-min(m) )
    
    rmse2 = f(rmse)
    rmse_se1 = f(rmse-rmse_se)
    rmse_se2 = f(rmse+rmse_se)
    
    points(barx,rmse2,pch=16, col='darkred')
    lines( barx, rmse2, col='darkred')
    arrows(barx , rmse_se1, barx , rmse_se2, angle=90, code=3, length=0.05, col='darkred')
    
    labels=seq(floor(min(rmse)),ceiling(min(rmse)),0.2)
    axis(side=2, at=f(labels), labels=labels, col='darkred' )
    
}


################################################################################################
########## load data
################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
load("MouseBrainHippo_spot_BIC_new_cells_on_spot.RData")

range = as.character( seq(0.5,7,0.5) )
result = result[range]

ratio=as.numeric(names(result))
score = lapply( result, function(x) calculate_BIC(x,ST_expr) )
rmse = sapply(score,function(x) mean(sqrt(x$mse)))
bic = sapply(score, function(x) mean(x$bic))
cell_count = sapply(score, function(x) mean(rowSums(x[[3]])) )

################################################################################################
############ draw barplot of RMSE
################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
pdf("spot_size_comparison.pdf",w=5,h=5)
barplot2(score)
dev.off()


################################################################################################
########## leaf plot
################################################################################################
prop_on_spot = lapply(score, function(x) {
    y = x$celltypes_on_spot
    t(apply(y,1,function(z) z/sum(z)))
} )

rmse = sapply(score,function(x) mean(sqrt(x$mse)))
ri = which.min(rmse)
ri = 6
ref = prop_on_spot[[ri]]

dist = lapply(prop_on_spot, function(x) {
    spots = intersect( rownames(x), rownames(ref) )
    x2 = x[match(spots,rownames(x)), ]
    ref2 = ref[match(spots, rownames(ref)), ]
    sapply( 1:length(spots), function(i) dist( rbind(x2[i,], ref2[i,]) ) )
    #sapply( 1:length(spots), function(i) cor(x2[i,], ref2[i,])^2 )
}  )


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hipocampus/results/STIE")
pdf('leaf_plot.pdf',w=6,h=6)
title = paste0(names(dist),"x")
leaf.plot(dist, upper=1:(ri-1), lower=(ri+1):length(dist), axis=TRUE, title=title, ref_title=title[ri])
leaf.plot(dist, upper=1:(ri-1), lower=(ri+1):length(dist), axis=TRUE)
leaf.plot(dist, upper=1:(ri-1), lower=(ri+1):length(dist), axis=FALSE)
dev.off()

