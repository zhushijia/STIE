deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

load("parameter_BreastCancer.R")

############################################################
# run on different spot size
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
    colors = c( "steelblue", "darkred", "black", "cyan", "yellow", "darkorange", "#4DAF4A")
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
    
    axis(side=4, at=seq(0,200,50), labels=seq(0,200,50) )
    
    ######################################################################
    #f = function(x) (x-min(x))/(max(x)-min(x))
    f = function(x) ( x-min(rmse) ) / ( max(rmse)-min(rmse) ) * ( max(m)-min(m) )
    
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
    
    labels=seq(12,18,1)
    axis(side=2, at=f(labels), labels=labels, col='darkred' )
    
}


################################################################################################
########## load data
################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/DWLS_on_Wu_etal_2021_BRCA_Signature")
x = load("Visium_FFPE_Human_Breast_Cancer_DWLS_on_Wu_etal_2021_BRCA_Signature.RData" )

prop_mat = prop_mat[, !colnames(prop_mat)%in%c("Myeloid","PVL")]
Signature = Signature[, match(colnames(prop_mat),colnames(Signature))]

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
z = load("BreastCancer_spot_BIC_new_get_cells_on_spot.RData")


################################################################################################
############ draw barplot of RMSE
################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
pdf("spot_size_comparison.pdf",w=5,h=5)
barplot2(score[ as.character(seq(0.5,7,0.5)) ])
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
ref = prop_on_spot[[ri]]

dist = lapply(prop_on_spot, function(x) {
    spots = intersect( rownames(x), rownames(ref) )
    x2 = x[match(spots,rownames(x)), ]
    ref2 = ref[match(spots, rownames(ref)), ]
    sapply( 1:length(spots), function(i) dist( rbind(x2[i,], ref2[i,]) ) )
    #sapply( 1:length(spots), function(i) cor(x2[i,], ref2[i,])^2 )
}  )

range = as.character( seq(0.5,7,0.5) )
dist = dist[range]

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
pdf('leaf_plot.pdf',w=6,h=6)
leaf.plot(dist, upper=1:5, lower=6:14, axis=TRUE, title=paste0(names(dist),"x"), ref_title="2.8x")
leaf.plot(dist, upper=1:5, lower=6:14, axis=TRUE)
leaf.plot(dist, upper=1:5, lower=6:14, axis=FALSE)
dev.off()

