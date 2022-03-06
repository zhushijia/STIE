setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/DWLS_on_Wu_etal_2021_BRCA_Signature")
x = load("Visium_FFPE_Human_Breast_Cancer_DWLS_on_Wu_etal_2021_BRCA_Signature.RData" )

prop_mat = prop_mat[, !colnames(prop_mat)%in%c("Myeloid","PVL")]
Signature = Signature[, match(colnames(prop_mat),colnames(Signature))]


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
z = load("BreastCancer_spot_BIC_new_get_cells_on_spot.RData")

rmse = sapply(score,function(x) mean(sqrt(x$mse)))
f = function(x) (x-min(x))/(max(x)-min(x))
plot( ratio, f(mse), type="l")
lines( ratio, f(cell_count) )


method = 'pearson'
c1 = cor(prop_mat, method=method)
c2 = cor(score[['6']]$celltypes_on_spot, method=method)
cs = cor(Signature, method=method)

index = upper.tri(c1)
wilcox.test( c(c1[index])^2, c(c2[index])^2 , 'greater', paired=T)
t.test( c(c1[index])^2, c(c2[index])^2 , 'greater', paired=T)

c1^2 - c2^2

par(mfrow=c(2,2))
plot( c(c2[index])^2, c(c1[index])^2)
abline(b=1,a=0)

boxplot( c(c1[index])^2 , c(c2[index])^2 )

plot( c1[index], cs[index] )
plot( c2[index], cs[index] )

# to test linear relationship
# both are associated with signature correlation, but c1 is driven by it,  
# while c2 is not
cor.test( c1[index], cs[index], method='pearson' )
cor.test( c2[index], cs[index], method='pearson' )


################################################################################################
prop_on_spot = lapply(score, function(x) {
    y = x$celltypes_on_spot
    t(apply(y,1,function(z) z/sum(z)))
} )

ri = which.min(rmse)
ref = prop_on_spot[[ri]]
#ref = prop_on_spot[[1]]

spots = rownames(prop_on_spot[[1]])
for(i in 1:length(prop_on_spot)) spots = intersect( spots, rownames(prop_on_spot[[i]]) )

dist = lapply(prop_on_spot, function(x) {
    spots = intersect( rownames(x), rownames(ref) )
    x2 = x[match(spots,rownames(x)), ]
    ref2 = ref[match(spots, rownames(ref)), ]
    sapply( 1:length(spots), function(i) dist( rbind(x2[i,], ref2[i,]) ) )
    #sapply( 1:length(spots), function(i) cor(x2[i,], ref2[i,])^2 )
}  )

dist2 = dist
dist2[ ri:length(dist) ] = lapply( dist[ ri:length(dist) ], function(x) -x )

plot(ratio, sapply(dist2,median), type='l')
plot(ratio, sapply(dist2,mean), type='l')

max( do.call(c,dist) )

plot_dist <- function(d, c='grey') {
    d = sort(d)
    i = seq( 1, (length(d)-1), 2 )
    j = rev(seq( 2, length(d), 2 ))
    d = d[ intersect( c(i,j), 1:length(d) ) ]
    angles <- seq(0, 2*pi, length.out = length(d))
    lines(d*cos(angles), d*sin(angles), col=c, lwd=1)
}

plot_dist <- function(d, c, s=0.5, title=NULL) {
    d = sort(d,decreasing=T)
    angles <- seq(0, pi, length.out = length(d))
    x = d*cos(angles)
    y = d*sin(angles)
    points(x, y, col=c, pch=16, cex=s)
    #data.frame(x=x,y=y)
    if( !is.null(title) )
    {
        i = which.max(abs(y))
        text(x[i],y[i],title)
    }
}

index = as.character( sort(c(2.8,seq(0.5,7,0.5))) )
dist3 = dist[index]
plot(NA,xlim=c(-0.8,1.5),ylim=c(-0.6,0.6))

all_dist = sort(do.call(c,dist3))
n = length(all_dist)
all_cols = rev( get_my_colors(n, 3) )
sapply(1:length(dist3), function(i) {
    d = sort(dist3[[i]],decreasing=T)
    col = all_cols[ match(d, all_dist) ]
    title = paste0(index[i],"X")
    if(as.numeric(index[i])<2.8) plot_dist (d, col, s=0.5, title=title)
    #if(as.numeric(index[i])==2.8) plot_dist (d, c, s=2)
    if(as.numeric(index[i])>2.8) plot_dist (-d, rev(col), s=0.5, title=title)
} )

points( all_dist, rep(0,length(all_dist)), col=all_cols, cex=1, pch=15 )

x = seq(0, max(all_dist), length=10)
x = x[-length(x)]
y = rep(0,length(x))
points( x, y, pch=3, cex=1 )
#points( -all_dist, rep(0,length(all_dist)), col=all_cols, cex=all_dist*2, pch=17 )
arrows( 0, 0, max(all_dist), 0 )



