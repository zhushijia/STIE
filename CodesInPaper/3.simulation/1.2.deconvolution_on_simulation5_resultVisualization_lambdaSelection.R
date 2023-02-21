errplot <- function(X, ylab="", line=T, ylim=NULL, col='black')
{
    nonNull = which(sapply(X,function(x)!is.null(x)))
    barx = c(1:length(X))[nonNull]
    X = X[nonNull]
    
    m <- sapply(X,function(x) mean(x) )
    se <- sapply(X,function(x) {
        sd(x)/sqrt(length(x))
    } )
    
    if( is.null(ylim) ) ylim=range( c( m-se, m+se) )
    
    if(line) {
        plot( barx, m, type="b", ylim=ylim, axes=FALSE, ylab=ylab, col='black', cex=0.5 )
    } else {
        plot( barx, m, type="l", ylim=ylim, axes=FALSE, ylab=ylab, col='white' )
    }
    
    axis(side=1,at=barx,label=names(X),las=3)
    axis(side=2)
    arrows(barx , m+se, barx , m, angle=90, code=3, length=0.05, col='grey' )
    arrows(barx , m-se, barx , m, angle=90, code=3, length=0.05, col='grey' )
    points( barx, m, pch=16, cex=1.5, col=col )
}


background = function(x) {
  if(class(x)=='list') x = do.call(c,x)
  r = range(x)
  d = r[2]-r[1]
  rect(0.5, r[1]-d, 2.5, r[2]+d, col=adjustcolor("red", 0.1),lty=0)
  rect(2.5, r[1]-d, 9.5, r[2]+d, col=adjustcolor("blue", 0.1),lty=0)
  #rect(9.5, r[1]-d, 14.5, r[2]+d, col=adjustcolor("orange", 0.1),lty=0)
}

##############################################################################
# sd_mr = 1; sd_er = 1
# sd_mr = 1; sd_er = 3
# sd_mr = 3; sd_er = 1
# sd_mr = 3; sd_er = 3

library(STIE)
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/Simulation/data")
fileName = paste0("SimulationData_sd_mr_",sd_mr, "_sd_er_", sd_er, ".RData")
x = load(fileName)
parentDir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/Simulation/results"
workdir=paste0(parentDir,"/SimulationData_sd_mr_",sd_mr, "_sd_er_", sd_er)
setwd(workdir)
x=load("STIE_result2.RData")
s = lapply(STIE_results,get_summary,ST_expr)
##############################################################################
STIE_corrs2 = c( list( corr_morph=corr_morph, corr_expr=corr_expr  ), 
                 STIE_corrs )
rmse <- lapply(s, function(x) sqrt(x$mse))
rmse2 <- c(list(mse_expr=sqrt(mse_expr)) , rmse)
cols = c( rep("grey75",2), rep('white',length(STIE_corrs)), "darkred","darkorange","darkgreen","steelblue","purple","magenta")


pdf( paste0("STIE_result2_sd_mr_",sd_mr, "_sd_er_", sd_er,".pdf"), w=10, h=15 )
#pdf( paste0("STIE_result2_sd_mr_",sd_mr, "_sd_er_", sd_er,".pdf") )
par(mfrow=c(3,3))
boxplot( STIE_corrs2, ylim=c(-1,1), outline=F, las=3 )
background(STIE_corrs2)
boxplot( STIE_corrs2, ylim=c(-1,1), outline=F, las=3, 
         col=cols[1:length(STIE_corrs2)], add=T )

#errplot( STIE_corrs, ylab='corr', ylim=c(-1,1))
errplot( lapply(s, function(x) sqrt(x$mse)), ylab="RMSE", col='darkred')
errplot( lapply(s, function(x) x$logLik_Morp ), ylab="logLik_Morp", col='darkgreen')
plot( sapply(s, function(x) mean(x$L2sum) ), ylab='L2', type='b', axes=FALSE, pch=16, cex=0.1 )
points( sapply(s, function(x) mean(x$L2sum) ), pch=16, col='darkorange', cex=1.5 )
axis(side=1)
axis(side=2)
sd_mr
sd_er
dev.off()


sapply(STIE_corrs2, mean)
sapply(STIE_corrs2, median)






