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
  
  #axis(side=1,at=barx,label=names(X),las=3)
  axis(side=1)
  axis(side=2)
  arrows(barx , m+se, barx , m, angle=90, code=3, length=0, col='grey' )
  arrows(barx , m-se, barx , m, angle=90, code=3, length=0, col='grey' )
  points( barx, m, pch=16, cex=0.5, col=col )
}



# sd_mr = 1; sd_er = 1
# sd_mr = 1; sd_er = 3
# sd_mr = 3; sd_er = 1
# sd_mr = 3; sd_er = 3

library(STIE)
library(quadprog)

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/Simulation/data")
fileName = paste0("SimulationData_sd_mr_",sd_mr, "_sd_er_", sd_er, ".RData")
x = load(fileName)

parentDir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/Simulation/results"
workdir=paste0(parentDir,"/SimulationData_sd_mr_",sd_mr, "_sd_er_", sd_er)
setwd(workdir)
x = load("STIE_result2.RData")
STIE_corrs = lapply(STIE_sums,function(x) x$prop_corr )
names(STIE_corrs) = c(0,10^c(1:6))

best_ind = which.max(sapply(STIE_corrs,median))
lai = as.numeric(names(STIE_corrs))[best_ind]

####################################################################################
## for additional feature simulation
####################################################################################
results = STIE_saveEachStep(ST_expr, Signature, cells_on_spot, features=features, 
                             lambda=lai, steps=100, morphology_steps=0,
                             known_signature=TRUE, known_cell_types=FALSE)

scores = lapply(results, function(result) get_summary(result,ST_expr) )
save(results, scores, file="STIE_result2_saveEachStep.RData")
##############################################################################

s = scores
par(mfrow=c(3,3))
plot(sapply(s,function(x) mean(x$logLik_Morp) ), type='l', main='morp')
plot(sapply(s,function(x) mean(x$logLik_Expr) ), type='l', main='expr')
plot(sapply(s,function(x) mean(x$L2sum) ), type='l', main='l2')
plot(sapply(s,function(x) x$rmse), type='l', main='mse')
plot(sapply(s,function(x) mean(x$logLik) ), type='l', main='logLik')
sd_mr
sd_er

mus = lapply(results,function(x)x$mu)
sds = lapply(results,function(x)x$sigma)

plotMus = function(mus) {
    d1 = nrow(mus[[1]])
    d2 = ncol(mus[[1]])
    cols = get_my_colors(d1,mode=2)
    plot(NA,xlim=c(1,length(mus)),ylim=range(do.call(c,mus)))
    for(i in 1:d1){
        for(j in 1:d2){
            x = sapply(mus,function(m)m[i,j])
            lines(x,col=cols[i],lty=j)
        }
    }
    
}

PE_on_spot_diff = lapply(2:length(results), function(i)  c(results[[i]]$PE_on_spot-results[[i-1]]$PE_on_spot)^2 )
mu_diff = lapply(2:length(results), function(i) c(results[[i]]$mu-results[[i-1]]$mu)^2 )
sd_diff = lapply(2:length(results), function(i) c(results[[i]]$sigma-results[[i-1]]$sigma)^2 )


pdf( paste0("STIE_result2_convergence_sd_mr_",sd_mr, "_sd_er_", sd_er,".pdf") )
par(mfrow=c(3,3))
plot( sapply(s,function(x) mean(x$logLik) ), type='l', main='logLik')
plot( sapply(PE_on_spot_diff,mean), type='l')
plot( sapply(mu_diff,mean), type='l')
plot( sapply(sd_diff,mean), type='l')
plotMus(mus)
plotMus(sds)

par(mfrow=c(3,2))
errplot( PE_on_spot_diff, ylab="PE_on_spot_diff", line=T, ylim=NULL, col='black')
errplot( mu_diff, ylab="mu_diff", line=T, ylim=NULL, col='black')
errplot( sd_diff, ylab="sd_diff", line=T, ylim=NULL, col='black')
plot(NA,xlim=c(1,10),ylim=c(1,10))
types = rownames(mus[[1]])
d1 = nrow(mus[[1]])
cols = get_my_colors(d1,mode=2)
legend('topleft',legend=types, lty=1, col=cols)
dev.off()


