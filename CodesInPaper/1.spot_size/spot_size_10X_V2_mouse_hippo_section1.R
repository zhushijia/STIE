deconvolution = TRUE
source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/parameters/parameters_10X_V2_mouse_hippo_section1.R")
############################################################
# run on different spot size
############################################################
ratios = seq(0.5,5,0.5)
lambdas = c(0,1e3,1e6)

scores_las = list()
for(i in 1:length(lambdas)) 
{
    results = list()
    for(j in 1:length(ratios))
    {
        lai = lambdas[i]
        spotj = ratios[j]*spot_radius
        
        cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates=spot_coordinates, spot_radius=spotj)
        results[[j]]  <- STIE(ST_expr, Signature, cells_on_spot, 
                              features=c('shape'), 
                              lambda=lai, steps=30, 
                              known_signature=TRUE, known_cell_types=FALSE )
    }
    scores_las[[i]] = lapply( results, function(x) get_summary(x, ST_expr) )
    names(scores_las[[i]]) = names(results) = ratios
}
mses = lapply(scores_las, function(scores) lapply(scores,function(s) sqrt(s$mse) ))

############################################################
##### Visualization
############################################################
errplot <- function(X, barx = c(1:length(X)), 
                    ylab="", line=T, ylim=NULL, col='black')
{
    m <- sapply(X,function(x) mean(x) )
    se <- sapply(X,function(x) {
        sd(x)/sqrt(length(x))
    } )
    
    if( is.null(ylim) ) ylim=range( c( m-se, m+se) )
    
    if(line) {
        plot( barx, m, type="b", ylim=ylim, xlim=range(barx)+c(0,1), axes=FALSE, ylab=ylab, col='white', cex=0.5, xlab='' )
        axis(side=1,at=barx,label=names(X),las=3)
        axis(side=2)
    } 
    
    lines( barx, m, col=col, cex=0.5 )
    points( barx, m, pch=1, cex=1, col=1 )
    points( barx, m, pch=16, cex=1, col=col )
    arrows(barx , m+se, barx , m, angle=0, code=3, length=0.05, col=col )
    arrows(barx , m-se, barx , m, angle=0, code=3, length=0.05, col=col )
    
}

cols = get_my_colors(length(scores_las), mode=2)
barx = as.numeric(names(scores_las[[1]]))*2
errplot(mses[[1]],barx=barx, col=cols[1], ylab="RMSE")
lapply( 2:length(mses), function(i) errplot(mses[[i]],barx=barx+(i-1)*0.1, line=F, col=cols[i]) )
legend( "topright", legend=lambdas, col=cols, pch=16 )
