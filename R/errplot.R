#' errplot
#'
#' errplot plot lines with error bar
#'
#' @param X a list of numeric values to plot
#' @param barx the coordinates of points in the plot.
#' @param xlab a title for the x axis 
#' @param ylab a title for the y axis
#' @param line a boolean value indicating whether to plot lines
#' @param ylim range for the y axis
#' @param col character values indicating the color of the plot
#'
#' @return
#' @export
#'
#' @examples
errplot <- function(X, barx = c(1:length(X)), xlab="", ylab="", line=T, ylim=NULL, col='black')
{
    m <- sapply(X,function(x) mean(x) )
    se <- sapply(X,function(x) {
        sd(x)/sqrt(length(x))
    } )
    
    if( is.null(ylim) ) ylim=range( c( m-se, m+se) )
    
    if(line) {
        plot( barx, m, type="b", ylim=ylim, xlim=range(barx)+c(0,1), axes=FALSE, xlab=xlab, ylab=ylab, col='white', cex=0.5 )
        axis(side=1,at=barx,label=names(X),las=3)
        axis(side=2)
    } 
    
    lines( barx, m, col=col, cex=0.5 )
    points( barx, m, pch=1, cex=1, col=1 )
    points( barx, m, pch=16, cex=1, col=col )
    arrows(barx , m+se, barx , m, angle=0, code=3, length=0.05, col=col )
    arrows(barx , m-se, barx , m, angle=0, code=3, length=0.05, col=col )
    
}

