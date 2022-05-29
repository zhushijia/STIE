#' get_my_colors
#'
#' get_my_colors generates a vector of colors
#'
#' @param n an integer value representing the number of colors
#' @param mode an integer value from 1 to 3 representing the category of colors
#'
#' @return a vector of colors
#' @export
#' 
#' @author Shijia Zhu, \email{shijia.zhu@@UTsouthwestern.edu}
#'
#' @seealso \code{\link{plot_sub_image}}; 
#'
#'
get_my_colors <- function(n, mode=1)
{
    if(mode==1)
    {
        myCol <- c( "green", "red", "black", "cyan", "yellow", "darkorange", "purple", "steelblue",
                   "blue", "darkred", "grey", "blue4", "chartreuse4", "burlywood1", "darkgoldenrod4" )
        stopifnot( n<=length(myCol) )
        colors <- myCol[1:n]
    }
    
    if(mode==2) {
        
        colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                        "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", 
                        "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", 
                        "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                        "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", 
                        "#8DD3C7", "#999999")
        
        if (n <= length(colorSpace)) {
            colors <- colorSpace[1:n]
        } else {
            colors <- (grDevices::colorRampPalette(colorSpace))(n)
        }
        
    }

    if(mode==3) {
        hues = seq(15, 375, length = n + 1)
        colors <- grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    return(colors)
    
}