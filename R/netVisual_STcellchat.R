#' netVisual_STcellchat
#'
#' @param cellchat 
#' @param cell_dist 
#' @param cell_types 
#'
#' @return
#' @export
#'
#' @examples
netVisual_STcellchat = function(cellchat, cell_dist, cell_types)
{
    
    weight = cellchat@net$weight
    
    cell_pair = weight*0
    
    for(k in 1:nrow(cell_dist))
    {
        cat(k,"\n")
        i = cell_dist$i[k]
        j = cell_dist$j[k]
        ci = cell_types[i]
        cj = cell_types[j]
        cell_pair[ci, cj] = cell_pair[cj, ci] = cell_pair[ci, cj] + 1
        #cell_pair[cj, ci] = cell_pair[cj, ci] + 1
    }
    
    groupSize <- as.numeric(table(cellchat@idents))
    cell_pair2 = cell_pair
    diag(cell_pair2) = 0
    
    mfrow = par()$mfrow
    mar = par()$mar
    
    par(mfrow=c(3,2))
    par(mar=c(2,2,2,2))
    netVisual_circle(weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    netVisual_circle(weight*cell_pair, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Weighted number of interactions")
    
    netVisual_circle( weight*(cell_pair-cell_pair2), vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Weighted number of interactions (Within)")
    netVisual_circle( weight*cell_pair2, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Weighted number of interactions (Between)")
    
    netVisual_circle(cell_pair, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Cell pairs")
    
    par( mfrow=mfrow )
    par( mar=mar)
    
}



