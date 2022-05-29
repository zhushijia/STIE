#' simulate_STdata
#'
#' simulate_STdata simulates the spot level ST data at a predefined spot size from the real ST data
#'
#' @param cells_coordinates 
#' @param cell_types a vector of character values with names representing the cell_id
#' @param Signature a matrix of numeric values representing the cell type transcriptomic signature with row being the gene and column being the cell type
#' @param spot_coordinates_ref 
#' @param spot_diameter_pixel_ref 
#' @param spot_size_ref 
#' @param spot_size_sim 
#' @param spot_dist_sim 
#' @param x_scale 
#' @param ncore 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
simulate_STdata = function( cells_coordinates, cell_types, Signature, 
                            spot_coordinates_ref, spot_diameter_pixel_ref=spot_radius*2, spot_size_ref=55, 
                            spot_size_sim, spot_dist_sim=spot_size_sim*2, x_scale=args$x_scale, ncore=30 ) {
    library(foreach)
    library(doMC)
    registerDoMC(ncore)
    
    unit_nm = spot_diameter_pixel_ref/spot_size_ref
    r = unit_nm*spot_size_sim/2
    l = unit_nm*spot_dist_sim
    
    uni_cells = cells_coordinates
    if(!is.null(cell_types)) {
        uni_cells$cell_types = cell_types[ match(as.character(uni_cells$cell_id), names(cell_types)) ]
        uni_cells = uni_cells[ !is.na(uni_cells$cell_types), ]
    }
    uni_cells = do.call(rbind, lapply( split( uni_cells, uni_cells$cell_id ), function(x) x[1,] ))
    
    intersection_area_two_circle = function(R, r, d)
    {
        r^2*acos( (d^2+r^2-R^2)/(2*d*r) ) + 
            R^2*acos( (d^2+R^2-r^2)/(2*d*R) ) - 
            sqrt( (-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R) )/2
    }
    
    
    min_pixel_x = min(spot_coordinates_ref$pixel_x) - spot_diameter_pixel_ref/2
    max_pixel_x = max(spot_coordinates_ref$pixel_x) + spot_diameter_pixel_ref/2
    min_pixel_y = min(spot_coordinates_ref$pixel_y) - spot_diameter_pixel_ref/2
    max_pixel_y = max(spot_coordinates_ref$pixel_y) + spot_diameter_pixel_ref/2
    
    M = ceiling((max_pixel_x-min_pixel_x)/l)
    N = ceiling((max_pixel_y-min_pixel_y)/l)
    
    
    sim = foreach(i = 1:M) %dopar%
    {
        cat(i, M, "\n")
        
        spot_coordinates_sim_i = NULL
        cells_on_spot_sim_i = NULL
        
        for(j in 1:N)
        {
            barcode = paste0(i,"x",j)
            
            x = min_pixel_x + (i-1)*l
            y = min_pixel_y + (j-1)*l
            
            d = sqrt( (uni_cells$pixel_x-x)^2 + (uni_cells$pixel_y-y)^2 )
            ind = which( d < (uni_cells$Major*x_scale+r) )
            n = length(ind)
            
            d2 = sqrt( (spot_coordinates_ref$pixel_x-x)^2 + (spot_coordinates_ref$pixel_y-y)^2 )
            if ( any(d2<(1.5*spot_diameter_pixel_ref)) | n>0 )
            {
                spot = data.frame(barcode=barcode, row=i, col=j, pixel_x=x, pixel_y=y, cell_count=n)
                spot_coordinates_sim_i = rbind(spot_coordinates_sim_i, spot)
            }
            
            if(n>0)
            {
                
                #spot = data.frame(barcode=barcode, row=i, col=j, pixel_x=x, pixel_y=y, cell_count=n)
                #spot_coordinates_sim_i = rbind(spot_coordinates_sim_i, spot)
                
                for( k in ind ) 
                {
                    R = uni_cells$Major[k]*x_scale
                    
                    if( abs(R-r)<d[k] )
                    {
                        area = intersection_area_two_circle(R, r, d[k])
                        percent = area / (pi*R^2) 
                    } else {
                        if(R>r) {
                            percent = (pi*r^2)  / (pi*R^2) 
                        } else {
                            percent = 1
                        }
                    }
                    
                    cell = uni_cells[k, ]
                    cell$spot = barcode
                    cell$dist = d[k]
                    cell$percent = percent 
                    
                    cells_on_spot_sim_i = rbind(cells_on_spot_sim_i, cell)
                }
                
            }
            
            cat(i, j, n, nrow(spot_coordinates_sim_i), "\n")
            
        }
        
        list( spot_coordinates_sim_i = spot_coordinates_sim_i, 
              cells_on_spot_sim_i = cells_on_spot_sim_i )
    }
    
    spot_coordinates_sim = do.call(rbind, lapply(sim, function(x) x$spot_coordinates_sim_i ))
    cells_on_spot_sim = do.call(rbind, lapply(sim, function(x) x$cells_on_spot_sim_i ))
    
    if( !is.null(Signature) & !is.null(cell_types) ) {
        
        split_cells_on_spot_sim = split(cells_on_spot_sim, cells_on_spot_sim$spot)
        
        ST_expr_sim = t( do.call(cbind, lapply( 1:length(split_cells_on_spot_sim), function(i) {
            cat(i, "\n")
            percent = split_cells_on_spot_sim[[i]]$percent
            cell = split_cells_on_spot_sim[[i]]$cell_types
            as.matrix(Signature[,cell]) %*% as.matrix(percent)
        })))
        rownames(ST_expr_sim) = names(split_cells_on_spot_sim)
        
    } else {
        
        ST_expr_sim = NULL
    }
    
    spot_coordinates_sim$imagerow = spot_coordinates_sim$pixel_y/x_scale
    spot_coordinates_sim$imagecol = spot_coordinates_sim$pixel_x/x_scale
    
    STdata_sim = list(spot_coordinates=spot_coordinates_sim, 
                      cells_on_spot=cells_on_spot_sim, 
                      ST_expr=ST_expr_sim)
    
    STdata_sim
    
}

