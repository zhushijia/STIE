The capture area has ~5,000 gene expression spots, each spot is ~55 microns with primers that include: Illumina TruSeq Read 1 (partial read 1 sequencing primer) 16 nt Spatial Barcode (all primers in a specific spot share the same Spatial Barcode) 12 nt unique molecular identifier (UMI)

simulate_STdata = function( cells_coordinates, cell_types, Signature,
                            spot_coordinates, spot_diameter_pixel=spot_radius*2, spot_size=55, 
                            spot_size_sim, spot_dist_sim=spot_size_sim*2 )
{
    unit_nm = spot_diameter_pixel/spot_size
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
    
    
    min_pixel_x = min(spot_coordinates$pixel_x) - spot_diameter_pixel/2
    max_pixel_x = max(spot_coordinates$pixel_x) + spot_diameter_pixel/2
    min_pixel_y = min(spot_coordinates$pixel_y) - spot_diameter_pixel/2
    max_pixel_y = max(spot_coordinates$pixel_y) + spot_diameter_pixel/2
    
    M = ceiling((max_pixel_x-min_pixel_x)/l)
    N = ceiling((max_pixel_y-min_pixel_y)/l)
    
    spot_coordinates_sim = NULL
    cells_on_spot_sim = NULL
    
    for(i in 1:M)
    {
        for(j in 1:N)
        {
            barcode = paste0(i,"x",j)
            
            x = min_pixel_x + (i-1)*l
            y = min_pixel_y + (j-1)*l
            
            d = sqrt( (uni_cells$pixel_x-x)^2 + (uni_cells$pixel_y-y)^2 )
            ind = which( d < (uni_cells$Major*args$x_scale+r) )
            n = length(ind)
            
            if(n>0)
            {
                spot = data.frame(barcode=barcode, row=i, col=j, pixel_x=x, pixel_y=y, cell_count=n)
                spot_coordinates_sim = rbind(spot_coordinates_sim, spot)
                
                for( k in ind ) 
                {
                    R = uni_cells$Major[k]*args$x_scale
                    
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
                    
                    cells_on_spot_sim = rbind(cells_on_spot_sim, cell)
                }
                
            }
            
            cat(i, j, n, nrow(spot_coordinates_sim), "\n")
            
        }
    }
    
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
    
    STdata_sim = list(spot_coordinates_sim=spot_coordinates, 
         cells_on_spot_sim=cells_on_spot, 
         ST_expr_sim=ST_expr)

    STdata_sim
    
}

STdata_sim = simulate_STdata(cells_coordinates=morphology_fts, cell_types=NULL, Signature=NULL,
                            spot_coordinates=spot_coordinates, spot_diameter_pixel=spot_radius*2, spot_size=55, 
                            spot_size_sim=10, spot_dist_sim=spot_size_sim*2 )

cells_coordinates=morphology_fts
cell_types=NULL
Signature=NULL
spot_coordinates=spot_coordinates
spot_diameter_pixel=spot_radius*2
spot_size=55
spot_size_sim=10
spot_dist_sim=spot_size_sim*2

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
save(STdata_sim, file=paste0("STdata_sim_radius",ri,"nm_STIEcluster",i,".RData") )




