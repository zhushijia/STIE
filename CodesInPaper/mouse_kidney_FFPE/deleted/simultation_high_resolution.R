deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
x = load("HumanBreastCancer_clustering.RData")


# https://www.calculator.net/triangle-calculator.html?vc=60&vx=10&vy=20&va=&vz=18&vb=&angleunits=d&x=124&y=18
# https://mathworld.wolfram.com/Circle-CircleIntersection.html


STIE_result = result[[6]]
ri=5

simulate_STdata = function( STIE_result, spot_coordinates, ri=5 )
{
    unit_nm = spot_radius*2/55
    #ri = 55 # diameter
    r = unit_nm*ri/2
    l = unit_nm*ri*2
    
    intersection_area_two_circle = function(R, r, d)
    {
        r^2*acos( (d^2+r^2-R^2)/(2*d*r) ) + 
            R^2*acos( (d^2+R^2-r^2)/(2*d*R) ) - 
            sqrt( (-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R) )/2
    }
    
    uni_cells = STIE_result$cells_on_spot
    uni_cells$cell_types = STIE_result$cell_types
    uni_cells = do.call(rbind, lapply( split( uni_cells, uni_cells$cell_id ), function(x) x[1,] ))
    dim(uni_cells)
    Signature = STIE_result$Signature
    
    min_pixel_x = min(spot_coordinates$pixel_x) - spot_radius
    max_pixel_x = max(spot_coordinates$pixel_x) + spot_radius
    min_pixel_y = min(spot_coordinates$pixel_y) - spot_radius
    max_pixel_y = max(spot_coordinates$pixel_y) + spot_radius
    
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
    
    split_cells_on_spot_sim = split(cells_on_spot_sim, cells_on_spot_sim$spot)
    
    ST_expr_sim = t( do.call(cbind, lapply( 1:length(split_cells_on_spot_sim), function(i) {
        cat(i, "\n")
        percent = split_cells_on_spot_sim[[i]]$percent
        cell = split_cells_on_spot_sim[[i]]$cell_types
        as.matrix(Signature[,cell]) %*% as.matrix(percent)
    })))
    rownames(ST_expr_sim) = names(split_cells_on_spot_sim)
    
    STdata_sim = list(spot_coordinates_sim=spot_coordinates_sim, 
         cells_on_spot_sim=cells_on_spot_sim, 
         ST_expr_sim=ST_expr_sim)

    STdata_sim
    
}


setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
save(STdata_sim, file=paste0("STdata_sim_radius",ri,"nm_STIEcluster",i,".RData") )

##########################################################################################
plot(NA, xlim=c(min_pixel_x, max_pixel_x), ylim=c(min_pixel_y, max_pixel_y) )
points(spot_coordinates_sim$pixel_x, spot_coordinates_sim$pixel_y)

barplot( table(spot_coordinates_sim$cell_count ))
plot( density(cells_on_spot_sim$percent ))

x=with(cells_on_spot_sim, table(spot, cell_types))
(y = table(apply(x,1,function(a)sum(a>0))))
barplot(y)

X = kmeans(ST_expr_sim, 6)
Signature_sim = t( apply(ST_expr_sim, 2, function(x) tapply(x,X$cluster,mean) ))
prop=lapply( 1:6, function(i) solveOLS( Signature, as.matrix(Signature_sim[,i]) ) )
data.frame(sapply(prop,which.max),sapply(prop,max))
lapply(prop,sort,decreasing=T)


real_cluster = with( cells_on_spot_sim, tapply(cell_types, spot, function(x) paste(sort(unique(x)),collapse="_") ))
sim_cluster = X$cluster
sim_cluster = sim_cluster[ match( names(real_cluster), names(sim_cluster) ) ]
table(real_cluster, sim_cluster)

##########################################################################################
ri = 10
i=6
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load( paste0("STdata_sim_radius",ri,"nm_STIEcluster",i,".RData") )

ST_expr_sim = STdata_sim$ST_expr_sim
cells_on_spot_sim = STdata_sim$cells_on_spot_sim

X = kmeans(ST_expr_sim, i)
Signature_sim = t( apply(ST_expr_sim, 2, function(x) tapply(x,X$cluster,mean) ))

lambda=0
steps=30
morphology_steps=10
known_signature=FALSE
known_cell_types=FALSE
min_cells=2

res = STIE(ST_expr_sim, Signature_sim, cells_on_spot_sim, features, 
           lambda=0, steps=30, 
           known_signature=FALSE, known_cell_types=FALSE)


X = kmeans(ST_expr_sim, 6)
Signature_sim = t( apply(ST_expr_sim, 2, function(x) tapply(x,X$cluster,mean) ))
prop=lapply( 1:6, function(i) solveOLS( Signature, as.matrix(Signature_sim[,i]) ) )
data.frame(sapply(prop,which.max),sapply(prop,max))
lapply(prop,sort,decreasing=T)


Signature = result[[6]]$Signature
prop=lapply( 1:6, function(i) solveOLS( Signature, as.matrix(res$Signature[,i]) ) )
data.frame(sapply(prop,which.max),sapply(prop,max))
lapply(prop,sort,decreasing=T)





