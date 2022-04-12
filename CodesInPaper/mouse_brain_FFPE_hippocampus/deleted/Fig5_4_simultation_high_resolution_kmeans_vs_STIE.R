deconvolution = TRUE
source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo3.R")

############################################################
## deconvolution on multiorgan
############################################################

cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)

STIE_result <- STIE(ST_expr, Signature, cells_on_spot, features, 
               lambda=0, steps=30, 
               known_signature=TRUE, known_cell_types=FALSE)

if(1) {
    
    cells_on_spot$cell_types = result$cell_types
    spot_id = as.character(cells_on_spot$spot)
    cell_id = as.character(cells_on_spot$cell_id)
    ST_expr2 = t( ST_expr[ rownames(ST_expr)%in%spot_id,  ] )
    
    coefs = table( as.character(cells_on_spot$spot), as.character(cells_on_spot$cell_types) )
    coefs = coefs[match( colnames(ST_expr2), rownames(coefs) ), ]
    Signature_sim = t( apply( ST_expr2, 1, function(x) solveNNLS( coefs, as.matrix(x), scaled=F ) ) )
    STIE_result$Signature = Signature_sim
}

STdata_sim <- simulate_STdata(cells_coordinates=STIE_result$cells_on_spot, 
                              cell_types=STIE_result$cell_types, 
                              Signature=STIE_result$Signature,
                              spot_coordinates_ref=spot_coordinates, 
                              spot_diameter_pixel_ref=spot_radius*2, 
                              spot_size_ref=55, 
                              spot_size_sim=5, 
                              x_scale=args$x_scale )



setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")
load("MouseBrainHippocampus_clustering_2.5Xspot.RData")



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

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/STIE")
save(STdata_sim, file=paste0("STdata_sim_radius",ri,"nm_STIEcluster",i,".RData") )

##########################################################################################
plot(NA, xlim=c(min_pixel_x, max_pixel_x), ylim=c(min_pixel_y, max_pixel_y) )
points(spot_coordinates_sim$pixel_x, spot_coordinates_sim$pixel_y)

barplot( table(spot_coordinates_sim$cell_count ))
plot( density(cells_on_spot_sim$percent ))

x=with(cells_on_spot_sim, table(spot, cell_types))
(y = table(apply(x,1,function(a)sum(a>0))))
barplot(y)

X = kmeans(ST_expr_sim, k)
Signature_sim = t( apply(ST_expr_sim, 2, function(x) tapply(x,X$cluster,mean) ))
prop=lapply( 1:k, function(i) solveNNLS( Signature, as.matrix(Signature_sim[,i]) ) )
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
prop=lapply( 1:6, function(i) solveNNLS( Signature, as.matrix(Signature_sim[,i]) ) )
data.frame(sapply(prop,which.max),sapply(prop,max))
lapply(prop,sort,decreasing=T)


Signature = result[[6]]$Signature
prop=lapply( 1:6, function(i) solveNNLS( Signature, as.matrix(res$Signature[,i]) ) )
data.frame(sapply(prop,which.max),sapply(prop,max))
lapply(prop,sort,decreasing=T)











##########################################################################################




k = 5
STIE_result = results[[k]]
Signature_sim = STIE_result$Signature
ST_expr_sim = STdata_sim$ST_expr_sim
cells_on_spot_sim = STdata_sim$cells_on_spot_sim

##########################################################################################
k=6
Km = kmeans(ST_expr_sim, k)
Signature_kmeans = t( apply(ST_expr_sim, 2, function(x) tapply(x,Km$cluster,mean) ))
prop_kmeans = lapply( 1:ncol(Signature_kmeans), function(i) solveNNLS( Signature_sim, Signature_kmeans[,i] ) )
data.frame(sapply(prop_kmeans,which.max),sapply(prop_kmeans,max))
lapply(prop_kmeans,sort,decreasing=T)

stie = STIE(ST_expr_sim, Signature_kmeans, cells_on_spot_sim, features, 
           lambda=1e3, steps=30, 
           known_signature=FALSE, known_cell_types=FALSE)
Signature_stie = stie$Signature
##########################################################################################

prop_kmeans = lapply( 1:ncol(Signature_kmeans), function(i) solveNNLS( Signature_sim, Signature_kmeans[,i] ) )
data.frame(sapply(prop_kmeans,which.max),sapply(prop_kmeans,max))
lapply(prop_kmeans,sort,decreasing=T)

prop_stie = lapply( 1:ncol(Signature_stie), function(i) solveNNLS( Signature_sim, Signature_stie[,i] ) )
data.frame(sapply(prop_stie,which.max),sapply(prop_stie,max))
lapply(prop_stie,sort,decreasing=T)


##########################################################################################
real_cluster = with( cells_on_spot_sim, tapply(cell_types, spot, function(x) paste(sort(unique(x)),collapse="_") ))
sim_cluster = Km$cluster
sim_cluster = sim_cluster[ match( names(real_cluster), names(sim_cluster) ) ]
table(real_cluster, sim_cluster)












