deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

dwls = function(S, B)
{
    library(Matrix)
    source("/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/DWLS/Deconvolution_functions.R")
    
    nb = ncol(B)
    ns = ncol(S)
    prop_mat <- as.data.frame(matrix(0, nrow = nb, ncol = ns) )
    rownames(prop_mat) <- colnames(B)
    colnames(prop_mat) <- colnames(S)
    for (ib in 1:nb) {
        print( sprintf("Estimating proportion for spot : %d / %d", ib, nb) )
        b <- B[,ib]
        tr <- trimData(S, b)
        tr$sig <- tr$sig[,colSums(tr$sig) > 0]
        is_pd <- eigen(t(tr$sig)%*%tr$sig)$values
        is_pd <- all(is_pd > 10e-6)
        if (!(is_pd)) { next}
        try(solDWLS <- solveDampenedWLS(tr$sig,tr$bulk),next)
        print("Proportions >> ")
        prop_mat[ib,names(solDWLS)] <- solDWLS
    }
    
    t(prop_mat)
    
}

##########################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
#x = load("HumanBreastCancer_clustering.RData")
x = load("HumanBreastCancer_clustering_2.5xSpot_lambda1000.RData")

k = 6
ri = 5
la = 0
result = results[[k]]

STdata_sim <- simulate_STdata(cells_coordinates=result$cells_on_spot, 
                              cell_types=result$cell_types, 
                              Signature=result$Signature,
                              spot_coordinates_ref=spot_coordinates, 
                              spot_diameter_pixel_ref=spot_radius*2, 
                              spot_size_ref=55, 
                              spot_size_sim=ri, 
                              x_scale=args$x_scale )

##########################################################################################
Signature_sim = result$Signature
ST_expr_sim = STdata_sim$ST_expr
cells_on_spot_sim = STdata_sim$cells_on_spot
spot_coordinates_sim = STdata_sim$spot_coordinates
spot_coordinates_sim = spot_coordinates_sim[ match( rownames(ST_expr_sim), as.character(spot_coordinates_sim$barcode) ) , ]
##########################################################################################

set.seed(123)
Km = kmeans(ST_expr_sim, k)
Signature_kmeans = t( apply(ST_expr_sim, 2, function(x) tapply(x,Km$cluster,mean) ))

stie = STIE(ST_expr_sim, Signature_kmeans, cells_on_spot_sim, features, 
           lambda=la, steps=30, 
           known_signature=FALSE, known_cell_types=FALSE)
Signature_stie = stie$Signature

##########################################################################################
kmeans_nnls = sapply( 1:ncol(Signature_kmeans), function(i) solveNNLS( Signature_sim, Signature_kmeans[,i] ) )
stie_nnls = sapply( 1:ncol(Signature_stie), function(i) solveNNLS( Signature_sim, Signature_stie[,i] ) )
kmeans_dwls = dwls(S=Signature_sim, B=Signature_kmeans)
stie_dwls = dwls(S=Signature_sim, B=Signature_stie)

index = apply(stie_nnls,2,which.max)

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
save( ST_expr_sim, cells_on_spot_sim, spot_coordinates_sim,
      Signature_sim, Signature_kmeans, Signature_stie, 
      stie, Km,
      kmeans_nnls, stie_nnls, kmeans_dwls, stie_dwls,
      file=paste0("Breastcancer_simulation_highresolution_",ri,"nm_from_STIE_clustering",k,"_reRunWtih_lambda",la,".RData") )


##########################################################################################
###### Visualization
##########################################################################################

k = 6
ri = 10
la = 0

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
x = load("HumanBreastCancer_clustering_2.5xSpot_lambda1000.RData")
y = load( paste0("Breastcancer_simulation_highresolution_",ri,"nm_from_STIE_clustering",k,"_reRunWtih_lambda",la,".RData") )
result = results[[k]]

table2df <- function(t) {
    data.frame(ID=rownames(t), do.call(rbind, lapply(1:nrow(t),function(i) t[i,])) )
}

index = apply(stie_nnls,2,which.max)
r = result$cell_types
s = stie$cell_types
n = unique(names(r))
r = r[match(n,names(r))]
s = s[match(n,names(s))]
t2 = table(r,s)[index,]
sum(diag(t2))/sum(t2)

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE/clustering_visualization")
write.table(table2df(t2), paste0("overlap_STIEvsKmeans_Breastcancer_simulation_highresolution_",ri,"nm_from_STIE_clustering",k,"_reRunWtih_lambda",la,".txt"),
            sep="\t",col.names=T,row.names=F,quote=F)

##########################################################################################
pdf( paste0("Breastcancer_simulation_highresolution_",ri,"nm_from_STIE_clustering",k,"_reRunWtih_lambda",la,".pdf"), w=10, h=8 )

myCol = c("steelblue",'darkgreen',"darkorange",'yellow','cyan','darkred')

layout(matrix(1))
par(mar=c(2,2,2,2))
par(mfrow=c(2,3))
barplot(kmeans_nnls[index, ], col=myCol, ylim=c(0,1), border='white')
barplot(stie_nnls[index, ], col=myCol, ylim=c(0,1), border='white')
plot.new()
barplot(kmeans_dwls[index, ], col=myCol, ylim=c(0,1), border='white')
barplot(stie_dwls[index, ], col=myCol, ylim=c(0,1), border='white')

# true low resolution STIE clustering
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
plot_sub_image(im=im, image_transparency=0,
               x_scale=args$x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=myCol, 
               plot_spot=F, plot_cell=T, 
               axis_tick=0, axis_col='grey'  )

# estimated high resolution STIE clustering
cell_types = stie$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
plot_sub_image(im=im, image_transparency=0,
               x_scale=args$x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=myCol[index], 
               plot_spot=F, plot_cell=T, 
               axis_tick=0, axis_col='grey'  )

# high resolution kmeans
spot_cols = myCol[ match( Km$cluster, index ) ]
plot_sub_image(im=im, 
               x_scale=args$x_scale, spot_coordinates=spot_coordinates_sim, 
               contour=contour2, plot_spot=T, plot_cell=F,
               spot_radius=spot_radius/args$x_scale/55*ri, 
               spot_cols=spot_cols, fill_spot=F  )

dev.off()








