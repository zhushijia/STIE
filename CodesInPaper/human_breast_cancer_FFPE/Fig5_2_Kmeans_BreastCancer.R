
source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")


for(i in 2:10)
{
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ match( as.character(spot_coordinates$barcode), as.character(cluster$Barcode)), ]
    
    if(i==6) {
        cc = data.frame( cluster = c(1,2,3,4,5,6), 
                         colors = c("steelblue",'darkgreen',"darkorange",'yellow','cyan','darkred') )
    } else if(i==7) {
        cc = data.frame( cluster = c(1,2,4,5,6,3,7), 
                         colors = c("steelblue",'darkgreen',"darkorange",'yellow','cyan','darkred',"black") )
    } else {
        cc = data.frame( cluster = c(1:i), 
                         colors = get_my_colors(i, mode=2) )
    }
    
    
    spot_cols = as.character(cc$colors)[ match( cluster$Cluster, cc$cluster ) ]
    
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/Kmeans")
    pdf(paste0("Breastcancer_spot_Kmeans_cluster",i,".pdf") )
    # png(paste0("Breastcancer_spot_Kmeans_cluster",i,".png"), w=2000, h=2000, res=300 )
    plot_sub_image(im=im, image_transparency=0, 
                   x_scale=args$x_scale, spot_coordinates=spot_coordinates, 
                   plot_spot=T, plot_cell=F, spot_cols=spot_cols, fill_spot=T  )
    
    dev.off()
}



