
source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_hippocampus/parameters_hippo3.R")

for(i in 2:10)
{
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ match( as.character(spot_coordinates$barcode), as.character(cluster$Barcode)), ]
    
    if(i==5) {
        cc = data.frame( cluster = c(1:5), 
                         colors = c("orange", "green", "black", "magenta", "blue") )
    } else {
        cc = data.frame( cluster = c(1:i), 
                         colors = get_my_colors(i, mode=2) )
    }
    
    spot_cols = as.character(cc$colors)[ match( cluster$Cluster, cc$cluster ) ]
    
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_hippocampus3/results/Kmeans")
    pdf(paste0("hippo_xoff8000_yoff10500_spot_Kmeans_cluster",i,".pdf") )
    # png(paste0("hippo_xoff8000_yoff10500_spot_Kmeans_cluster",i,".png"), w=2000, h=2000, res=300 )
    plot_sub_image(im=im, w=6000, h=5000, xoff=8000, yoff=10500, 
                   x_scale=0.1, spot_coordinates=spot_coordinates, 
                   plot_spot=T, plot_cell=F, spot_cols=spot_cols, fill_spot=T  )
    
    dev.off()
}



