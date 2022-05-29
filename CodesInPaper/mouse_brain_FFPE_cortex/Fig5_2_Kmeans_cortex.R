
source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/mouse_brain_FFPE_cortex/parameters_cortex.R")

for(i in 2:10)
{
    cluster = read.csv(paste0(args$spaceranger_count_dir, "/outs/analysis/clustering/kmeans_",i,"_clusters/clusters.csv") )
    cluster = cluster[ match( as.character(spot_coordinates$barcode), as.character(cluster$Barcode)), ]
    
    if(i==6) {
        cc = data.frame( cluster = c(1:6), 
                         colors = c("red", "steelblue", "blue", "black", "yellow", "cyan") )
    } else {
        cc = data.frame( cluster = c(1:i), 
                         colors = get_my_colors(i, mode=2) )
    }
    
    spot_cols = as.character(cc$colors)[ match( cluster$Cluster, cc$cluster ) ]
    
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/AdultMouseBrain_FFPE/count_cortex/results/Kmeans")
    pdf(paste0("hippo_xoff14000_yoff7000_spot_Kmeans_cluster",i,".pdf") )
    plot_sub_image(im=im, w=14000, h=7000, xoff=6000, yoff=8000, 
                   x_scale=0.05, spot_coordinates=spot_coordinates, 
                   plot_spot=T, plot_cell=F, spot_cols=spot_cols, fill_spot=T  )
    
    dev.off()
}

