
run_imageJ_plugin <- function( imageJ, plugin_macro, split_image_dir, feature_dir=split_image_dir, pattern="jpg$" )
{
    
    #imageJ = "/archive/SCCC/Hoshida_lab/s184554/Code/github/ImageJ/Fiji.app/ImageJ-linux64" 
    #plugin_macro = "/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/Code/STIE/v1.1/data/DeepImageJ_plugin_multi_organ_3000_3000.fiji.ijm"
    
    #split_image_dir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/MouseLiverNatComms/sample_1/img/high_resolution/CN73_Liver_HE_C1_0.1_split_images"

    setwd(split_image_dir)
    
    images <- paste0( split_image_dir, "/", list.files(split_image_dir,pattern=pattern) )
    
    images <- images[ !grepl('.mask.jpg!',images) ]
    done <- sapply(images, function(im) file.exists( paste0(im,".feature")) )
    images <- images[!done]
    
    setwd(feature_dir)
    for(i in 1:length(images))
    {
        cmd <- paste( imageJ, "-macro", plugin_macro, images[i] )
        system(cmd)
    }

    cell_info <- merge_feature( feature_dir )
    setwd(feature_dir)
    save(cell_info, file="cell_info.RData")
    
}







