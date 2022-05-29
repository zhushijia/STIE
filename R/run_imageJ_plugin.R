#' run_imageJ_plugin
#'
#' run_imageJ_plugin runs the imageJ plugin in R
#'
#' @param imageJ a character value representing the path to imageJ excutable file
#' @param plugin_macro a character value representing the path to the imageJ macro file (.ijm)
#' @param split_image_dir a character value representing the path to save the split images
#' @param feature_dir a character value representing the path to save the cell segmentation and extracted morphological features
#' @param pattern a regular expression to select the image files of interest
#'
#' @return
#' @export
#' 
#' @author Shijia Zhu, \email{shijia.zhu@@utsouthwestern.edu}
#'
#' @seealso \code{\link{split_image}}; \code{\link{merge_feature}}; 
#' 
#' 
run_imageJ_plugin <- function( imageJ, plugin_macro=NULL, split_image_dir, feature_dir=split_image_dir, pattern="jpg$" )
{
    
    #imageJ = "/archive/SCCC/Hoshida_lab/s184554/Code/github/ImageJ/Fiji.app/ImageJ-linux64" 
    #plugin_macro = "/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/Code/STIE/v1.1/data/DeepImageJ_plugin_multi_organ_3000_3000.fiji.ijm"
    
    #split_image_dir="/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/MouseLiverNatComms/sample_1/img/high_resolution/CN73_Liver_HE_C1_0.1_split_images"

    if( is.null(plugin_macro) )
    {
        STIE.dir = system.file(package = "STIE")
        plugin_macro = paste0(STIE.dir , "/data/DeepImageJ_plugin_multi_organ_3000_3000.fiji.ijm")
    }
    
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







