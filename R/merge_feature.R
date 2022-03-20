#' merge_feature
#'
#' @param feature_dir 
#' @param tag 
#'
#' @return
#' @export
#'
#' @examples
merge_feature <- function( feature_dir, tag=NULL )
{
  
  feature_files <- paste0( feature_dir, "/", list.files(feature_dir,pattern=".feature$") )
  mask_files <- gsub("feature$","mask.txt",feature_files)
  
  features_noext <- paste0( feature_dir, "/", sapply( strsplit(basename(feature_files), "[.]"), function(x) x[1] ) )
  
  xoff <- sapply( strsplit(features_noext, "_"), function(x) {
      i = which(x=='xoff')
      as.integer(x[i+1])
  })

  yoff <- sapply( strsplit(features_noext, "_"), function(x) {
      i = which(x=='yoff')
      as.integer(x[i+1])
  })
  
  features <- do.call(rbind, lapply(1:length(feature_files), function(i) {
    cat("loading morphology feature for image split:", i,"\n")
    feature <- read.delim(feature_files[i],sep="\t",header=T)
    if( !"Eccentricity"%in%colnames(feature) ) feature$Eccentricity <- with(feature, sqrt( 1-Minor^2/Major^2 ) )
    feature$X = feature$X + xoff[i]
    feature$Y = feature$Y + yoff[i]
    feature
  }))

  contours <- do.call(c, lapply(1:length(feature_files), function(i) {
      cat("extracting cell contour for image split:", i,"\n")
      mask_mat <- as.matrix(fread(mask_files[i]))
      contour <- ocontour(Image(mask_mat))
      
      #### contour: x = 2nd column, y=1st col
      contour_adjusted = lapply( contour, function(x) data.frame(x=x[,2]+xoff[i],y=x[,1]+yoff[i]) )
      contour_adjusted
  }))
  
  cat("checking duplicates ... ...\n")
  duplicates <- do.call(c, lapply(1:nrow(features), function(i) {
    disti <- sqrt( (features$X-features$X[i])^2 + (features$Y-features$Y[i])^2 )
    disti[i] <- Inf
    j = which( disti<(features$Width/2) | disti<(features$Height/2) )
    du = NULL
    if(length(j)>0) {
      ind = c(i,j)
      #du = ind[ which.min(features$Area[ind]) ]
      du = ind[ -which.max(features$Area[ind]) ]
    } 
    du
  } ))
  
  cat("removing duplicates ... ...\n")
  features = features[-duplicates, ]  
  contours = contours[-duplicates]
  
  cell_id = 1:nrow(features)
  features =  data.frame(cell_id,features)
  names(contours) = cell_id
  
  info <- list(cell_feature=features, cell_contour=contours)
  
  info
}



# https://imagej.net/plugins/morphological-segmentation

# labelsToROIs
#https://www.youtube.com/watch?v=G_vGQCzZ3xg

#
#https://www.unige.ch/medecine/bioimaging/files/3714/1208/5964/CellCounting.pdf

#https://forum.image.sc/t/creating-labeled-image-from-rois-in-the-roi-manager/4256/3

