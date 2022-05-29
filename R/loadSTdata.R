#' load_STdata
#'
#' @param sample_names a character value representing the sample name
#' @param sample_dir a character value representing the path where spaceranger generate the read count
#' @param umi_cutoff an integer value representing the UMI count to filter spots
#' @param is_normalized a boolean value indicating whether to normalize the raw read count
#'
#' @return A list containing the follow components:
#' \itemize{
#'  \item {matrix} a matrix of raw read count on the spot
#'  \item {normalized_matrix} a matrix of normalized read count on the spot
#'  \item {bcs_merge} a data frame representing the spot information
#'  \item {images_tibble} an image tibble object for the ST matched image
#' }
#' 
#' @author Shijia Zhu, \email{shijia.zhu@@utsouthwestern.edu}
#' @export
#'
#' @references
#'
#' 
#' 
load_STdata <- function( sample_names, sample_dir, umi_cutoff=0, is_normalized=TRUE ) {
    
    image_paths <- paste( sample_dir, "outs/spatial/tissue_lowres_image.png", sep="/")
    scalefactor_paths <- paste( sample_dir, "outs/spatial/scalefactors_json.json", sep="/")
    tissue_paths <- paste( sample_dir, "outs/spatial/tissue_positions_list.csv", sep="/")
    cluster_paths <- paste( sample_dir, "outs/analysis/clustering/graphclust/clusters.csv", sep="/")
    matrix_paths <- paste( sample_dir, "outs/filtered_feature_bc_matrix.h5", sep="/")
    
    images_cl <- list()
    for (i in 1:length(sample_names)) {
        images_cl[[i]] <- read.bitmap(image_paths[i])
    }
    
    height <- list()
    for (i in 1:length(sample_names)) {
        height[[i]] <-  data.frame(height = nrow(images_cl[[i]]))
    }
    height <- bind_rows(height)
    
    width <- list()
    for (i in 1:length(sample_names)) {
        width[[i]] <- data.frame(width = ncol(images_cl[[i]]))
    }
    width <- bind_rows(width)
    
    grobs <- list()
    for (i in 1:length(sample_names)) {
        grobs[[i]] <- rasterGrob(images_cl[[i]], width=unit(1,"npc"), height=unit(1,"npc"))
    }
    
    images_tibble <- tibble(sample=factor(sample_names), grob=grobs)
    images_tibble$height <- height$height
    images_tibble$width <- width$width
    
    
    scales <- list()
    for (i in 1:length(sample_names)) {
        scales[[i]] <- rjson::fromJSON(file = scalefactor_paths[i])
    }
    
    clusters <- list()
    for (i in 1:length(sample_names)) {
        clusters[[i]] <- read.csv(cluster_paths[i])
    }
    
    bcs <- list()
    
    for (i in 1:length(sample_names)) {
        bcs[[i]] <- read.csv(tissue_paths[i],col.names=c("barcode","tissue","row","col","imagerow","imagecol"), header = FALSE)
        bcs[[i]]$imagerow <- bcs[[i]]$imagerow * scales[[i]]$tissue_lowres_scalef    # scale tissue coordinates for lowres image
        bcs[[i]]$imagecol <- bcs[[i]]$imagecol * scales[[i]]$tissue_lowres_scalef
        bcs[[i]]$tissue <- as.factor(bcs[[i]]$tissue)
        bcs[[i]] <- merge(bcs[[i]], clusters[[i]], by.x = "barcode", by.y = "Barcode", all = TRUE)
        bcs[[i]]$height <- height$height[i]
        bcs[[i]]$width <- width$width[i]
    }
    
    names(bcs) <- sample_names
    
    matrix <- list()
    for (i in 1:length(sample_names)) {
        matrix[[i]] <- as.data.frame(t(as.matrix(Read10X_h5(matrix_paths[i]))))
    }
    
    umi_sum <- list() 
    for (i in 1:length(sample_names)) {
        umi_sum[[i]] <- data.frame(barcode =  row.names(matrix[[i]]),
                                   sum_umi = Matrix::rowSums(matrix[[i]]))
        
    }
    names(umi_sum) <- sample_names
    umi_sum <- bind_rows(umi_sum, .id = "sample")
    
    gene_sum <- list() 
    for (i in 1:length(sample_names)) {
        gene_sum[[i]] <- data.frame(barcode =  row.names(matrix[[i]]),
                                    sum_gene = Matrix::rowSums(matrix[[i]] != 0))
    }
    names(gene_sum) <- sample_names
    gene_sum <- bind_rows(gene_sum, .id = "sample")
    
    bcs_merge <- bind_rows(bcs, .id = "sample")
    bcs_merge <- merge(bcs_merge,umi_sum, by = c("barcode", "sample"))
    bcs_merge <- merge(bcs_merge,gene_sum, by = c("barcode", "sample"))
    
    names(matrix) = sample_names
    
    
    if(umi_cutoff>0) {
        
        bcs_merge = subset(bcs_merge, sum_umi>=umi_cutoff)
        
        for(i in 1:length(matrix)) {
            bcs_selected = subset( bcs_merge, sample==names(matrix)[i] & barcode%in%rownames(matrix[[i]]) )
            matrix[[i]] = matrix[[i]][ rownames(matrix[[i]])%in%as.character(bcs_selected$barcode), ]
        }
    }
    
    if( is_normalized ) {
        normalized_matrix <- list()
        for (i in 1:length(sample_names)) {
            cat('normalizing.. ', names(matrix)[i] ,'\n')
            counts <- matrix[[i]][ rowSums(matrix[[i]])>0, ]
            seurat_obj <- CreateSeuratObject(counts = t(counts) )
            seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
            normalized_data <- as.matrix(seurat_obj$SCT@data)
            normalized_matrix[[i]] <- data.frame(t(normalized_data))
        }
        names(normalized_matrix) = names(matrix)
    } else {
        normalized_matrix = NULL
    }
    
    ST_data = list(matrix=matrix, normalized_matrix=normalized_matrix, bcs_merge=bcs_merge, images_tibble=images_tibble)
    
    ST_data
}










