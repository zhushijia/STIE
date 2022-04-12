deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

################################################################################################
########## load data
################################################################################################
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/Wu_etal_2021_BRCA_scRNASeq")
y = load("Wu_etal_2021_BRCA_scRNASeq_DWLS_Signature.RData")
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/Visium_FFPE_Human_Breast_Cancer/DWLS_on_Wu_etal_2021_BRCA_Signature")
x = load("Visium_FFPE_Human_Breast_Cancer_DWLS_on_Wu_etal_2021_BRCA_Signature.RData" )

prop_mat = prop_mat[, !colnames(prop_mat)%in%c("Myeloid","PVL")]
Signature = Signature[, match(colnames(prop_mat),colnames(Signature))]

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
z = load("BreastCancer_spot_BIC_new_get_cells_on_spot.RData")

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
load("BreastCancer_spot_new_get_cells_on_spot_full_signature.RData")


###########################################################
######### STIE result
############################################################

result = results[['2.8']]

stie_cell_types = result$cell_types
stie_cell_types = stie_cell_types[ match(unique(names(stie_cell_types)), names(stie_cell_types)) ]
#stie_cell_types = stie_cell_types[!stie_cell_types%in%c("Myeloid")]

stie_cell_feature = morphology_fts[ match( names(stie_cell_types), rownames(morphology_fts) ), ]
stie_cell_contour = cell_info$cell_contour[ match( names(stie_cell_types), names(cell_info$cell_contour) ) ]

cells_on_spot <- result$cells_on_spot
Signature = result$Signature

STdata = load_STdata( "Visium_FFPE_Human_Breast_Cancer", args$spaceranger_count_dir, umi_cutoff=0, is_normalized=TRUE ) 
normalized_count = STdata$normalized_matrix[[1]]
count = STdata$matrix[[1]]
Signature = Signature[rownames(Signature) %in% colnames(normalized_count),]
#Signature = Signature[, !colnames(Signature)%in%c("Myeloid")]

############################################################
########## deep learning features
############################################################
library(data.table)
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/")
cells = data.frame(fread("Breast_Cancer_FFPE_information_model40_image20.csv"))
cells = subset(cells, probability>0.5)

dl_cell_feature = with(cells, data.frame(
                    celltypes=nucleus_class,
                    probability=probability,
                    X=centroid_x0,
                    Y=centroid_y0,
                    Area=area,
                    Major=major_axis_length,
                    Minor=minor_axis_length,
                    Eccentricity=eccentricity,
                    Orientation=orientation,
                    Solidity=solidity) )
dl_cell_feature = subset(dl_cell_feature, celltypes!="blood" )
dl_cell_feature = dl_cell_feature[order(dl_cell_feature$probability,decreasing=T),]
#dl_cell_feature = subset( dl_cell_feature, celltypes!="macrophage")
#> dim(dl_cell_feature)
#[1] 65398    10

################################################################################################
########## Visualization
################################################################################################
if(0)
{
    im_scaled <- read_image(image=args$image, x_scale=args$x_scale)
    
    cells$pixel_x = cells$centroid_x0 * args$x_scale
    cells$pixel_y = cells$centroid_y0 * args$x_scale
    uni_celltypes = unique(cells$nucleus_class)
    # c("stroma", "necrosis", "lymphocyte", "blood", "tumor", "ductal epithelium", "macrophage")
    cols = c("purple", "yellow", "green", "red", "blue", "black","cyan")
    
    setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
    pdf("Breast_Cancer_FFPE_DL_celltypes_information_model40_image20.pdf")
    display(im_scaled, method='raster')
    for(i in 1:length(uni_celltypes))
    {
        X = subset( cells, nucleus_class==uni_celltypes[i] )[,c('pixel_x','pixel_y')]
        points(X,cex=0.1,col=cols[i])
    }
    plot(NA,xlim=c(0,10),ylim=c(0,10))
    legend('topright',legend=uni_celltypes,col=cols,pch=1)
    
    dev.off()
    
}

###########################################################
######### STIE and DL overlap
############################################################
overlap_index = lapply( stie_cell_contour, function(x) {
    i = with( dl_cell_feature, which( X>min(x[,1]) & X<max(x[,1]) & Y>min(x[,2]) & Y<max(x[,2]) ) )
    cat(i,"\n")
    i
} )

dl_cell_types = sapply(overlap_index, function(i) {
    cc = ""
    if(length(i)>0) {
        j = i[which.max(dl_cell_feature$probability[i])]
        cc = as.character(dl_cell_feature$celltypes[j])
    }
    cc
}) 

table2df <- function(t) {
    data.frame(ID=rownames(t), do.call(rbind, lapply(1:nrow(t),function(i) t[i,])) )
}

t1 = table(dl_cell_types, stie_cell_types)
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
#write.table(table2df(t1),"DL_STIE_concordance.txt",sep="\t",col.names=T,row.names=F,quote=F)

###########################################################
######### assign larger category
############################################################

dl_cell_types2 = dl_cell_types
dl_cell_types2[dl_cell_types%in%c("tumor", "necrosis", "ductal epithelium")] = "Epithelia"
dl_cell_types2[dl_cell_types%in%c("lymphocyte", "macrophage")] = "Immune"
dl_cell_types2[dl_cell_types%in%c("stroma")] = "Stroma"

stie_cell_types2 = stie_cell_types
stie_cell_types2[stie_cell_types%in%c("CancerEpithelial", "NormalEpithelial")] = "Epithelia"
stie_cell_types2[stie_cell_types%in%c("Tcells", "Bcells", "Plasmablasts","Myeloid")] = "Immune"
stie_cell_types2[stie_cell_types%in%c("CAFs", "Endothelial","PVL")] = "Stroma"

t2 = table(dl_cell_types2, stie_cell_types2) 
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
#write.table(table2df(t2),"DL_STIE_concordance_larger_category.txt",sep="\t",col.names=T,row.names=F,quote=F)

t2
sum(diag(t2[-1,]))/sum(t2[-1,])
sum(diag(t2[-1,]))/sum(t2)

#> sum(diag(t2[-1,]))/sum(t2[-1,])
# [1] 0.7248732
#> sum(diag(t2[-1,]))/sum(t2)
# [1] 0.6384792

############################################################
###########################################################
############################################################

stie_cell_types3 = paste0("DeepLearning_",dl_cell_types2, "_vs_STIE_",stie_cell_types2)
names(stie_cell_types3) = names(stie_cell_types2)
#stie_cell_types3 = stie_cell_types3[ stie_cell_types2!=dl_cell_types2 ]

cells_on_spot3 = cells_on_spot[ as.character(cells_on_spot$cell_id) %in% names(stie_cell_types3), ]
cells_on_spot_celltypes3 = stie_cell_types3[ match( as.character(cells_on_spot3$cell_id), names(stie_cell_types3) ) ]

t3 = table(cells_on_spot3$spot, cells_on_spot_celltypes3)
t3 = t3[rownames(t3)%in%rownames(normalized_count), ]
ST_expr = normalized_count[match( rownames(t3), rownames(normalized_count) ),]
ST_expr = ST_expr[ , match( rownames(Signature), colnames(ST_expr) ) ]


cc = do.call(cbind, lapply(1:ncol(t3),function(i) {
    top10 = order(t3[,i],decreasing=T)[1:10]
    merged_expr = t(ST_expr[top10,]) %*% apply(as.matrix.noquote(t3[top10,]),2,as.numeric)
    ci = cor(Signature,as.numeric(merged_expr[,i]),method='spearman')
    cat(colnames(merged_expr)[i],"\n")
    ci
} ) )
cc = data.frame( overlap=colnames(t3), t(cc) )

data.frame( a=colnames(t3), b=apply(cc[,-1],1,function(x)colnames(cc)[which.max(x)+1]) )

setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
#write.table(cc,"DL_STIE_scRNA_Signature_concordance_larger_category.txt",sep="\t",col.names=T,row.names=F,quote=F)


top10_cells = sapply(1:ncol(t3),function(i) {
    top10 = order(t3[,i],decreasing=T)[1:10]
    sum(t3[top10,i])
} )

all_cells = sapply( strsplit( colnames(t3), "_vs_" ), function(x) {
    i = which( rownames(t2) == gsub("DeepLearning_","",x[1]) )
    j = which( colnames(t2) == gsub("STIE_","",x[2]) )
    t2[i,j]
} )

c2 = data.frame(cc,top10_cells,all_cells)
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
#write.table(c2,"DL_STIE_scRNA_Signature_concordance_larger_category_allinfo.txt",sep="\t",col.names=T,row.names=F,quote=F)



c3 = c2[,-1]
rownames(c3) = as.character(c2[,1])
data.frame( apply(c3,1,function(x) colnames(c3)[which.max(x[1:8])] ), 
            c3$all_cells)



