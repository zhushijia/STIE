deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source("/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/human_breast_cancer_FFPE/parameter_BreastCancer.R")

################################################################################################
########## load data
################################################################################################
STdata = load_STdata( "Visium_FFPE_Human_Breast_Cancer", args$spaceranger_count_dir, umi_cutoff=0, is_normalized=TRUE ) 
normalized_count = STdata$normalized_matrix[[1]]
count = STdata$matrix[[1]]

thres = 0.5
setwd("/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/HumanBreastCancer_FFPE/count/results/STIE")
x = load( paste0("BreastCancer_spot_DL_CellSegmentation_thres",thres,".RData") )

0.5 + 2.5: ####
# spearman: DeepLearning_Epithelia_vs_STIE_Stroma             Myeloid  1400
# pearson: #
#    DeepLearning_Epithelia_vs_STIE_Immune                CAFs   700
#    DeepLearning_Stroma_vs_STIE_Immune                   CAFs  2500

0.5 + 3: 
# spearman: 
#    DeepLearning_Epithelia_vs_STIE_Immune    CancerEpithelial   700
#    DeepLearning_Epithelia_vs_STIE_Stroma             Myeloid  1400 ###
#    DeepLearning_Stroma_vs_STIE_Immune                   CAFs  2500
# pearson: 
#    DeepLearning_Epithelia_vs_STIE_Immune                CAFs   700
#    DeepLearning_Stroma_vs_STIE_Immune                   CAFs  2500

0 + 2.5: 
# spearman: #
DeepLearning_Epithelia_vs_STIE_Immune    CancerEpithelial  1900
DeepLearning_Epithelia_vs_STIE_Stroma             Myeloid  4200
# pearson: 
DeepLearning_Epithelia_vs_STIE_Immune                CAFs  1900
DeepLearning_Stroma_vs_STIE_Immune                   CAFs 12000

0 + 3
# spearman: 
DeepLearning_Stroma_vs_STIE_Immune                   CAFs 12000
DeepLearning_Epithelia_vs_STIE_Immune    CancerEpithelial  1900
# pearson: 
DeepLearning_Stroma_vs_STIE_Immune                   CAFs 12000
DeepLearning_Epithelia_vs_STIE_Immune    CancerEpithelial  1900

###########################################################
######### STIE result
############################################################
spot = '2.5'
result = results[[spot]]
corr_method = 'spearman'

dl_cell_types = as.character(result$cells_on_spot$cell_types)
stie_cell_types = result$cell_types
names(dl_cell_types) = names(stie_cell_types) = as.character(result$cells_on_spot$cell_id)

uni_cell_ids = unique(as.character(result$cells_on_spot$cell_id))
ind = match( uni_cell_ids , names(dl_cell_types)  )
dl_cell_types = dl_cell_types[ind]
stie_cell_types = stie_cell_types[ind]


cells_on_spot <- result$cells_on_spot
Signature = result$Signature
Signature = Signature[rownames(Signature) %in% colnames(normalized_count),]
#Signature = Signature[, !colnames(Signature)%in%c("Myeloid")]

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
sum(diag(t2))/sum(t2)

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
    ci = cor(Signature,as.numeric(merged_expr[,i]),method=corr_method)
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



c3 = signif(c2[,-1],2)
rownames(c3) = as.character(c2[,1])
n = ncol(c3)-2
data.frame( x=apply(c3,1,function(x) colnames(c3)[which.max(x[1:n])] ), 
            y=c3$all_cells)



