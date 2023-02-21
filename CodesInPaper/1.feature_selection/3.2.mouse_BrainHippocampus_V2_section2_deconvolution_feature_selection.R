deconvolution = TRUE
clustering = FALSE
signature_learning = FALSE

source( "/archive/SCCC/Hoshida_lab/s184554/Code/github/STIE/CodesInPaper/parameters/parameters_10X_V2_mouse_hippo_section2.R")
cells_on_spot <- get_cells_on_spot( cell_coordinates=morphology_fts, spot_coordinates, 2.5*spot_radius)
outputDir = "/archive/SCCC/Hoshida_lab/shared/fastq/SpatialTranscriptome/10X_public_dataset/V2_MouseBrainCoronalSection2_FFPE/CytAssist_FFPE_Mouse_Brain_Rep2/results/STIE_search"

criterion = "rmse"
outputDir = paste0(outputDir, "_", criterion)
lambdas=c(0,1e1,1e2,1e3,1e4,1e5,1e6)

############################################################
## run deconvolution
############################################################
dir.create(outputDir)
setwd(outputDir)

pdf("V2_MouseBrainHippocampus_STIE_searchPaths_deconvolution.pdf")
paths = STIE_search(ST_expr, Signature, cells_on_spot, 
                    steps=30, known_signature=TRUE, known_cell_types=FALSE, 
                    lambdas=lambdas, min_cells=5,
                    criterion = criterion )
dev.off()

save(paths, file="V2_MouseBrainHippocampus_STIE_searchPaths_deconvolution.RData")

############################################################
## plot 
############################################################
setwd(outputDir)
load("V2_MouseBrainHippocampus_STIE_searchPaths_deconvolution.RData")
result = paths[[1]]$best_result[[2]]
cell_types = result$cell_types
contour2 = cell_info$cell_contour[ match(names(cell_types), names(cell_info$cell_contour)) ]
colors = c("magenta", "blue", "green", "black", "orange", "cyan")
plot_sub_image(im=im,  w=8700, h=5200, xoff=10500, yoff=7000,
               x_scale=x_scale, spot_coordinates=spot_coordinates, 
               contour=contour2, cell_types=cell_types, color_use=colors, plot_spot=T, plot_cell=T)

