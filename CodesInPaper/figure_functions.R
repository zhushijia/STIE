get_PME_diff = function(results) 
{   
    lapply( results, function(res) {
        PE_on_spot = res$PE_on_spot
        PM_on_spot = apply(res$PM_on_cell, 2, function(x) tapply(x, as.character(res$cells_on_spot$spot), mean) )
        PM_on_spot = t( apply(PM_on_spot,1,function(x)x/sum(x)) )
        PM_on_spot = PM_on_spot[ match( rownames(PE_on_spot), rownames(PM_on_spot) ) , ]
        rowSums((PE_on_spot-PM_on_spot)^2)
    } )
}


errplot <- function(X, ylab="", line=T, ylim=NULL)
{
    m <- sapply(X,function(x) mean(x) )
    se <- sapply(X,function(x) {
        sd(x)/sqrt(length(x))
    } )
    
    barx = 1:length(m)
    
    if( is.null(ylim) ) ylim=range( c( m-se, m+se) )
    
    if(line) {
        plot( m, type="l", ylim=ylim, axes=FALSE, ylab=ylab )
    } else {
        plot( m, type="l", ylim=ylim, axes=FALSE, ylab=ylab, col='white' )
    }
    
    axis(side=1,at=barx,label=names(X),las=3)
    axis(side=2)
    points( 1:length(m), m, pch=16, cex=1.5 )
    arrows(barx , m+se, barx , m, angle=90, code=3, length=0.05)
    arrows(barx , m-se, barx , m, angle=90, code=3, length=0.05)
}


leaf.plot = function(dist, xlim=NULL, ylim=NULL, upper=1:length(dist), lower=NULL, 
                     axis=TRUE, title=NULL, ref_title=NULL, main="leaf.plot")
{
    
    plot_dist <- function(d, col, size=0.5, title=NULL) {
        d = sort(d,decreasing=T)
        angles <- seq(0, pi, length.out = length(d))
        x = d*cos(angles)
        y = d*sin(angles)
        points(x, y, col=col, pch=16, cex=size)
        if( !is.null(title) )
        {
            i = which.max(abs(y))
            text(x[i], y[i], title)
        }
    }
    
    get_dist <- function(d) {
        d = sort(d,decreasing=T)
        angles <- seq(0, pi, length.out = length(d))
        x = d*cos(angles)
        y = d*sin(angles)
        data.frame(x,y)
    }
    
    all_dist = sort(do.call(c,dist))
    n = length(all_dist)
    
    hues = seq(15, 375, length = n + 1)
    colors <- grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
    all_cols = rev( colors )
    
    values = lapply(dist, get_dist)
    if(is.null(xlim)) xlim = range( do.call(c,lapply(values,function(z)z$x)) )
    if(is.null(ylim)) {
        ylim1 = range( do.call(c,lapply(values[upper],function(z) z$y )) )
        ylim2 = range( do.call(c,lapply(values[lower],function(z) -z$y )) )
        ylim = range(c(ylim1,ylim2))
    }
    #plot(NA, xlim=xlim, ylim=ylim, xlab="", ylab="")
    plot(NA, axes=FALSE, frame.plot=FALSE, xlim=xlim, ylim=ylim, xlab="",ylab="", main=main)
    tmp = sapply(length(dist):1, function(i) {
        d = sort(dist[[i]],decreasing=T)
        col = all_cols[ match(d, all_dist) ]
        if( i %in% upper ) plot_dist (d, col, s=0.5, title=title[i])
        if( i %in% lower ) plot_dist (-d, rev(col), s=0.5, title=title[i])
    } )
    
    if(axis) {
        points( all_dist, rep(0,length(all_dist)), col=all_cols, cex=1, pch=15 )
        x = seq(0, max(all_dist), length=10)
        x = x[-length(x)]
        y = rep(0,length(x))
        points( x, y, pch=3, cex=1 )
        arrows( 0, 0, max(all_dist), 0 )
    } else {
        points( all_dist, rep( min(ylim)-0.1,length(all_dist)), col=all_cols, cex=1, pch=15 )
    }
    
    points(0,0,pch=16,col=all_cols[3000],cex=3)
    if( !is.null(ref_title) ) text(0,0,ref_title)
    
}

barplot1 <- function(score, colors, celltypes)
{
    colors = col2rgb(colors)
    colors = apply(colors,2,function(x) rgb(x[1],x[2],x[3],max=255,alpha=0.7*255) )
    
    t2f = function(t) {
            as.data.frame( do.call(rbind,lapply(1:nrow(t),function(i)t[i,])), 
                        row.names=rownames(t) )
    }
      
    celltypes_on_spot2 = lapply(score, function(x) {
        y = t2f(x$celltypes_on_spot)
        missed = setdiff(celltypes, colnames(y)) 
        if(length(missed)>0) {
            for(i in 1:length(missed))
            {
                y = data.frame(0,y)
                colnames(y)[1] = missed[i]
            }
        }
        y = y[,match(celltypes,colnames(y))]
        y
    })
    
    data <- do.call(cbind, lapply(celltypes_on_spot2,colMeans ) )
    #data = data[order(rownames(data)),]
    
    m <- sapply(celltypes_on_spot2,function(x) mean(rowSums(x)) )
    se <- sapply(celltypes_on_spot2,function(y) sd(rowSums(y))/sqrt(nrow(y)) )
    
    # Get the stacked barplot
    barx <- barplot(data, 
                    col=colors , 
                    border="white", 
                    space=0.04, 
                    font.axis=2, 
                    xlab="spot size",
                    las=3,
                    axes=FALSE,
                    ylim=c(-10,200))
    
    lines( barx, m )
    arrows(barx , m+se, barx , m, angle=90, code=3, length=0.05)
    
    labels = seq(0,150,50)
    axis(side=2, at=labels, labels=labels )
    axis(side=1, at=barx )
    
}



barplot2 <- function(score, colors)
{
    colors = col2rgb(colors)
    colors = apply(colors,2,function(x) rgb(x[1],x[2],x[3],max=255,alpha=0.7*255) )
    
    data <- do.call(cbind, lapply(score,function(x) colMeans(x$celltypes_on_spot)) )
    data = data[order(rownames(data)),]
    
    m <- sapply(score,function(x) mean(rowSums(x$celltypes_on_spot)) )
    se <- sapply(score,function(x) {
        y = x$celltypes_on_spot
        sd(rowSums(y))/sqrt(nrow(y))
    } )
    
    # Get the stacked barplot
    barx <- barplot(data, 
                    col=colors , 
                    border="white", 
                    space=0.04, 
                    font.axis=2, 
                    xlab="spot size",
                    las=3,
                    axes=FALSE,
                    ylim=c(-10,200))
    
    lines( barx, m )
    arrows(barx , m+se, barx , m, angle=90, code=3, length=0.05)
    
    labels = seq(0,150,50)
    axis(side=4, at=labels, labels=labels )
    
    ######################################################################
    #f = function(x) (x-min(x))/(max(x)-min(x))
    f = function(x) ( x-min(rmse)+min(rmse_se) ) / ( max(rmse)-min(rmse)+2*min(rmse_se) ) * ( max(m)-min(m) )
    
    rmse = sapply(score,function(x) mean(sqrt(x$mse)))
    rmse_se = sapply(score,function(x) {
        y = sqrt(x$mse)
        sd(y)/sqrt(length(y))
    })
    rmse2 = ( rmse-min(rmse) ) / ( max(rmse)-min(rmse) ) * ( max(m)-min(m) )
    rmse_se1 = ( rmse-rmse_se-min(rmse) ) / ( max(rmse)-min(rmse) ) * ( max(m)-min(m) )
    rmse_se2 = ( rmse+rmse_se-min(rmse) ) / ( max(rmse)-min(rmse) ) * ( max(m)-min(m) )
    
    rmse2 = f(rmse)
    rmse_se1 = f(rmse-rmse_se)
    rmse_se2 = f(rmse+rmse_se)
    
    points(barx,rmse2,pch=16, col='darkred')
    lines( barx, rmse2, col='darkred')
    arrows(barx , rmse_se1, barx , rmse_se2, angle=90, code=3, length=0.05, col='darkred')
    
    labels=seq(floor(min(rmse)),ceiling(min(rmse)),0.2)
    axis(side=2, at=f(labels), labels=labels, col='darkred' )
    
}


ST_scatterpie_overlay <- function( prop_data, coordinates, x_scale, cell_types_all, 
                                   img_path=NULL, scatterpie_alpha=1, pie_scale = 0.8, colors=NULL )
{
    if(0)
    {
        slice <- names(se_obj@images)[1]
        x_scale = se_obj@images[[slice]]@scale.factors$lowres
        coordinates = data.frame(se_obj@images[[slice]]@coordinates)
        prop_data = se_obj@meta.data
    }
    
    prop_data = data.frame(barcodeID=rownames(prop_data), prop_data)
    coordinates = data.frame(barcodeID=rownames(coordinates), coordinates)
    spatial_coord = merge(coordinates, prop_data)
    spatial_coord$imagerow_scaled = spatial_coord$imagerow*x_scale
    spatial_coord$imagecol_scaled = spatial_coord$imagecol*x_scale
    
    if( !is.null(img_path) ) {
        
        img_frmt <- base::tolower(stringr::str_sub(img_path, -4, -1))
        
        if (img_frmt %in% c(".jpg", "jpeg")) {
            img <- jpeg::readJPEG(img_path)
        } else if (img_frmt == ".png") {
            img <- png::readPNG(img_path)
        }
        
        img_grob <- grid::rasterGrob(img, interpolate = FALSE, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
        scatterpie_plt <- suppressMessages(ggplot2::ggplot() + 
                                               ggplot2::annotation_custom(grob = img_grob, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) + 
                                               scatterpie::geom_scatterpie(data = spatial_coord, 
                                                                           ggplot2::aes(x = imagecol_scaled, y = imagerow_scaled), 
                                                                           cols = cell_types_all, color = NA, 
                                                                           alpha = scatterpie_alpha, pie_scale = pie_scale) ) 
        if(!is.null(colors)) {
            scatterpie_plt <- suppressMessages(scatterpie_plt + 
                                                   scale_fill_manual(values=colors ) +
                                                   scale_x_continuous(expand=c(0,0), lim=c(0,3)) +
                                                   scale_y_continuous(expand=c(0,0), lim=c(0,6)) )
        }
        
        scatterpie_plt <- suppressMessages( scatterpie_plt + 
                                                ggplot2::scale_y_reverse() + 
                                                ggplot2::ylim(nrow(img), 0) + 
                                                ggplot2::xlim(0, ncol(img)) + 
                                                cowplot::theme_half_open(11, rel_small = 1) + 
                                                ggplot2::theme_void() + 
                                                ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on"))
        
        
    } else {
        
        spatial_coord$imagerow_adjusted = max(spatial_coord$imagerow) + min(spatial_coord$imagerow) - spatial_coord$imagerow
        spatial_coord$imagecol_adjusted = max(spatial_coord$imagecol) + min(spatial_coord$imagecol) - spatial_coord$imagecol
        
        scatterpie_plt <- suppressMessages(ggplot2::ggplot() + 
                                               scatterpie::geom_scatterpie(data = spatial_coord,
                                                                           ggplot2::aes(x = imagerow_adjusted, y = imagecol_adjusted), 
                                                                           cols = cell_types_all, color = NA, # if NA, then colors set by default
                                                                           alpha = scatterpie_alpha, pie_scale = pie_scale) )
        if(!is.null(colors)) {
            scatterpie_plt <- suppressMessages(scatterpie_plt + 
                                                   scale_fill_manual(values=colors ) +
                                                   scale_x_continuous(expand=c(0,0), lim=c(0,3)) +
                                                   scale_y_continuous(expand=c(0,0), lim=c(0,6)) )
        }
        
        scatterpie_plt <- suppressMessages( scatterpie_plt + 
                                                ggplot2::scale_y_reverse() + 
                                                ggplot2::xlim( range(spatial_coord$imagerow_adjusted)*c(0.99, 1.01) ) + 
                                                ggplot2::ylim( range(spatial_coord$imagecol_adjusted)*c(0.99, 1.01) ) + 
                                                cowplot::theme_half_open(11, rel_small = 1) + 
                                                ggplot2::theme_void() + 
                                                ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") )
    }
    
    scatterpie_plt
    
}


BayesSpaceVisium_scatterpie_overlay <- function( prop_data, coordinates, x_scale, cell_types_all, 
                                   img_path=NULL, scatterpie_alpha=1, pie_scale = 0.8, colors=NULL )
{
    prop_data = data.frame(barcodeID=rownames(prop_data), prop_data)
    coordinates = data.frame(barcodeID=rownames(coordinates), coordinates)
    spatial_coord = merge(coordinates, prop_data)
    
    rr = max(coordinates$imagerow) - min(coordinates$imagerow)
    cr = max(coordinates$imagecol) - min(coordinates$imagecol)
    
    fs = function(x) (x-min(x))/(max(x)-min(x))
    spatial_coord$imagerow_scaled = fs(spatial_coord$imagerow)
    spatial_coord$imagecol_scaled = fs(spatial_coord$imagecol)
    spatial_coord = spatial_coord[order(as.character(spatial_coord$barcodeID)), ]
    
    spatial_coord_1 = spatial_coord[ grep("[.]1",as.character(spatial_coord$barcodeID)), ]
    distMat <- dist(spatial_coord_1[, c("imagerow_scaled", "imagecol_scaled")])
    distMat <- as.matrix(distMat)
    diag(distMat) <- Inf
    dr <- mean(apply(distMat, 2, min))
    
    group = gsub( "[.][0-9]", "", as.character(spatial_coord$barcodeID) )
    spatial_coord_list = split(spatial_coord, group)
    
    spatial_coord = do.call(rbind, lapply( spatial_coord_list, function(x) {
        xr = x$imagerow_scaled
        xc = x$imagecol_scaled
        cr = mean(xr)
        cc = mean(xc)
        
        dx <- dist(x[, c("imagerow_scaled", "imagecol_scaled")])
        ratio = dr/max(dx)/2*1.5
        
        x$imagerow_scaled = cr + ratio*(xr - cr)
        x$imagecol_scaled = cc + ratio*(xc - cc)
        x
    } ))
    
    with(spatial_coord, plot( imagerow_scaled, imagecol_scaled ))
    
    scatterpie_plt <- suppressMessages(ggplot2::ggplot() + 
                                           scatterpie::geom_scatterpie(data = spatial_coord,
                                                                       ggplot2::aes(x = imagerow_scaled, y = imagecol_scaled), 
                                                                       cols = cell_types_all, color = NA, # if NA, then colors set by default
                                                                       alpha = scatterpie_alpha, pie_scale = pie_scale) )
    if(!is.null(colors)) {
        scatterpie_plt <- suppressMessages(scatterpie_plt + 
                                               scale_fill_manual(values=colors ) +
                                               scale_x_continuous(expand=c(0,0), lim=c(0,3)) +
                                               scale_y_continuous(expand=c(0,0), lim=c(0,6)) )
    }
    
    scatterpie_plt <- suppressMessages( scatterpie_plt + 
                                            ggplot2::scale_y_reverse() + 
                                            ggplot2::xlim( range(spatial_coord$imagerow_scaled)+c(-0.1, 0.1) ) + 
                                            ggplot2::ylim( range(spatial_coord$imagecol_scaled)+c(-0.1, 0.1) ) + 
                                            cowplot::theme_half_open(11, rel_small = 1) + 
                                            ggplot2::theme_void() + 
                                            ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") )

    
    scatterpie_plt
    
}


dwls = function(S, B)
{
    library(Matrix)
    source("/archive/SCCC/Hoshida_lab/s184554/Project/stRNAseq/DWLS/Deconvolution_functions.R")
    
    nb = ncol(B)
    ns = ncol(S)
    prop_mat <- as.data.frame(matrix(0, nrow = nb, ncol = ns) )
    rownames(prop_mat) <- colnames(B)
    colnames(prop_mat) <- colnames(S)
    for (ib in 1:nb) {
        print( sprintf("Estimating proportion for spot : %d / %d", ib, nb) )
        b <- B[,ib]
        tr <- trimData(S, b)
        tr$sig <- tr$sig[,colSums(tr$sig) > 0]
        is_pd <- eigen(t(tr$sig)%*%tr$sig)$values
        is_pd <- all(is_pd > 10e-6)
        if (!(is_pd)) { next}
        try(solDWLS <- solveDampenedWLS(tr$sig,tr$bulk),next)
        print("Proportions >> ")
        prop_mat[ib,names(solDWLS)] <- solDWLS
    }
    
    t(prop_mat)
    
}


my_col2rgb = function(col, transparency)
{
    color = col2rgb(col)
    rgb(color[1],color[2],color[3],max=255,alpha=transparency*255)
}





