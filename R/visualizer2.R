visualizer2 <-
  function(x, y, z, col = NULL, ncolor = 21,
           xlim = NULL, ylim = NULL, zlim = NULL, zline = 1.5,
           xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, 
           xticks = NULL, xlabels = NULL, yticks = NULL, ylabels = NULL, ...){
    
    oldplt <- par()$plt
    oldnew <- par()$new
    on.exit(par(plt = oldplt, new = oldnew))
    
    if(is.null(col)){
      col <- c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
               "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
    }
    col <- colorRampPalette(col)(ncolor)
    
    z <- matrix(z, length(x), length(y))
    if(is.null(xlim)) xlim <- range(x)
    if(is.null(ylim)) ylim <- range(y)
    if(is.null(zlim)) zlim <- range(z)
    if(is.null(xlab)) xlab <- "x"
    if(is.null(ylab)) ylab <- "y"
    if(is.null(zlab)) zlab <- "z"
    
    par(plt = c(0.85, 0.9, 0.15, 0.85))
    zseq <- seq(zlim[1], zlim[2], length.out = ncolor)
    graphics::image(x = 1, y = zseq,
                    z = matrix(zseq, nrow = 1, ncol = ncolor), 
                    xlab = "", ylab = "", col = col, axes = FALSE)
    ticks <- pretty(zlim)
    labels <- ticks
    midpoint <- floor(median(1:length(ticks)))
    labels[midpoint] <- NA
    axis(4, at = ticks, labels = labels)
    box()
    mtext(zlab, side = 4, line = zline)
    
    par(plt = c(0.15, 0.8, 0.15, 0.85), new = TRUE)
    graphics::image(x = x, y = y, z = z, col = col, 
                    xlim = xlim, ylim = ylim, zlim = zlim,
                    xlab = xlab, ylab = ylab, main = main, 
                    axes = FALSE, ...)
    axis(1, at = xticks, labels = xlabels)
    axis(2, at = yticks, labels = ylabels)
    box()
    
    par(plt = oldplt, new = oldnew)
    
  } # end visualizer2