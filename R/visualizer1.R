visualizer1 <-
  function(x, y, bars = FALSE, bw = 0.02, lty = 1, lwd = 2, col = "black", 
           lwr = NULL, upr = NULL, ci.lty = 2, ci.lwd = 1.25, ci.col = "black",
           zero = TRUE, zero.lty = 3, xlim = NULL, ylim = NULL, 
           xlab = NULL, ylab = NULL, main = NULL, ...){
    
    n <- length(x)
    y <- as.numeric(y)
    if(length(y) != n) stop("Inputs 'x' and 'y' must have the same length")
    if(!is.null(lwr) && length(lwr) != n) stop("Inputs 'x' and 'lwr' must have the same length")
    if(!is.null(upr) && length(upr) != n) stop("Inputs 'x' and 'upr' must have the same length")
    if(is.null(ylim)) ylim <- extendrange(c(y, lwr, upr))
    
    if(bars){
      
      xr <- range(x)
      xr <- xr[2] - xr[1]
      plot(x, y, type = "n", pch = 19, col = col, 
           xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
      if(zero) abline(h = 0, lty = zero.lty)
      if(!is.null(lwr) & !is.null(upr)){
        for(i in 1:n){
          segments(x0 = x[i], y0 = lwr[i], y1 = upr[i], lwd = lwd, lty = lty, col = ci.col)
          segments(x0 = x[i] - xr * bw, y0 = lwr[i], x1 = x[i] + xr * bw, lwd = lwd, lty = lty, col = ci.col)
          segments(x0 = x[i] - xr * bw, y0 = upr[i], x1 = x[i] + xr * bw, lwd = lwd, lty = lty, col = ci.col)
        }
        points(x, y, pch = 19, col = col)
      }
      
    } else {
      
      plot(x, y, type = "n", lty = lty, lwd = lwd, col = col, 
           xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
      if(zero) abline(h = 0, lty = zero.lty)
      if(!is.null(lwr)) lines(x, lwr, lty = ci.lty, lwd = ci.lwd, col = ci.col)
      if(!is.null(upr)) lines(x, upr, lty = ci.lty, lwd = ci.lwd, col = ci.col)
      lines(x, y, lty = lty, lwd = lwd, col = col)
      
    } # end if(bars)
    
    
  } # end visualizer1