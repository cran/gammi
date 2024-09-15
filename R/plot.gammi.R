plot.gammi <-
  function(x, terms = x$term.labels, conf.int = TRUE, n = 400, 
           intercept = FALSE, random = TRUE, ask = dev.interactive(), 
           xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, ...){
    # plot method for gammi objects
    # Nathaniel E. Helwig (helwig@umn.edu)
    # updated: 2024-07-22
    
    
    
    #########***#########   INITIAL CHECKS   #########***#########
    
    ### check x
    if(!inherits(x, "gammi")) stop("Input 'x' must be of class 'gammi'")
    if(is.null(x$formula)) stop("Plotting method not available for gammi.default objects")
    
    ### check terms
    if(is.null(terms)){
      terms <- x$term.labels
      tid <- 1:length(x$terms)
    } else {
      terms <- as.character(terms)
      tid <- match(terms, x$term.labels)
      natid <- is.na(tid)
      if(any(natid)) stop("Input 'terms' contains terms not included in the model:\n  ", paste(terms[which(natid)], collapse = ", "))
    }
    nterms0 <- length(x$term.labels)
    
    ### check conf.int
    conf.int <- as.logical(conf.int[1])
    if(!any(conf.int == c(TRUE, FALSE))) stop("Input 'conf.int' must be TRUE or FALSE")
    
    ### check n
    n <- as.integer(n[1])
    if(n <= 3) stop("Inpit 'n' must be a integer greater than or equal to 3.")
    
    ### check intercept
    intercept <- as.logical(intercept[1])
    if(!any(intercept == c(TRUE, FALSE))) stop("Input 'intercept' must be TRUE or FALSE")
    int <- ifelse(intercept, x$coefficients[1], 0)
    
    ### reset par on exit
    oldplt <- par()$plt
    oldnew <- par()$new
    on.exit(par(plt = oldplt, new = oldnew))
    
    ### effect type (main or interaction?)
    nterms <- length(terms)
    effect <- rep(NA, nterms)
    effvar <- vector("list", nterms)
    for(i in 1:nterms){
      effvar[[i]] <- strsplit(terms[i], split = ":")[[1]]
      effect[i] <- length(effvar[[i]])
    } # end for(i in 1:nterms)
    
    ### add fixed formula
    if(!is.null(x$fixed)){
      fixform <- strsplit(as.character(x$fix), "~", fixed = T)[[2]]
      ranform <- strsplit(as.character(x$formula), "~", fixed = T)
      newform <- as.formula(paste0(ranform[[2]], " ~ ", fixform, " + ", ranform[[3]]))
    } else {
      newform <- x$formula
    }
    
    ### variable names, classes, and ranges
    xnames <- all.vars(newform)[-1]
    nvars <- length(xnames)
    xclass <- rep(NA, nvars)
    xrange <- matrix(0, 2, nvars)
    xlevel <- vector("list", nvars)
    xisfac <- rep(FALSE, nvars)
    for(j in 1:nvars) {
      xclass[j] <- class(x$data[,xnames[j]])[1]
      if(xclass[j] %in% c("factor", "character")){
        xisfac[j] <- TRUE
        xlevel[[j]] <- levels(as.factor(x$data[,xnames[j]]))
        xrange[,j] <- c(1L, length(xlevel[[j]]))
      } else {
        xrange[,j] <- range(x$data[,xnames[j]])
      }
    }
    
    ### interactive plot setup
    total.num.terms <- nterms + (random && !is.null(x$random.coefficients))
    if(total.num.terms > 1L && ask){
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    
    ### response name
    yname <- all.vars(x$formula)[1]
    
    
    
    #########***#########   PLOT FIXED/SMOOTH EFFECTS   #########***#########
    
    ### plotting...
    for(i in 1:nterms){
      
      ## main effect
      if(effect[i] == 1L){
        
        varid <- match(effvar[[i]], xnames)
        if(xisfac[varid]){
          
          newdata <- x$data[1:length(xlevel[[i]]),,drop=FALSE]
          newseq <- xlevel[[i]]
          newdata[,effvar[[i]]] <- newseq
          fit <- predict(x, newdata = newdata, type = "terms", conf.int = conf.int)
          dev.hold()
          if(conf.int){
            visualizer1(x = as.integer(factor(newseq, levels = newseq)), 
                        y = fit$fit[,effvar[[i]]] + int, 
                        bars = TRUE, 
                        lwr = fit$lwr[,effvar[[i]]] + int, 
                        upr = fit$upr[,effvar[[i]]] + int,
                        xlab = if(is.null(xlab)) effvar[[i]] else xlab,
                        ylab = if(is.null(ylab)) yname else ylab, 
                        main = if(is.null(main)) paste(effvar[[i]], "effect") else main,
                        axes = FALSE, ...)
            axis(1, at = 1:length(xlevel[[i]]), labels = xlevel[[i]], ...)
            axis(2, ...)
          } else {
            visualizer1(x = as.integer(factor(newseq, levels = newseq)), 
                        y = fit[,effvar[[i]]] + int,
                        bars = TRUE, 
                        xlab = if(is.null(xlab)) effvar[[i]] else xlab,
                        ylab = if(is.null(ylab)) yname else ylab, 
                        main = if(is.null(main)) paste(effvar[[i]], "effect") else main,
                        axes = FALSE, ...)
            axis(1, at = 1:length(xlevel[[i]]), labels = xlevel[[i]], ...)
            axis(2, ...)
          }
          box()
          dev.flush()
          
        } else {
          
          newdata <- x$data[1:n,]
          newseq <- seq(xrange[1,varid], xrange[2,varid], length.out = n)
          newdata[,effvar[[i]]] <- newseq
          fit <- predict(x, newdata = newdata, type = "terms", conf.int = conf.int)
          dev.hold()
          if(conf.int){
            visualizer1(x = newseq, 
                        y = fit$fit[,effvar[[i]]] + int, 
                        lwr = fit$lwr[,effvar[[i]]] + int, 
                        upr = fit$upr[,effvar[[i]]] + int,
                        xlab = if(is.null(xlab)) effvar[[i]] else xlab,
                        ylab = if(is.null(ylab)) yname else ylab, 
                        main = if(is.null(main)) paste(effvar[[i]], "effect") else main, ...)
          } else {
            visualizer1(x = newseq, 
                        y = fit[,effvar[[i]]] + int,
                        xlab = if(is.null(xlab)) effvar[[i]] else xlab,
                        ylab = if(is.null(ylab)) yname else ylab, 
                        main = paste(effvar[[i]], "effect"), ...)
          }
          dev.flush()
          
        } # end if(xisfac[varid])
        
      } else if(effect[i] == 2L){
        
        # which variables?
        varid <- match(effvar[[i]], xnames)
        
        # which types of variables?
        if(xisfac[varid[1]] & xisfac[varid[2]]){
          nlev1 <- length(xlevel[[varid[1]]])
          nlev2 <- length(xlevel[[varid[2]]])
          newdata <- x$data[1:(nlev1 * nlev2),]
          newdata[,xnames[varid]] <- expand.grid(xlevel[[varid[1]]], xlevel[[varid[2]]])
          fit <- predict(x, newdata = newdata, type = "terms")[,paste(effvar[[i]], collapse = ":")]
          dev.hold()
          visualizer2(x = 1:nlev1, y = 1:nlev2,
                      z = matrix(fit, nlev1, nlev2), 
                      xticks = 1:nlev1, xlabels = xlevel[[varid[1]]],
                      yticks = 1:nlev2, ylabels = xlevel[[varid[2]]],
                      xlim = c(0.5, nlev1 + 0.5), ylim = c(0.5, nlev2 + 0.5),
                      xlab = if(is.null(xlab)) effvar[[i]][1] else xlab,
                      ylab = if(is.null(ylab)) effvar[[i]][2] else ylab, 
                      zlab = if(is.null(zlab)) yname else zlab, 
                      main = if(is.null(main)) paste(paste(effvar[[i]], collapse = ":"), "effect") else main,
                      ...)
          dev.flush()
        } else if(xisfac[varid[1]] & !xisfac[varid[2]]){
          
          nlev1 <- length(xlevel[[varid[1]]])
          for(j in 1:nlev1){
            
            effname <- paste(effvar[[i]], collapse = ":")
            newdata <- x$data[1:n,]
            newdata[,effvar[[i]][1]] <- xlevel[[varid[1]]][j]
            newseq <- seq(xrange[1,varid[2]], xrange[2,varid[2]], length.out = n)
            newdata[,effvar[[i]][2]] <- newseq
            fit <- predict(x, newdata = newdata, type = "terms", conf.int = conf.int)
            dev.hold()
            if(conf.int){
              visualizer1(x = newseq, 
                          y = fit$fit[,effname] + int, 
                          lwr = fit$lwr[,effname] + int, 
                          upr = fit$upr[,effname] + int,
                          xlab = if(is.null(xlab)) effvar[[i]][2] else xlab,
                          ylab = if(is.null(ylab)) yname else ylab, 
                          main = if(is.null(main)) paste(effname, "effect") else main, ...)
            } else {
              visualizer1(x = newseq, 
                          y = fit[,effname] + int,
                          xlab = if(is.null(xlab)) effvar[[i]][2] else xlab,
                          ylab = if(is.null(ylab)) yname else ylab, 
                          main = if(is.null(main)) paste(effname, "effect") else main, ...)
            }
            title(paste0(effvar[[i]][1], " = ", xlevel[[varid[1]]][j]), line = 0.5)
            dev.flush()
            
          } # end for(j in 1:nlev1)
          
          
        } else if(!xisfac[varid[1]] & xisfac[varid[2]]){
          
          nlev2 <- length(xlevel[[varid[2]]])
          for(j in 1:nlev2){
            
            effname <- paste(effvar[[i]], collapse = ":")
            newdata <- x$data[1:n,]
            newseq <- seq(xrange[1,varid[1]], xrange[2,varid[1]], length.out = n)
            newdata[,effvar[[i]][1]] <- newseq
            newdata[,effvar[[i]][2]] <- xlevel[[varid[2]]][j]
            fit <- predict(x, newdata = newdata, type = "terms", conf.int = conf.int)
            dev.hold()
            if(conf.int){
              visualizer1(x = newseq, 
                          y = fit$fit[,effname] + int, 
                          lwr = fit$lwr[,effname] + int, 
                          upr = fit$upr[,effname] + int,
                          xlab = if(is.null(xlab)) effvar[[i]][1] else xlab,
                          ylab = if(is.null(ylab)) yname else ylab, 
                          main = if(is.null(main)) paste(effname, "effect") else main, ...)
            } else {
              visualizer1(x = newseq, 
                          y = fit[,effname] + int,
                          xlab = if(is.null(xlab)) effvar[[i]][1] else xlab,
                          ylab = if(is.null(ylab)) yname else ylab, 
                          main = if(is.null(main)) paste(effname, "effect") else main, ...)
            }
            title(paste0(effvar[[i]][2], " = ", xlevel[[varid[2]]][j]), line = 0.5)
            dev.flush()
            
          } # end for(j in 1:nlev2)
          
        } else {
          nsqrt <- ceiling(sqrt(n))
          newdata <- x$data[1:(nsqrt^2),]
          xseq1 <- seq(xrange[1,varid[1]], xrange[2,varid[1]], length.out = nsqrt)
          xseq2 <- seq(xrange[1,varid[2]], xrange[2,varid[2]], length.out = nsqrt)
          newdata[,xnames[varid]] <- expand.grid(xseq1, xseq2)
          fit <- predict(x, newdata = newdata, type = "terms")[,paste(effvar[[i]], collapse = ":")]
          dev.hold()
          visualizer2(x = xseq1, y = xseq2,
                      z = matrix(fit, nsqrt, nsqrt), 
                      xlab = if(is.null(xlab)) effvar[[i]][1] else xlab,
                      ylab = if(is.null(ylab)) effvar[[i]][2] else ylab, 
                      zlab = if(is.null(zlab)) yname else zlab, 
                      main = if(is.null(main)) paste(paste(effvar[[i]], collapse = ":"), "effect") else main)
          dev.flush()
        } # end if(xisfac[varid[1]] & xisfac[varid[2]])
        
        
      } else {
        
        warning("Three-way and higher order interactions are not supported.")
        
      } # end if(effect[i] == 1L)
      
    } # end for(i in 1:nterms)
    
    
    
    #########***#########   RANDOM EFFECTS   #########***#########
    
    if(!is.null(x$random.coefficients) && random){
      
      rnames <- names(x$random.coefficients)
      for(i in 1:length(x$random.coefficients)){
        if(ncol(x$random.coefficients[[i]]) == 1L){
          dev.hold()
          qqnorm(unlist(x$random.coefficients[[i]]), ...)
          title(main = rnames[i], line = 0.5)
          qqline(unlist(x$random.coefficients[[i]]))
          dev.flush()
        } else {
          inames <- colnames(x$random.coefficients)
          for(j in 1:ncol(x$random.coefficients[[i]])){
            dev.hold()
            qqnorm(x$random.coefficients[[i]][,j], ...)
            title(main = paste0(rnames[i], ".", inames[j]), line = 0.5)
            qqline()
            dev.flush(x$random.coefficients[[i]][,j])
          }
        } # end if(ncol(x$random.coefficients) == 1L)
      } # end for(i in 1:length(x$random.coefficients))
      
    } # end if(!is.null(x$random.coefficients))
    
    
  } # end plot.gammi