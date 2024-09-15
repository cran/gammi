summary.gammi <-
  function(object, ...){
    # summary method for gammi objects
    # Nathaniel E. Helwig (helwig@umn.edu)
    # updated: 2024-07-17
    
    
    ### initial check
    if(!inherits(object, "gammi")) stop("Input 'object' should be of class 'gammi'.")
    
    ### get group info
    group <- as.factor(object$group)
    ngrps <- nlevels(group)
    gsize <- table(group)
    group <- as.integer(group)
    
    ### make x?
    if(is.null(object$x)){
      if(is.null(object$fixed)){
        object$x <- spline.model.matrix(object$formula, data = object$data,
                                        knots = object$spline.info$knots,
                                        m = object$spline.info$m,
                                        periodic = object$spline.info$periodic,
                                        xlev = object$spline.info$xlev)[,-1,drop=FALSE]
      } else {
        x0 <- model.matrix(object$fixed, data = object$data, contrasts.arg = object$contrasts)
        x1 <- spline.model.matrix(object$formula, data = object$data,
                                  knots = object$spline.info$knots,
                                  m = object$spline.info$m,
                                  periodic = object$spline.info$periodic,
                                  xlev = object$spline.info$xlev)
        object$x <- cbind(x0[,-1], x1[,-1])
      }
    }
    
    ### get edf for each term
    edf <- rep(0, ngrps + 1)
    names(edf) <- c("(Intercept)", object$term.labels)
    edf[1] <- 1
    for(k in 1:ngrps){
      id <- which(group == k)
      xin <- object$x[,id,drop=FALSE]
      xout <- cbind(1, object$x[,-id,drop=FALSE])
      xsvd <- svd(xout, nv = 0)
      xproj <- xin - xsvd$u %*% crossprod(xsvd$u, xin)
      edf[k+1] <- sum((xproj %*% (object$vcovchol[id+1,,drop=FALSE]))^2) / object$dispersion
    }
    
    ### test statistic
    Tstat <- rep(0, ngrps+1)
    Tstat[1] <- object$coefficients[1]^2 / sum(object$vcovchol[1,]^2)
    for(k in 1:ngrps){
      id <- which(group == k)
      temp <- svd(object$vcovchol[id+1,,drop=F])
      temp <- temp$v %*% diag(1/temp$d, nrow = length(id), ncol = length(id)) %*% t(temp$u)
      Tstat[k+1] <- sum(( temp %*% object$coefficients[id+1] )^2)
      #Tstat[k+1] <- t(object$coefficients[id+1]) %*% solve(tcrossprod(object$vcovchol[id+1,,drop=F])) %*% object$coefficients[id+1]
    }
    
    ### coefficients table
    if(object$family$family %in% c("gaussian", "Gamma", "inverse.gaussian")){
      SS <- Tstat * object$dispersion
      Tstat <- Tstat / edf
      pvalue <- 1 - pf(Tstat, edf, object$df.residual)
      tablenames <- c("")
      statname <- "F value"
      pvalname <- "Pr(>F)"
    } else {
      SS <- Tstat
      pvalue <- 1 - pchisq(Tstat, edf)
      statname <- "Chi value"
      pvalname <- "Pr(>Chi)"
    }
    
    ### table results
    coeftab <- cbind(edf, SS, Tstat, pvalue)
    colnames(coeftab) <- c("edf", "Sum Sq", statname, pvalname)
    
    ### calculate deviance residuals
    dev.resid <- object$family$dev.resids(object$mer@resp$y, object$fitted.values, object$mer@resp$weights)
    dev.resid <- sign(object$mer@resp$y - object$fitted.values) * sqrt(dev.resid)
    
    ### calculate importance
    object$formula <- NULL
    fit <- scale(predict(object, newx = object$x, type = "terms"), scale = FALSE)
    fit0 <- rowSums(fit)
    ss0 <- sum(fit0^2)
    imp <- rep(0, ncol(fit))
    names(imp) <- object$term.labels
    for(k in 1:ncol(fit)){
      imp[k] <- sum(fit[,k] * fit0) / ss0
    }
    
    ### calculate variance inflation factor
    Cmat <- crossprod(fit) / tcrossprod(sqrt(colSums(fit^2)))
    vif <- diag(solve(Cmat))
    
    ### collect results
    res <- list(call = object$call,
                term.labels = object$term.labels,
                family = object$family,
                logLik = object$logLik,
                aic = object$aic,
                deviance = object$deviance,
                deviance.resid = dev.resid,
                r.squared = object$r.squared,
                
                df = object$edf + object$df.random,
                significance = coeftab,
                importance = imp * 100,
                vif = vif)
    class(res) <- "summary.gammi"
    return(res)
    
  } # end summary.gammi


print.summary.gammi <-
  function(x, digits = max(3, getOption("digits") - 3),
           symbolic.cor = x$symbolic.cor,
           signif.stars = getOption("show.signif.stars"),
           show.residuals = FALSE, ...){
    
    # set digits
    oldoptions <- options(digits = digits)
    on.exit(options(oldoptions))
    
    # print call
    cat("\nCall:\n")
    print(x$call)
    
    # print residuals
    cat("\nResiduals:\n")
    rsum <- summary(x$deviance.resid)[-4]
    names(rsum) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(rsum)
    
    # check if significance codes are needed
    if(signif.stars){
      pcodes <- rev(c("***", "**", "*", ".", ""))
    }
    
    # smooth effect table
    cat("\nSignificance:\n")
    smoo <- as.data.frame(x$significance)
    if(signif.stars){
      pstars <- symnum(x$significance[,4],
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " "), na = NA)
      smoonames <- names(smoo)
      smoo <- cbind(smoo, pstars)
      colnames(smoo) <- c(smoonames, "")
    }
    print(as.matrix(smoo), na.print = "", quote = FALSE, right = TRUE)
    if(signif.stars){
      cat("---\n")
      cat("Signif. codes: ", attr(pstars, "legend"), "\n")
    }
    
    # diagnostics
    cat("\nDiagnostics:\n")
    print(cbind(edf = x$significance[-1,1], 
                imp = x$importance,
                vif = x$vif))
    cat("\n")
    
  } # end print.summary.gammi