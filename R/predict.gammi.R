predict.gammi <-
  function(object, 
           newx,
           newdata, 
           se.fit = FALSE,
           type = c("link", "response", "terms"),
           conf.int = FALSE, 
           conf.level = 0.95,
           ...){
    # predict method for gammi objects
    # Nathaniel E. Helwig (helwig@umn.edu)
    # updated: 2024-07-17
    
    
    ### check object
    if(!inherits(object, "gammi")) stop("Input 'object' must be of class 'gammi'")
    
    ### default or formula?
    default <- ifelse(is.null(object$formula), TRUE, FALSE)
    
    ### check newx or newdata
    if(default){
      if(missing(newx)) stop("Input 'newx' must be provided when using 'default' gammi method")
      if(ncol(newx) != (length(object$coefficients) - 1)){
        stop("Input 'newx' must satisfy:  ncol(newx) == (length(object$coef)-1)")
      }
      newnobs <- nrow(newx)
    } else {
      #if(missing(newdata)) stop("Input 'newdata' must be provided when using 'formula' gammi method")
      if(missing(newdata)) newdata <- object$data
      newnobs <- nrow(newdata)
      if(is.null(object$fixed)){
        newx <- spline.model.matrix(object$formula, data = newdata,
                                    Boundary.knots = lapply(object$spline.info$knots, 
                                                            function(x) x[c(1, length(x))]),
                                    knots = object$spline.info$knots,
                                    m = object$spline.info$m,
                                    periodic = object$spline.info$periodic,
                                    xlev = object$spline.info$xlev)[,-1,drop=FALSE]
      } else {
        x0 <- model.matrix(object$fixed, data = newdata, contrasts.arg = object$contrasts)
        x1 <- spline.model.matrix(object$formula, data = newdata,
                                  Boundary.knots = lapply(object$spline.info$knots, 
                                                          function(x) x[c(1, length(x))]),
                                  knots = object$spline.info$knots,
                                  m = object$spline.info$m,
                                  periodic = object$spline.info$periodic,
                                  xlev = object$spline.info$xlev)
        newx <- cbind(x0[,-1], x1[,-1])
      }
    } # end if(default)
    
    ### check type
    type <- as.character(type[1])
    types <- c("link", "response", "terms")
    tyid <- pmatch(type, types)
    if(is.na(tyid)) stop("Input 'type' must be one the following: 'link', 'response', or 'terms'")
    type <- types[tyid]
    
    ### check conf.int
    conf.int <- as.logical(conf.int[1])
    if(!any(conf.int == c(TRUE, FALSE))) stop("Input 'conf.int' must be TRUE or FALSE")
    
    ### check conf.level
    conf.level <- as.numeric(conf.level[1])
    if(conf.level <= 0 | conf.level >= 1) stop("Input 'conf.level' must be between 0 and 1")
    
    ### terms requested?
    if(type == "terms"){
      
      group <- as.integer(as.factor(object$group))
      coefs <- object$coefficients[-1]
      if(se.fit | conf.int) vcovchol <- object$vcovchol[-1,]
      nterms <- length(object$term.labels)
      fit <- matrix(0, newnobs, nterms)
      se <- lwr <- upr <- NULL
      if(se.fit | conf.int) se <- fit
      if(conf.int) {
        lwr <- upr <- fit
        crit <- qnorm(p = 1 - (1 - conf.level)/2)
      } 
      colnames(fit) <- object$term.labels
      if(!is.null(se)) colnames(se) <- object$term.labels
      if(!is.null(lwr)) colnames(lwr) <- object$term.labels
      if(!is.null(upr)) colnames(upr) <- object$term.labels
      for(k in 1:nterms){
        kid <- which(group == k)
        fit[,k] <- newx[,kid,drop=FALSE] %*% coefs[kid]
        if(se.fit | conf.int) se[,k] <- sqrt(rowSums((newx[,kid,drop=FALSE] %*% vcovchol[kid,,drop=FALSE])^2))
        if(conf.int){
          lwr[,k] <- fit[,k] - se[,k] * crit
          upr[,k] <- fit[,k] + se[,k] * crit
        }
      }
      
      res <- list(fit = fit, se = se, lwr = lwr, upr = upr)
      res <- res[!sapply(res, is.null)]
      if(length(res) == 1L) res <- res[[1]]
      return(res)
      
    } # end if(type == "terms")
    
    ### linear predictors
    fit <- as.numeric(object$coefficients[1] + newx %*% object$coefficients[-1])
    pred <- data.frame(fit = fit)
    
    ### standard errors required?
    if(se.fit | conf.int){
      se <- sqrt(rowSums( (cbind(1, newx) %*% object$vcovchol[1:length(object$coefficients),,drop=FALSE])^2 ))
      pred <- cbind(pred, se = se)
    } 
    
    ### confidence interval?
    if(conf.int){
      crit <- qnorm(p = 1 - (1 - conf.level)/2)
      lwr <- fit - se * crit
      upr <- fit + se * crit
      pred <- cbind(pred, lwr = lwr, upr = upr)
    }
    
    ### transform to response scale?
    if(type == "response") {
      pred <- object$family$linkinv(pred)
      if(se.fit | conf.int){
        pred$se <- sqrt( (se^2) * object$family$mu.eta(fit)^2 )
      }
    }
    
    ### return
    if(ncol(pred) == 1L) pred <- unname(unlist(pred))
    return(pred)
    
  } # end predict.gammi