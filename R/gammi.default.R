gammi.default <-
  function(x,
           y,
           group,
           family = gaussian,
           fixed = NULL,
           random = NULL,
           data = NULL,
           REML = TRUE,
           control = NULL,
           start = NULL,
           verbose = 0L,
           nAGQ = 10L,
           subset, 
           weights, 
           na.action, 
           offset, 
           mustart = NULL, 
           etastart = NULL,
           ...){
    # generalized additive mixed model interface (default method)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # updated: 2024-09-11
    
    
    
    ######***######   INITIALIZATION   ######***######
    
    ### get call
    gammi.call <- match.call()
    
    ### check x and y
    x <- as.matrix(x) + 0.0
    nobs <- nrow(x)
    nvars <- ncol(x)
    ny <- if(is.matrix(y)) nrow(y) else length(y)
    if(ny != nobs) stop("Inputs 'x' and 'y' must satisfy:\nnrow(x) == length(y)  or  nrow(x) == nrow(y)")
    
    ### x and y names
    xnames <- colnames(x)
    if(is.null(xnames)) xnames <- paste0("x", 1:nvars)
    yname <- "y"
    
    ### check group
    if(missing(group)) {
      #group <- ingroup <- 1:nvars
      #gsize <- rep(1L, nvars)
      #names(gsize) <- 1:nvars
      #ngrps <- nvars
      group <- ingroup <- rep(1L, nvars)
      gsize <- nvars
      names(gsize) <- 1L
      ngrps <- 1L
    } else {
      if(length(group) != nvars) stop("Inputs 'x' and 'group' must satisfy:  ncol(x) == length(group)")
      ingroup <- group
      group <- as.factor(group)
      ngrps <- nlevels(group)
      if(ngrps > nvars) stop("Inputs 'x' and 'group' must satisfy:  ncol(x) >= nlevels(group)")
      gsize <- table(group)
      group <- as.integer(group)
    }
    term.labels <- names(gsize)
    
    ### check term.labels
    oterm.labels <- term.labels
    term.labels <- gsub(":", "..x..", term.labels)
    term.labels <- gsub("~", "", term.labels)
    term.labels <- gsub("!", "", term.labels)
    term.labels <- gsub("@", "", term.labels)
    term.labels <- gsub("#", "", term.labels)
    term.labels <- gsub("$", "", term.labels)
    term.labels <- gsub("%", "", term.labels)
    term.labels <- gsub("^", "", term.labels)
    term.labels <- gsub("&", "", term.labels)
    term.labels <- gsub("*", "", term.labels)
    term.labels <- gsub("-", "", term.labels)
    term.labels <- gsub("+", "", term.labels)
    term.labels <- gsub("=", "", term.labels)
    term.labels <- gsub("<", "", term.labels)
    term.labels <- gsub(">", "", term.labels)
    term.labels <- gsub("/", "", term.labels)
    term.labels <- gsub("?", "", term.labels)
    
    ### check family (from glm function)
    if(is.character(family)){
      family <- get(family, mode = "function", envir = parent.frame())
    } else if(is.function(family)){
      family <- family()
    } else if(is.null(family$family)){
      stop("'family' not recognized")
    }
    
    ### check fixed
    nfix <- 0L
    if(!is.null(fixed)) {
      fixed <- as.character(unique(fixed))
      nfix <- length(fixed)
      if(nfix > ngrps) stop("Input 'fixed' contains more terms than 'group'")
      ofixed <- fixed
      fixed <- gsub(":", "..x..", fixed)
      fixed <- gsub("~", "", fixed)
      fixed <- gsub("!", "", fixed)
      fixed <- gsub("@", "", fixed)
      fixed <- gsub("#", "", fixed)
      fixed <- gsub("$", "", fixed)
      fixed <- gsub("%", "", fixed)
      fixed <- gsub("^", "", fixed)
      fixed <- gsub("&", "", fixed)
      fixed <- gsub("*", "", fixed)
      fixed <- gsub("-", "", fixed)
      fixed <- gsub("+", "", fixed)
      fixed <- gsub("=", "", fixed)
      fixed <- gsub("<", "", fixed)
      fixed <- gsub(">", "", fixed)
      fixed <- gsub("/", "", fixed)
      fixed <- gsub("?", "", fixed)
      for(j in 1:nfix){
        if(!(ofixed[j] %in% oterm.labels)) stop("Input 'fixed' contains ", fixed[j], " which is not in 'group'")
      }
    }
    
    ### check random
    if(!is.null(random)) random <- as.formula(random)
    
    ### check REML
    REML <- as.logical(REML[1])
    if(!any(REML == c(TRUE, FALSE))) stop("Input 'REML' must be TRUE or FALSE")
    
    ### check control
    if(is.null(control)){
      if(family$family == "gaussian"){
        control <- lme4::lmerControl()
      } else {
        control <- lme4::glmerControl()
      }
    } 
    
    ### check verbose
    verbose <- as.integer(verbose[1])
    if(verbose < 0L) stop("Input verbose should be a non-negative integer")
    
    ### check nAGQ
    nAGQ <- as.integer(nAGQ[1])
    if(nAGQ < 1L) stop("Input nAGQ should be a positive integer")
    
    ### check weights
    if(missing(weights)) {
      weights <- rep(1.0, nobs) 
    } else {
      weights <- as.numeric(weights)
      if(length(weights) != nobs) stop("Input 'weights' must have the same length as reponse.")
      if(any(weights <= 0)) stop("Input 'weights' must be positive")
    }
    
    ### check offset
    if(missing(offset)) {
      offset <- rep(0.0, nobs) 
    } else {
      offset <- as.numeric(offset)
      if(length(offset) != nobs) stop("Input 'offset' must have the same length as reponse.")
    }
    
    ### check subset
    if(missing(subset)) {
      subset <- rep(TRUE, nobs) 
    } else {
      subset <- as.numeric(subset)
      if(length(subset) != nobs) stop("Input 'subset' must have the same length as reponse.")
    }
    
    ### group z-score x
    xorig <- x
    xmean <- colMeans(x)
    x <- x - matrix(xmean, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
    xscale <- rep(NA, ngrps)
    names(xscale) <- term.labels
    tvec <- rep(NA, nvars)
    for(k in 1:ngrps){
      id <- which(group == k)
      xscale[k] <- sqrt(mean(x[,id]^2))
      x[,id] <- x[,id] / xscale[k]
      tvec[id] <- xscale[k]
    }

    ### reverse transformation
    tmat <- diag(c(1, 1/tvec))
    tmat[1,] <- c(1, -xmean/tvec)
    
    
    
    ######***######   LMER FORMULA AND DATA   ######***######
    
    ### make lmer formula
    fixid <- NULL
    lmer.form <- paste0("(1 | ", term.labels, ")")
    if(!is.null(fixed)){
      fixid <- pmatch(fixed, term.labels)
      lmer.form[fixid] <- fixed
    } 
    lmer.form <- paste(lmer.form, collapse = " + ")
    
    ### any random effects?
    if(!is.null(random)){
      random.char <- as.character(random)
      if(length(random.char) != 2L | random.char[1] != "~"){
        stop("Input 'random' must be a one-sided formula, e.g.,\n ~ (1+x|subject) + (0+z|stimulus)")
      }
      lmer.form <- paste(lmer.form, random.char[2], sep = " + ")
    }
    
    ### add response
    lmer.form <- as.formula(paste0("y ~ ", lmer.form))
    
    ### make data
    lmer.data <- vector("list", ngrps + 1L)
    names(lmer.data) <- c(term.labels, "y")
    for(k in 1:ngrps){
      if(k %in% fixid){
        lmer.data[[k]] <- x[,(group == k)]
      } else {
        lmer.data[[k]] <- factor(rep(1:gsize[k], length.out = nobs))
      }
    }
    lmer.data[[ngrps+1]] <- y
    
    ### any data?
    if(!is.null(data)) lmer.data <- c(lmer.data, as.list(data))
    
    
    ######***######   LMER MODEL FITTING   ######***######
    
    ### process formula to build model matrices
    if(family$family == "gaussian"){
      lform <- lme4::lFormula(lmer.form, 
                              data = lmer.data, 
                              REML = REML,
                              subset = subset,
                              weights = weights,
                              na.action = na.action,
                              offset = offset,
                              control = control)
    } else {
      lform <- lme4::glFormula(lmer.form, 
                               data = lmer.data,
                               family = family,
                               subset = subset,
                               weights = weights,
                               na.action = na.action,
                               offset = offset,
                               start = start,
                               mustart = mustart,
                               control = control)
    } # end if(family$family == "gaussian")
    
    ### replace dummy coded variables
    lterm <- names(lform$reTrms$cnms)
    for(k in 1:ngrps){
      if(k %in% fixid) next
      id <- which(lterm == term.labels[k])
      xindex <- which(group == k)
      zindex <- (lform$reTrms$Gp[id]+1):lform$reTrms$Gp[id+1]
      lform$reTrms$Zt[zindex,] <- as(t(x[,xindex]),"dgCMatrix")
    } # end for(i in 1:length(lterm))
    
    ### make and optimize deviance function
    if(family$family == "gaussian"){
      devfun <- do.call(lme4::mkLmerDevfun, lform)
      opt <- lme4::optimizeLmer(devfun, control = control$optCtrl)
    } else {
      devfun <- do.call(lme4::mkGlmerDevfun, lform)
      opt <- lme4::optimizeGlmer(devfun, control = control$optCtrl)
    }
    
    
    
    ######***######   COLLECT RESULTS   ######***######
    
    ### package merMod results
    mer <- lme4::mkMerMod(environment(devfun), opt, lform$reTrms, fr = lform$fr)
    
    ### extract blups
    random.coefs <- lme4::ranef(mer)
    
    ### extract fixed effects
    fixed.coefs <- lme4::fixef(mer)
    fixed.names <- names(fixed.coefs)
    
    ### put coefficients in order
    coefs <- rep(0.0, nvars)
    names(coefs) <- xnames
    fxnames <- NULL
    for(k in 1:ngrps){
      xid <- which(group == k)
      if(k %in% fixid){
        fxnames <- c(fxnames, colnames(x[,xid,drop=FALSE]))
        if(gsize[k] == 1L){
          zid <- which(fixed.names == term.labels[k])
        } else {
          zid <- pmatch(paste0(term.labels[k], xnames[xid]), fixed.names)
        }
        coefs[xid] <- fixed.coefs[zid]
      } else {
        zid <- which(lterm == term.labels[k])
        coefs[xid] <- unlist(random.coefs[[zid]])
      }
    } # end for(k in 1:ngrps)
    if(fixed.names[1] == "(Intercept)") {
      coefs <- c(fixed.coefs[1], coefs)
      fxnames <- c("(Intercept)", fxnames)
    }
    
    ### get blups (if applicable)
    blups <- NULL
    if(!is.null(random)){
      blups <- random.coefs[which(!(lterm %in% term.labels))]
    }
    
    ### linear predictors
    lp <- as.numeric(coefs[1] + x %*% coefs[-1])
    mu <- family$linkinv(lp)
    
    
    ######***######   COVARIANCE MATRIX   ######***######
    
    ### extract blocks of cross-product matrix
    L <- mer@pp$L()
    if(inherits(L, "CHMsimpl")){
      L <- as(mer@pp$L(), "dtCMatrix")
    } else if(inherits(L, "CHMsuper")){
      L <- as(mer@pp$L(), "dgCMatrix")
    } 
    RZX <- mer@pp$RZX
    RX <- mer@pp$RX()
    
    ### evaluate Linv and RXinv
    Linv <- as.matrix(Matrix::solve(L))
    RXinv <- base::solve(RX)
    
    ### create cholesky factor of vcov matrix
    zeros <- matrix(0, nrow(RXinv), nrow(Linv))
    vcovchol <- rbind(crossprod(mer@pp$Lambdat, cbind(t(Linv), -t(Linv) %*% RZX %*% RXinv)), 
                      cbind(zeros, RXinv))
    
    ### assign rownames to vcovchol
    lmer.names <- NULL
    for(j in 1:length(lterm)){
      k <- which(term.labels == lterm[j])
      if(length(k) >= 1L){
        # term entered through x
        lmer.names <- c(lmer.names, xnames[group == k])
      } else {
        # term entered through random
        cnames <- colnames(random.coefs[[j]])
        for(i in 1:ncol(random.coefs[[j]])){
          lmer.names <- c(lmer.names, paste0(rownames(random.coefs[[j]]), ".", cnames[i]))
        }
      } # end if(length(k) >= 1L)
    } # end for(j in 1:length(lterm))
    lmer.names <- c(lmer.names, fxnames)
    rownames(vcovchol) <- lmer.names
    
    ### reorder rows to match ordering of coefs
    oid <- pmatch(lmer.names, names(coefs))
    naid <- is.na(oid)
    if(any(naid)) oid[naid] <- max(oid, na.rm = TRUE) + 1
    neword <- sort(oid, index = TRUE)$ix
    vcovchol <- as.matrix(vcovchol[neword,])
    
    
    
    ######***######   RETURN RESULTS   ######***######
    
    ### retransform coefs and vcovchol
    coefs <- as.numeric(tmat %*% coefs)
    names(coefs) <- c("(Intercept)", xnames)
    indx <- 1:length(coefs)
    vcovchol[indx,] <- tmat %*% vcovchol[indx,,drop=FALSE]
    
    ### extract and retransform variance components
    temp <- as.data.frame(lme4::VarCorr(mer))
    for(k in 1:ngrps){
      id <- which(temp$grp == term.labels[k])
      if(length(id) > 0){
        temp$vcov[id] <- temp$vcov[id] / xscale[k]^2 
        temp$sdcor[id] <- sqrt(temp$vcov[id])
      }
    }
    
    ### calculate leverages and df
    lev <- rowSums((cbind(1, xorig) %*% vcovchol[1:length(coefs),,drop=FALSE])^2)
    if(family$family == "gaussian" && family$link == "identity"){
      edf <- sum(lev * weights)
    } else {
      edf <- sum(lev * weights * family$mu.eta(lp)^2 / family$variance(mu))
    }
    df.random <- 0.0
    if(!is.null(random)){
      df.random <- length(mer@theta) - (length(term.labels) - nfix)
    }
    
    ### calculate fit information
    aic.scale.penalty <- ifelse(family$family %in% c("gaussian", "Gamma", "inverse.gaussian"), 1, 0)
    if(is.null(random)){
      dev <- sum(family$dev.resids(y, mu, weights))
      wtdmu <- sum(y * weights) / sum(weights)
      null.dev <- sum(family$dev.resids(y, wtdmu, weights))
      LL <- family$aic(y, nobs, mu, weights, dev)
      LL <- (-1/2)*(LL - 2*aic.scale.penalty)
    } else {
      dev <- mer@devcomp$cmp[[ifelse(family$family == "gaussian" && REML, "REML", "dev")]]
      null.dev <- NA
      LL <- -dev/2      # -2*logLik = devCrit
    } # end if(is.null(random))
    
    ### change ".." back to the ":" 
    term.labels <- gsub("..x..", ":", term.labels)
    
    ### dispersion parameter
    if(family$family %in% c("gaussian", "Gamma", "inverse.gaussian")){
      if(is.null(random)){
        disp <- dev / (nobs - edf)
      } else {
        disp <- ifelse(REML, mer@devcomp$cm['sigmaREML']^2, mer@devcomp$cm['sigmaML']^2)
      }
    } else {
      disp <- 1.0
    }
    
    ### collect results
    res <- NULL
    res$fitted.values <- family$linkinv(lp)
    res$linear.predictors <- lp
    res$coefficients <- coefs
    res$random.coefficients <- blups
    res$term.labels <- term.labels
    res$dispersion <- disp
    res$vcovchol <- sqrt(res$dispersion) * vcovchol
    res$family <- family
    res$logLik <- LL
    res$aic <- 2 * (edf + df.random - LL)
    res$deviance <- dev
    res$null.deviance <- ifelse(is.null(random), null.dev, NA)
    res$r.squared <- ifelse(is.null(random), 1 - dev / null.dev, NA)
    res$nobs <- nobs
    res$leverages <- lev
    res$edf <- edf
    res$df.random <- df.random
    res$df.residual <- nobs - edf
    res$x <- xorig
    res$group <- ingroup
    res$scale <- xscale
    res$fixed <- fixed
    res$random <- random
    res$mer <- mer
    res$VarCorr <- temp
    res$call <- gammi.call
    
    ### return
    class(res) <- "gammi"
    return(res)
    
  } # end gammi.default