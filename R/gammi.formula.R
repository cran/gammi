gammi.formula <-
  function(formula, 
           data, 
           family = gaussian,
           fixed = NULL, 
           random = NULL, 
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
    # generalized additive mixed model interface (formula method)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # updated: 2024-07-16
    
    
    ### get call
    gammi.call <- match.call()
    
    ### check formula
    formula <- as.formula(formula)
    formula0 <- stats::update(formula, . ~ . - 1)
    
    ### any data?
    if(missing(data)) data <- model.frame(formula)
    
    ### extract response
    allvars <- all.vars(formula)
    if(length(allvars) < 2L) stop("Input 'formula' must have 2 or more variables (response and 1+ predictor)")
    yname <- allvars[1]
    if(any(grepl("..x..", allvars[-1], fixed = TRUE))) stop("Input 'formula' must not contain variables with '..x..' appearing in their name.")
    y <- data[[yname]]
    
    ### y included in data
    if(is.null(y)){
      tempdata <- model.frame(formula, data = data)
      y <- tempdata[[yname]]
      dnames <- names(data)
      data <- cbind(y, data)
      names(data) <- c(yname, dnames)
    } 
    
    ### all variables in data?
    xvars <- c(all.vars(fixed), all.vars(formula)[-1])
    check <- match(xvars, names(data))
    if(any(is.na(check))) stop("Input 'data' is missing variables referenced in 'fixed' and/or 'formula'.")
    
    ### create x and group
    ssi <- NULL
    if(is.null(fixed)){
      contr.list <- vector()
      x <- spline.model.matrix(formula0, data = data, ...)
      group <- factor(attr(x, "term.labels")[attr(x, "assign")], levels = attr(x, "term.labels"))
      terms.fixed <- NULL
      ssi <- list(knots = attr(x, "knots"),
                  m = attr(x, "m"),
                  periodic = attr(x, "periodic"),
                  xlev = attr(x, "xlev"))
    } else {
      fixvars <- all.vars(fixed)
      contr.list <- vector("list", length(fixvars))
      names(contr.list) <- fixvars
      for(j in 1:length(fixvars)){
        jc <- class(data[,fixvars[j]])[1]
        if(jc == "factor" | jc == "character"){
          contr.list[[j]] <- "contr.sum"
          data[,fixvars[j]] <- as.factor(data[,fixvars[j]])
          contrasts(data[,fixvars[j]]) <- contr.sum(nlevels(data[,fixvars[j]]))
        }
      }
      contr.list <- contr.list[!sapply(contr.list, is.null)]
      terms.fixed <- attr(terms(fixed), "term.labels")
      x0 <- model.matrix(fixed, data = data)
      g0 <- attr(x0, "assign")
      if(g0[1] == 0L){
        g0 <- g0[-1]
        x0 <- x0[,-1]
      }
      x1 <- spline.model.matrix(formula0, data = data, ...)
      g1 <- attr(x1, "assign")
      x <- cbind(x0, x1)
      group <- c(terms.fixed[g0], attr(x1, "term.labels")[g1])
      group <- factor(group, levels = c(terms.fixed, attr(x1, "term.labels")))
      ssi <- list(knots = attr(x1, "knots"),
                  m = attr(x1, "m"),
                  periodic = attr(x1, "periodic"),
                  xlev = attr(x1, "xlev"))
    }
    
    ### call gammi.default
    res <- gammi.default(x = x, y = y, group = group, family = family,
                         fixed = terms.fixed, random = random, data = data,
                         REML = REML, control = control, start = start, 
                         verbose = verbose, nAGQ = nAGQ, subset = subset, 
                         weights = weights, na.action = na.action, 
                         offset = offset, mustart = mustart, etastart = etastart)
    
    ### replace x with data
    res$x <- NULL
    res$data <- data
    
    ### update call and fixed
    res$call <- gammi.call
    res$fixed <- fixed
    
    ### add spline info and formula
    res$contrasts <- if(length(contr.list) > 0L) contr.list else NULL
    res$spline.info <- ssi
    res$formula <- formula
    
    ### return
    return(res)
    
  } # end gammi.formula