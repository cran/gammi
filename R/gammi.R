gammi <-
  function(x, ...){
    UseMethod("gammi")
  } # end gammi

print.gammi <-
  function(x, ...){
    cat("\n")
    cat("  Call:  ")
    print(x$call)
    cat("\n")
    cat("Family: ", x$family$family,"\n")
    cat("  Link: ", x$family$link,"\n\n")
    if(!is.na(x$r.squared)) cat("   R^2: ", x$r.squared, "\n")
    cat("logLik: ", x$logLik, "\n")
    cat("   AIC: ", x$aic, "\n")
    cat("   edf: ", sum(x$edf), "\n\n")
  } # end print.gammi