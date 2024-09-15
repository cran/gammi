# Some code (for startup message) is re-purposed from the R package
# "mclust" (Scrucca et al., 2016) https://cran.r-project.org/package=mclust

gammiStartupMessage <- 
  function(){
    msg <- c(paste0("                                    
      __ _  __ _ _ __ ___  _ __ ___ (_)
     / _` |/ _` | '_ ` _ \\| '_ ` _ \\| |
    | (_| | (_| | | | | | | | | | | | |
     \\__, |\\__,_|_| |_| |_|_| |_| |_|_|
     |___/                 version ",       
                    packageVersion("gammi"), "\n"),
             "\nType 'citation(\"gammi\")' to cite this package.\n")
    return(msg)
  }

.onAttach <- 
  function(lib, pkg){
    msg <- gammiStartupMessage()
    if(!interactive()) msg[1] <- paste("Package 'gammi' version", packageVersion("gammi"))
    packageStartupMessage(msg)      
    invisible()
  }