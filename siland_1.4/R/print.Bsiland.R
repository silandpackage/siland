print.Bsiland<-function(x,...){
  cat("Model: ")
  print(x$formula,showEnv=F)
  cat("\n")
  
  cat("Landscape variables: ")
  cat(colnames(x$buffer))
  cat("\n \n")
  
  cat("Coefficients:\n")
  print(round(x$coefficients,3))
  if(x$modelType=="LMM" ||x$modelType=="GLMM")
  {
    cat("\nRandom effects:\n")
    if(x$modelType=="LMM"||x$modelType=="GLMM")
      print(x$rand.StdDev)
  }
  cat("\n")
  cat("standard error: ")
  cat(x$sd.error)
  cat("\n")
  cat(paste("AIC: ",round(x$AIC,2), "\t", "AIC (no landscape): ",round(x$AIC0,2),sep=""))
  cat("\n")
  cat("(No landscape effect) p-value: ")
  if(x$pval0 < 1e-16) cat("<1e-16") else cat(x$pval0)
  
  
  
}
