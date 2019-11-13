AIC.Fsiland<-function(object,...,k=2)
{
  cat(paste("AIC = ",signif(object$AIC,digits=6)))
  invisible(object$AIC)
}
