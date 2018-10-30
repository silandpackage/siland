BIC.siland<-function(object,...)
{
  N=length(object$fitted)
  loglik=object$loglik
  nparam=length(object$coefficients)
  
  if(object$family!="gaussian")
  BIC  =-2*loglik +nparam*log(N)
  else
  BIC=-2*loglik +(nparam+1)*log(N)
  
  cat(paste("BIC = ",signif(BIC,digits=6)))
  invisible(BIC)
}