summary.Bsiland<-function(object,...)
{
  #data farame summary
  x<-object
  summaryx=summary(x$result)
  summaryx$call=x$formula
  
  
  cat("Buffer sizes:\n")
  print(round(x$parambuffer,4))
  cat("\n")
  cat("-- Tests are given conditionnaly to the best estimated buffer sizes --" )
  cat("\n")
  print(summaryx)
}

