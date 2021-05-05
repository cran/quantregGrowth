logLik.gcrq<- function (object, summ=TRUE, ...) {
  taus<- object$taus
  val<- object$rho
  n <- if(length(taus)>1) nrow(object$residuals) else length(object$residuals) 
  edf <- colSums(object$edf.j)
  val <- n * (log(taus * (1 - taus)) - 1 - log(object$rho/n))
  if(summ) {
    edf<-sum(edf)
    val<-sum(val)
  }
  attr(val, "df")<- paste(round(edf,3), collapse = " ") #round(edf,3) #
  attr(val, "n") <-  n   
  class(val)<-"logLik"
  val
}




