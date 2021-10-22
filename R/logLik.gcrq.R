logLik.gcrq<- function (object, summ=FALSE, ...) {
  taus<- object$taus
  val<- object$rho
  n.tau<-length(taus)
  n <- if(length(taus)>1) nrow(object$residuals) else length(object$residuals) 
  edf <- colSums(object$edf.j)
  if(summ) {
    val <- n * (sum(log(taus * (1 - taus))) - n.tau - log(sum(object$rho)/(n*n.tau)) )
    edf<-sum(edf)
    } else {
      val <- n * (log(taus * (1 - taus)) - 1 - log(object$rho/n))
  }

  #sic<-log(sum(x$rho)/(n*n.tau))+sum(x$edf.j)*log(n*n.tau)/(2*n*n.tau)
  
  attr(val, "df")<- paste(round(edf,3), collapse = " ") #round(edf,3) #
  attr(val, "n") <-  n   
  class(val)<-"logLik"
  val
}




