AIC.gcrq<-function (object, ..., k = 2, bondell=FALSE) {
  if(bondell){
    ll<-logLik.gcrq(object, summ=TRUE)
    n <- attributes(ll)$n
    df <- as.numeric(attributes(ll)$df)
    RhoTot<- sum(object$rho)
    r<- log(RhoTot/n) + log(n)*df/(2*n)
    # attr(r,"df")<-paste(df)
    # attr(r,"n") <- n
    # attr(r,"class") <- "logLik"
  } else {
    ll <- logLik.gcrq(object, ...)
    edf <- as.numeric(strsplit(attr(ll, "df"), " ")[[1]]) #as.numeric(attr(ll, "df"))
    if (k < 0) k <- log(attr(ll, "n"))
    r<- -2*ll + k*edf
    attr(r,"df")<-NULL
    attr(r,"n") <- NULL
    attr(r,"class") <- NULL
    
  }
    r
}
