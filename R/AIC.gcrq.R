AIC.gcrq<-function (object, ..., k = 2) {
  ll <- logLik.gcrq(object, ...)
  edf <- as.numeric(strsplit(attr(ll, "df"), " ")[[1]]) #as.numeric(attr(ll, "df"))
  if (k < 0) k <- log(attr(ll, "n"))
  r<- -2*ll + k*edf
  #r<-unclass(r)
  r
}
