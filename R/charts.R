charts<-function(fit, k, file=NULL, digits=2, ...){
  if(length(fit$BB)>1) stop("charts works just with a single smooth variable")
  if(length(k)<=1) {
    m <- min(attr(fit$BB[[1]], "covariate.35"))
    M <- max(attr(fit$BB[[1]], "covariate.35"))
    values<- seq(m, M, length=k)
  } else {
    values<-k
    } 
  d<-data.frame(values)
  names(d) <-attr(fit$BB[[1]], "smoothName")
  r<-predict.gcrq(fit, newdata=d, se.fit=FALSE)
  r<-round(r, digits)
  rownames(r)<- paste(attr(fit$BB[[1]], "smoothName"), round(values,2))
  if(is.null(file)){
    return(r)
  } else {
    if(!is.character(file)) stop("If provided, 'file' should be a character")
    write.csv(r, file=file, ...)
  }
  invisible(NULL)
}