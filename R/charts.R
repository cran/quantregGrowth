charts<-function(fit, k, file=NULL, digits=2, conf.level=0, dataframe=FALSE, 
                 transf=NULL, se.type=c("sandw","boot"), ...){
  if(conf.level<0 || conf.level>=1) stop(" 'conf.level' is wrong.")
  if(length(fit$BB)>1) stop("charts works just with a single smooth variable")
  if(length(k)<=1 && is.integer(k)) {
    m <- min(attr(fit$BB[[1]], "covariate.35"))
    M <- max(attr(fit$BB[[1]], "covariate.35"))
    values<- seq(m, M, length=k)
  } else {
    values<-k
    } 
  d<-data.frame(values)
  names(d) <-attr(fit$BB[[1]], "smoothName")
  if(conf.level>0){
    r<-predict.gcrq(fit, newdata=d, se.fit=TRUE, transf=transf, type=se.type)
    est <- r$fit
    se.est <- r$se.fit
    zalpha<- abs(qnorm((1-conf.level)/2))
    id <-rep(1:ncol(est), each=3)
    r <- matrix(, length(values), ncol(est)*3)
    for(i in 1:ncol(est)){
      r[,(1:ncol(r))[id==i]]<-cbind(est[,i], est[,i]-zalpha*se.est[,i], est[,i]+zalpha*se.est[,i])
    }
    #browser()
    #rownames(r)<- paste(attr(fit$BB[[1]], "smoothName"), round(values,2))
    rownames(r)<- paste(attr(fit$BB[[1]], "smoothName"), round(values,2), sep="=")
    colnames(r) <- as.vector(sapply(colnames(est), function(.x) paste(.x, c("est", "inf", "sup"), sep=" ")))
  } else {
    r<-predict.gcrq(fit, newdata=d, se.fit=FALSE, transf=transf)
    if(!is.matrix(r)) r<- matrix(r, nrow=nrow(d), dimnames=list(NULL, names(r)))
    rownames(r)<- paste(attr(fit$BB[[1]], "smoothName"), round(values,2), sep="=")
  }
  r<-round(r, digits)
  
  if(dataframe && conf.level>0){
    #m <- charts(mm, k=4L, conf.level=.9)
    ntaus<- length(fit$taus)
    
    id.est <- seq(1, ncol(r), by=3)  
    id.inf <- id.est +1   
    id.sup <- id.est +2   
    
    M<-cbind(as.vector(r[,id.est]), as.vector(r[,id.inf]), as.vector(r[,id.sup]))
    out<-data.frame(matrix(, nrow=nrow(M), ncol=1)) 
    for(i in 1:nrow(M)){
      ci <- paste("(", paste(M[i,2:3], collapse=", "), ")", sep ="")
      out[i,1]<-paste(M[i,1], ci, sep=" ")
    }
    
    colnames(out)<-  paste("Est. (", conf.level*100, "% Conf. Interv.)", sep="") 
    r <- cbind(Percentile= paste("tau =", rep(fit$taus, each=length(values))), Values=rep(rownames(r), ntaus), out)
  }
  
  if(is.null(file)){
    return(r)
  } else {
    if(!is.character(file)) stop("If provided, 'file' should be a character")
    write.csv(r, file=file, ...)
  }
  invisible(NULL)
}