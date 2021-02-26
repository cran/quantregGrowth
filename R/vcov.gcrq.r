vcov.gcrq <- function(object, term, type=c("boot","sandw"),...){
    type<-match.arg(type)
    if(type=="boot"){
        if(!("boot.coef" %in% names(object))) {
              object <- update(object, n.boot=100)
              warning("The fit misses boot replicates.. 100 runs have been performed now")
              }
        #if(is.null(object$boot.coef)) stop(" 'vcov.gcrq' works only with boot")
        coef.boot<-object$boot.coef
        n.boot<-dim(coef.boot)[3]
        n.tau<- dim(coef.boot)[2]
        n.coef<-dim(coef.boot)[1]
        VCOV<-NULL
        for(j in 1:n.tau) {
            m<-var(t(coef.boot[,j,]))
            colnames(m)<-rownames(m)<-rownames(as.matrix(object$coefficients))
            VCOV[[length(VCOV)+1]]<-m
            }
        if(!missing(term)){
            if(!term %in% all.vars(formula(object))[-1]) stop(" 'term' is not a variable in the model")
            f.smooth<-function(.x){
                id<-startsWith(colnames(.x),paste(term, "ps.",sep="."))
                r<-.x[id,id,drop=TRUE]
                r
                }
            f.lin<-function(.x){
                id<-grep(term,colnames(.x))
                r<-.x[id,id,drop=TRUE]
                r
                }
            f.ok<- if(term %in% names(object$BB)) f.smooth else f.lin

            VCOV <- lapply(VCOV, f.ok)
            }
        names(VCOV)<-paste(object$tau) #colnames(object$coefficients)
        return(VCOV)
    } else {
        stop("only 'type=boot' allowed ")
        #.. uses a kernel estimate of the sandwich as proposed by Powell(1990).
        tau<-object$taus
        if(length(tau)>1) stop("multiple taus not yet allowed")
        uhat <- object$residuals
        n<-length(uhat)
        x.tilde<-object$x[1:n,]
        #x<- x.tilde[1:n,]
        p<-ncol(x.tilde)
        h <- bandwidth.rq(tau, n)
        h <- (qnorm(tau + h) - qnorm(tau - h)) * min(sqrt(var(uhat)), 
                                                  (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34)
        f <- dnorm(uhat/h)/h
        f<-c(f, rep(1, nrow(x.tilde)-n))
        
        A<- solve(t(x.tilde)%*% diag(f) %*% x.tilde)
        cov <- tau * (1 - tau) * A %*% crossprod(x.tilde) %*% A
        return(cov) #se<- sqrt(rowSums((X %*% v) * X)) #sqrt(diag(X%*%v%*%t(X)))
        
        fxxinv <- diag(p)
        fxxinv <- backsolve(qr(sqrt(f) * x.tilde)$qr[1:p, 1:p, drop = FALSE], fxxinv)
        fxxinv <- fxxinv %*% t(fxxinv)
        VCOV <- tau * (1 - tau) * fxxinv %*% crossprod(x.tilde) %*% fxxinv 
        return(VCOV)
        }
}
        
        
        