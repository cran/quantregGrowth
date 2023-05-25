vcov.gcrq <- function(object, term, type=c("sandw","boot"),...){
    type<-match.arg(type)
    if(type=="boot"){
        if(!("boot.coef" %in% names(object))) {
              #object <- update(object, n.boot=100)
              #warning("The fit misses boot replicates.. 100 runs have been performed now")
              stop("the object fit misses the boot replicates")  
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
            stop(" 'term' is not allowed. Only the full covariance matrix may be returned.")
            if(!term %in% all.vars(formula(object))[-1]) stop(" 'term' is not a variable in the model")
            f.smooth<-function(.x){
                #id<-startsWith(colnames(.x),paste(term, "ps.",sep="."))# non funziona con vc i cui nomi sono.. "x1:x2.ps.1"  "x1:x2.ps.2"
                id<- grep(paste(term, "ps.",sep="."),colnames(VCOV[[1]]))
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
        #stop("only 'type=boot' allowed ")
        #.. uses a kernel estimate of the sandwich as proposed by Powell(1990).
        taus<-object$taus
        n.tau<-length(taus) #if(length(tau)>1) stop("multiple taus not yet allowed")
        n<-nrow(as.matrix(object$fitted.values))
        ##################
        x.tildeUp<-object$x[1:n,,drop=FALSE] 
        
        #browser()
        
        if(!is.null(object$lambda)){
            n.coefs<- sapply(object$BB, function(.x)length(attr(.x,"coef.names"))) #n. di coef per ogni termine smooth
            n.coefs<-n.coefs[grep("ps\\(", names(n.coefs))]
            lambda<- matrix(as.matrix(object$lambda), ncol=n.tau, nrow=length(n.coefs)) 
            Dt <- t(object$D.matrix)
            nrowD<- ncol(Dt)
        }
        
        p <- ncol(x.tildeUp)
        VCOV<-NULL
        for(j in 1:n.tau) {
            tau<- taus[j]
            uhat <- as.matrix(object$residuals)[,j]
            h <- bandwidth.rq(tau, n)
            h <- (qnorm(tau + h) - qnorm(tau - h)) * min(sqrt(var(uhat)), 
                                                         (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34)
            f <- dnorm(uhat/h)/h
            if(!is.null(object$lambda)){
                D.lambda <- t(c(rep(0,object$pLin), rep(lambda[,j], n.coefs))*Dt)
                x.tilde<-rbind(x.tildeUp, D.lambda)
                f<-c(f, rep(1, nrowD))
            } else {
                x.tilde<-rbind(x.tildeUp)
            }

            A<- solve(crossprod(sqrt(f)* x.tilde)) #aggiungere + 1e-08? #solve(t(x.tilde)%*% diag(f) %*% x.tilde)

            # U<-x.tilde[1:n,]*(tau-1*(uhat<0))
            # RIS<-matrix(NA, 2000, p)
            # for(i in 1:2000){
            #     id<-sample(n, replace=TRUE)
            #     RIS[i,]<-  colSums(U[id,])
            # }
            # v<- var(RIS) 
            # m <- A %*% v %*% A
            m <- tau * (1 - tau) * A %*% crossprod(x.tilde) %*% A #x.tilde[1:n,]
            colnames(m)<-rownames(m)<-rownames(as.matrix(object$coefficients))
            VCOV[[length(VCOV)+1]]<-m
        }
        names(VCOV)<-paste(object$tau) #colnames(object$coefficients)
        return(VCOV) #se<- sqrt(rowSums((X %*% v) * X)) #sqrt(diag(X%*%v%*%t(X)))
        #stesso risultato.. E' piu veloce?
        #fxxinv <- diag(p)
        #fxxinv <- backsolve(qr(sqrt(f) * x.tilde)$qr[1:p, 1:p, drop = FALSE], fxxinv)
        #fxxinv <- fxxinv %*% t(fxxinv)
        #VCOV <- tau * (1 - tau) * fxxinv %*% crossprod(x.tilde) %*% fxxinv 
        }
}
        
        
        