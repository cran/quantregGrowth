ncross.rq.fitXB <-
  function(y, x, B=NULL, X=NULL, taus, monotone=FALSE, concave=FALSE, #tolto argomento interc=FALSE,
           nomiBy=NULL, byVariabili=NULL,
           ndx=10, deg=3, dif=3, lambda=0, eps=.0001, var.pen=NULL, penMatrix=NULL,
           lambda.ridge=0,  dropcList=FALSE, decomList=FALSE, vcList=FALSE, dropvcList=FALSE, centerList=FALSE, 
           ridgeList=FALSE, ps.matrix.list=FALSE, colmeansB=NULL, Bconstr=NULL, adjX.constr=TRUE, 
           adList=FALSE, it.j=10, myeps=NULL, ...){
    #usare spam??

    plott=0 #prima era un argomento
    adj.middle=FALSE
    err.rho=TRUE #prima era un argomento
    u.rho=TRUE ##prima era un argomento
    if(is.null(myeps)) myeps<-1e-06
    #--------------------------------------------------------
    Rho <- function(u, tau) u * (tau - (u < 0))
    #--------------------------------------------------------
    bspline <- function(x, ndx, xlr = NULL, knots=NULL, deg = 3, deriv = 0, outer.ok=FALSE) {
      # x: vettore di dati
      # xlr: il vettore di c(xl,xr)
      # ndx: n.intervalli in cui dividere il range
      # deg: il grado della spline
      # Restituisce ndx+deg basis functions per ndx-1 inner nodi
      #ci sono "ndx+1" nodi interni + "2*deg" nodi esterni
      #    require(splines)
      if(is.null(knots)) {
        if (is.null(xlr)) {
          xl <- min(x) - 0.01 * diff(range(x))
          xr <- max(x) + 0.01 * diff(range(x))
        }
        else {
          if (length(xlr) != 2)
            stop("quando fornito, xlr deve avere due componenti")
          xl <- xlr[1]
          xr <- xlr[2]
        }
        dx <- (xr - xl)/ndx
        knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
      }
      #else {
      #if(length(knots)!=(ndx+1+2*deg)) stop("errore nel numero di nodi fornito")
      #}
      B <- splineDesign(knots, x, ord = deg + 1, derivs = rep(deriv, length(x)), outer.ok=outer.ok)
      B
    }
    #-------------------------------------
    blockdiag <- function(...) {
      args <- list(...)
      nc <- sapply(args,ncol)
      cumnc <- cumsum(nc)
      ##  nr <- sapply(args,nrow)
      ## NR <- sum(nr)
      NC <- sum(nc)
      rowfun <- function(m,zbefore,zafter) {
        cbind(matrix(0,ncol=zbefore,nrow=nrow(m)),m,
              matrix(0,ncol=zafter,nrow=nrow(m)))
      }
      ret <- rowfun(args[[1]],0,NC-ncol(args[[1]]))
      for (i in 2:length(args)) {
        ret <- rbind(ret,rowfun(args[[i]],cumnc[i-1],NC-cumnc[i]))
      }
      ret
    }
    #-------------------------------------
    build.D<-function(var.pen.ok, p.ok, dif.ok, lambda.ok, penMatrix.ok=NULL, dropc=FALSE, decomp=FALSE, dropvc=FALSE, ridge=FALSE){
      #Build the penalty matrix to be appended below the design matrix
      #var.pen.ok: a string to mean the varying penalty
      #p.ok: no. of coefficients
      #dif.ok: difference order (if dif.ok<=0 the identity matrix is returned)
      #lambda.ok: spar parameter
      #dropvc: TRUE if the penalty refers to VC models wherein the basis has been buit with identif. constr. (dropping the 1 column basis)
      if(!is.null(penMatrix.ok)){
        D.ok<- penMatrix.ok
        return(D.ok) 
      } else {
        if(dropvc) p.ok<-p.ok-1
        if(is.null(var.pen.ok)){
          xx.var.pen <- rep(1,(p.ok-dif.ok))
        } else {
          f.var.pen <- function(k) eval(parse(text = var.pen.ok))
          xx.var.pen <- 1:(p.ok - dif.ok)
          xx.var.pen <- f.var.pen(max(xx.var.pen))
        }
        
        #D.ok<-if(decomp || ridge) xx.var.pen*diag(p.ok) else xx.var.pen*diff(diag(p.ok), diff=dif.ok)
        D.ok<-if(dif.ok<=0) xx.var.pen*diag(p.ok) else xx.var.pen*diff(diag(p.ok), diff=dif.ok)
      }
      if(dif.ok>=0 && dropc) D.ok<-D.ok[,-1] # D%*%C infatti deve risultare t(C)%*%Penalty%*%C
      if(dropvc) D.ok<- cbind(0,D.ok) #rbind(0,cbind(0,D.ok))
      D.ok<-lambda.ok*D.ok
      D.ok
    }
    #-------------------------------------
    edf.rqXB <-function(obj, tau, id.coef, vMonot, vConc, pesiL1, output, 
                        df.nc=TRUE, DFvalues=NULL, D.matrix, Bconstr, D1, D2, Y, X){
      #versione per il calcolo edf in ..fitXB()
      #restituisce una lista dove la prima componente e' df.all, la seconda e' la somma per gruppi...
      #obj$DF.NEG e obj$DF.POS.. per i df accounting for noncross..
      #output: which df should be returned?
      #   4 = come 3, ma tutto (anche pen) attraverso su approx quadratica basata sui valori assoluti
      #   3 = the trace of the pseudo-hat matrix based on the smooth approximation (for REML-like estimate)
      #   2 = the number of non-null (penalized) coefficients (for REML-like estimate)
      #   1 = the number of all (penalized) coefficients (for ML estimate)
      # dalla versione 1.0-0 solo le opzioni 2 o 4 sono contemplate.. 
      #Restituisce i df attraverso diag(HatMatrix). tau deve essere uno scalare 
      #Se return.all=TRUE restituisce diag(HatMatrix), cioe' un vettore p-dimensionale. (p=n.variabili)
      #   se return.all=FALSE restituisce un vettore pari al numero dei "gruppi" di variabili
      #   ogni gruppo e' o l'insieme dei termini parametrici o un termine smooth.
      #g: percentuale per calcolare il valore soglia (solo se output=3)
      #id.coef: identif di appartenenza di coeff deve contenere un attributo "nomi (cio? attr(id.coef,"nomi")
      make.positive.definite<-function (m, tol) {
        #preso da  package ?lqmm? 1.5
        #Se M=XtWX+P non ? invertibile calcolare solve(make.positive.definite(M))? da provare
        # vedi comunque qualche altra funzione nel package corpcor
        if (!is.matrix(m)) m = as.matrix(m)
        d = dim(m)[1]
        if (dim(m)[2] != d) stop("Input matrix is not square!")
        es = eigen(m)
        esv = es$values
        if (missing(tol)) tol = d * max(abs(esv)) * .Machine$double.eps
        delta = 2 * tol
        tau1 = pmax(0, delta - esv)
        dm = es$vectors %*% diag(tau1, d) %*% t(es$vectors) #tau1 ? gi? un vettore, diag(tau1) dovrebbe essere sufficiente..
        #dm = crossprod(tau1*es$vectors) 
        return(m + dm)
      }
      #-----------
      blockdiag <- function(...) {
        args <- list(...)
        nc <- sapply(args,ncol)
        cumnc <- cumsum(nc)
        ##  nr <- sapply(args,nrow)
        ## NR <- sum(nr)
        NC <- sum(nc)
        rowfun <- function(m,zbefore,zafter) {
          cbind(matrix(0,ncol=zbefore,nrow=nrow(m)),m,
                matrix(0,ncol=zafter,nrow=nrow(m)))
        }
        ret <- rowfun(args[[1]],0,NC-ncol(args[[1]]))
        for (i in 2:length(args)) {
          ret <- rbind(ret,rowfun(args[[i]],cumnc[i-1],NC-cumnc[i]))
        }
        ret
      }
      #--------------
      n<- obj$n
      interc<-"(Intercept)"%in% rownames(as.matrix(obj$coef))
      nomiGruppi<-attr(id.coef,"nomi") #"Xlin"(eventuale) + nomi smooth
      b<-if(is.matrix(obj$coef)) obj$coef[,paste(tau)] else obj$coef
      n.par<-tapply(id.coef, id.coef, length) #num. dei parametri di ciascun "gruppo": lineari, 1st smooth, 2nd smooth,..
      n.Gruppi<-length(n.par)#length(n.par) e' il numero dei gruppi= gruppoLin + n.termini smooth
      if(nomiGruppi[1]=="Xlin") {
        pLin<- n.par[1]
        n.par<-n.par[-1] 
      } else {
        pLin<-0
      }
      H <- length(n.par) #n. smooths
      list.b <-tapply(b, id.coef, function(.x).x) 
      names(list.b)<- nomiGruppi
      if(any(vConc!=0)) list.b1 <-list.b #creato list.b1 per l'eventuale concav (infatti list.b viene sovrascritto dalla monot)
      id.Xlin <- match("Xlin", nomiGruppi, 0)
      P<- if(nrow(D.matrix)==0) crossprod(D.matrix) else crossprod(drop(pesiL1)*D.matrix[(1:length(pesiL1)),,drop=FALSE])
      #browser()
      
      if(any(vMonot!=0)){
        wMon<- lapply(1:H, function(.x) matrix(0,n.par[.x]-1,n.par[.x]))
        for(j in 1:H){ 
          if(attr(Bconstr[[j]], "constr.fit")) {
            list.b[[j + id.Xlin]] <-  drop(Bconstr[[j]] %*% list.b[[j + id.Xlin]])
            wMon[[j]]  <- crossprod((10^6*(abs(diff(list.b[[j+ id.Xlin]], diff=1))<.00001))*D1[[j]])
          } else {
            wMon[[j]]  <- crossprod((10^6*(abs(drop(D1[[j]]%*%list.b[[j+ id.Xlin]]))<.00001))*D1[[j]])
          }
        }
        wMon<-if(length(wMon)==1) wMon[[1]] else do.call("blockdiag", wMon)
        if(nomiGruppi[1]=="Xlin") wMon <-blockdiag(matrix(0, pLin, pLin), wMon) 
        P <- P + wMon
      }
      
      if(any(vConc!=0)){
        wConc<-lapply(1:H, function(.x) matrix(0,n.par[.x]-2,n.par[.x]))
        for(j in 1:H){ 
          if(attr(Bconstr[[j]], "constr.fit")) list.b1[[j+ id.Xlin]] <-  drop(Bconstr[[j]] %*% list.b1[[j+ id.Xlin]])
          wConc[[j]] <- crossprod((10^6*(abs(diff(list.b1[[j+ id.Xlin]], diff=2))<.00001))*D2[[j]])
        }
        if(any(vConc!=0)) wConc <-if(length(wConc)==1) wConc[[1]] else do.call("blockdiag", wConc) 
        if(nomiGruppi[1]=="Xlin") wConc <-blockdiag(matrix(0, pLin, pLin), wConc) 
        P <- P + wConc
      }
      
      #vincoli nc, dalla >=1.2-0
      if(df.nc) diag(P) <- diag(P)+10^6*DFvalues
      #if(df.nc){
      #DF.all<-cbind(obj$DF.NEG, obj$DF.POS) #NB DF.NEG non parte dal tau piu' piccolo, ma non e' un problema perche' si prende un tau alla volta..
      #if(paste(tau)%in%colnames(DF.all)) diag(P) <- diag(P) + 10^6*DF.all[1:length(b),paste(tau)]
      #}
      
#      fittedvalues <-drop(obj$fitted.values)[1:n]
      e<- drop(obj$residuals)[1:n]
      #Y<- obj$y[1:n]
      #X<-obj$x[1:n,,drop=FALSE] #design matrix
      #w<-ifelse(Y > fittedvalues, tau, 1-tau)
      w<-ifelse(e > 0, tau, 1-tau)
      w<-w/(abs(e)+.00001)
      XtWX<-crossprod(X*sqrt(w))
      A<-try(solve(XtWX+P, XtWX), silent=TRUE)
      if(class(A)[1]=="try-error") A<-solve(make.positive.definite(XtWX+P), XtWX)
      df.all<-diag(A)
      df.j<-tapply(df.all, id.coef, sum)
      if(is.matrix(df.j)) rownames(df.j)<-attr(id.coef, "nomi") else names(df.j)<-attr(id.coef, "nomi")
      r<-list(edf.all=df.all, edf.j=df.j)
      r
    }
    #-------------------------------------
    
    taus<-sort(taus)
    n<-length(y)

    #===========================
    #Se manca la base B, costruiscila
    #===========================
    if(is.null(B)) {
      #if(missing(x) || !is.matrix(x)) stop("'B' (list) or 'x' (matrix) have to be supplied")
      if(missing(x) || !(is.matrix(x)||is.data.frame(x))) stop("'B' (list) or 'x' (matrix) have to be supplied")
      B<-vector("list", length=ncol(x))      
      for(j in 1:ncol(x)) {
        B[[j]]<-bspline(x[,j], ndx=ndx[j], deg=deg[j])
        if(dropcList[j]) B[[j]]<-B[[j]][,-1]
        #C=contr.sum(ncol(B)) (o C=contr.treatment(ncol(B)))
        #C<-contr.treatment(ncol(B))
        #Bnew<-B%*%C
        #Pnew<- t(C)%*%Penalty%*%C
        if(decomList[j]) {
          Dj<-diff(diag(ncol(B[[j]])), diff = dif[j])
          B[[j]]<- B[[j]]%*%t(solve(tcrossprod(Dj),Dj)) #sono ncol(B[[j]])-vDiff[j]
          if(monotone[j]!=0) stop("something wrong in monotone/spline..") 
          #non c'? bisogno di costruire i termini fissi della spline perche' gi? sono in X
          #BFixed[[j]]<-poly(x[,j], degree=dif[j]-1, raw=TRUE)
          #colnames(BFixed[[j]])<- paste(nomiPS[j],"fix",1:(dif[j]-1),sep=".")
          #vDiff[j] <-0 #serve per ottenere una penalit? ridge quando si passa a build.R().
        }
        #per costruire una base per i vcm
        if(!is.null(nomiBy[j])) B[[j]]<- cbind(byVariabili[,j], byVariabili[,j]*B[[j]]) #modificato il 10/3/21. Se sopra togli una colonna della base devi metter l'interc
      }
    }
    #===========================
    if(missing(x)) plott<-0
    
    #deve diventare una matrice..
    all.p<-sapply(B, ncol)
    pSmooth <- sum(all.p)
    H <-length(B) #no. of smooth terms
    p <- pSmooth
    if(!is.null(X)) p<- p + ncol(X)
    # B <-matrix(unlist(B), n, pSmooth) #matrix(unlist(B), nrow=3,byrow = FALSE)
    # XB<-cbind(X,B)
    # p<-ncol(XB) #all coeffs, lineari + smooth
    if(!is.null(X)) {
      pLin <- ncol(X)
      id.interc <- match("(Intercept)", colnames(X), 0)
    } else {
      pLin <- 0
      id.interc <- 0
      }
    
    #browser()
    #attenzione. Qua se la matrice B comprende variabili da selezionare allora devi comunque consideraarle variabili lineari 
    #e quindi modificare i vincoli..
      
    if((pLin - id.interc)>0){ #se ci sono variabili lineari (oltre all'eventuale intercetta)
      if(adjX.constr){        
        colnamesX <- colnames(X)
        if(id.interc>0){ #se c'e' intercetta
          minX <- apply(X[,-1, drop=FALSE] , 2, min)
          names(minX)<- colnamesX[-1]
          X <- cbind(X[,1], apply(X[,-1, drop=FALSE], 2, function(.x) .x- min(.x)))
          } else {
          minX <- apply(X , 2, min)
          names(minX) <- colnamesX
          X <- apply(X, 2, function(.x) .x- min(.x))
          }
        all.max <- apply(X, 2, max)
        colnames(X) <- colnamesX
      } else {
          if( (length(taus)>1)) warning("shifting of linear covariates is requested for noncross")
          all.max <- apply(X, 2, max)
          minX <- NULL
      }
    } else {
      minX<-NULL
      }

    if(pLin!=(p-pSmooth)) stop("errore") #pSmooth=p1 #n. termini lineari

    all.edf<-matrix(NA,p,length(taus))
    group.edf<- matrix(NA, H+(pLin>0),length(taus))
    colnames(all.edf)<-colnames(group.edf)<-paste(taus)
    
    #browser()
    #un controllo se penMatrix e' fornita
    if(!is.null(penMatrix)){
      for(j in 1:H) {
        if(!is.null(penMatrix[[j]])){
          if(ncol(penMatrix[[j]])!=all.p[j]) stop(paste("wrong dimension (ncol) of pen.matrix for smooth #", j, sep=""))
        }
      }
    }
    
    #quando le basi sono non centrate la media e' 0.
    m.colmeansB <- lapply(colmeansB, function(.x) -.x)
    
    if(length(dropvcList)!=length(all.p)) stop("errore..") #non dovrebbe servire..
    id.intercVC <- tapply(rep(1:length(all.p), all.p), rep(1:length(all.p), all.p), function(.x)c(TRUE, rep(FALSE, length(.x)-1)))
    id.intercVC <- unlist(id.intercVC) & rep(dropvcList, all.p)
    
    id.intc.VC<- which((id.intercVC)) #le posizioni delle intercette dei VC (quando ci sono gruppi)
    exist.VC <-sum(vcList) #se >0 ci sono VC
    #n.righe = p
    #n.colonne = n. di termini lineari + n.smooth. Se non ci sono termini lineari la 1st colonna e' una colonna di FALSE..
    id.Matr <- sapply(0:length(all.p), function(.x) I(c(rep(0, pLin), rep(1:length(all.p), all.p))==.x))
    all.maxB<-as.list(rep(1,H))
    
    for(j in 1:H){
      if(ridgeList[j] && ps.matrix.list[j]) {
        cln <- colnames(B[[j]])
        B[[j]] <- apply(B[[j]], 2, function(.x) .x- min(.x))
        all.maxB[[j]] <-apply(B[[j]], 2, max)
        colnames(B[[j]]) <- cln
      }
    }
    #H <-length(B) #no. of smooth terms
    B <-matrix(unlist(B), n, pSmooth) #matrix(unlist(B), nrow=3,byrow = FALSE)
    XB<-cbind(X,B)
    
    id.start.tau<-which.min(abs(taus-0.5))
    lambda.tau.spec<-FALSE
    if(is.matrix(lambda)) {
      lambda.tau.spec<-TRUE
      lambda.matrix<-lambda
      colnames(lambda.matrix)<-paste(taus)
      lambda<-lambda.matrix[,id.start.tau]
    }

    lambdaM<-matrix(lambda, nrow=nrow(as.matrix(lambda)), ncol=length(taus)) #serve per il calcolo degli edf..
    colnames(lambdaM)<-paste(taus)
    
    D.list<-D.list.lambda<-D1<-D2<-vector("list", length=H)
    #browser()
    for(j in 1:H){
      #monot?
      if(monotone[j]!=0) {
        if(attr(Bconstr[[j]],"constr.fit")) { 
          D1[[j]] <-sign(monotone[j])*diff(diag(nrow(Bconstr[[j]])), diff=1) %*% Bconstr[[j]]
        } else {
          D1[[j]]<-if(dropcList[j]) sign(monotone[j])*(diff(diag(all.p[j]+1), diff=1)[,-1]) else sign(monotone[j])*diff(diag(all.p[j]), diff=1)
          #browser()
          #D1[[j]]<- sign(monotone[j])*diff(diag(all.p[j]), diff=1)
        }
      } else { 
        D1[[j]]<-matrix(0, all.p[j]-1, all.p[j]) 
      }
      #concave?
      if(concave[j]!=0) {
        if(attr(Bconstr[[j]],"constr.fit")) { 
          D2[[j]] <-sign(concave[j])*diff(diag(nrow(Bconstr[[j]])), diff=2) %*% Bconstr[[j]]
        } else {
          D2[[j]]<-if(dropcList[j]) sign(-concave[j])*(diff(diag(all.p[j]+1), diff=2)[,-1]) else sign(-concave[j])*diff(diag(all.p[j]), diff=2)
          #D2[[j]] <- sign(-concave[j])*diff(diag(all.p[j]), diff=2)
        }
      } else { 
        D2[[j]]<- matrix(0, all.p[j]-2, all.p[j]) 
      }
      
      D.list[[j]] <-build.D(var.pen[j], all.p[j]+1*dropcList[j], dif[j], lambda.ok=1, penMatrix[[j]], 
                            dropc=dropcList[j], decomp=decomList[j], dropvc = dropvcList[j], ridge = ridgeList[j])
      D.list.lambda[[j]]<- D.list[[j]]*lambda[j]
    }      
    
    R.monot<-if(length(D1)<=1) D1[[1]] else do.call("blockdiag",D1)
    R.monot<-cbind(matrix(0,nrow=nrow(R.monot), ncol=pLin), R.monot)
    if(any(concave!=0)) {
      R.conc<-if(length(D2)<=1) D2[[1]] else do.call("blockdiag",D2)
      R.conc<-cbind(matrix(0,nrow=nrow(R.conc), ncol=pLin), R.conc)
      R.monot<-rbind(R.monot, R.conc)
    }
    r.monot<-rep(0, nrow(R.monot))
    #matrice D *NON* comprende lambda
    D.matrix<-if(length(D.list)<=1) D.list[[1]] else do.call("blockdiag",D.list) 
    D.matrix<- blockdiag(diag(rep(0,pLin), ncol=pLin),D.matrix)
    #matrice D che comprende lambda. D.list.lambda prima era DD
    D.matrix.lambda<-if(length(D.list.lambda)<=1) D.list.lambda[[1]] else do.call("blockdiag",D.list.lambda) 
    D.matrix.lambda<- blockdiag(diag(rep(0,pLin), ncol=pLin), D.matrix.lambda)
    
    start.tau<-taus[id.start.tau]
    pos.taus<-taus[(taus-start.tau)>0]
    neg.taus<-taus[(taus-start.tau)<0]
    n.pos.taus<-length(pos.taus)
    n.neg.taus<-length(neg.taus)
    XB.orig<-XB
    
    if(any(lambda>0)){
      XB<- rbind(XB.orig, D.matrix.lambda)
    }
    if(lambda.ridge>0) XB<-rbind(XB, lambda.ridge*diag(ncol(XB))) #a small ridge penalty
    y<-c(y, rep(0,nrow(XB)-n))
    
    
    #apprx=TRUE
    
    #browser()
    
    if(it.j<=2){
      #o<-rq.fit.fnb(y=y, x=XB, tau=taus) #restituisce solo coef p x ntau
      #colnames(o$coefficients)<- paste(taus)
      #o$fitted.values <- XB.orig%*%o$coefficients
      #o$residuals<- sapply(1:length(taus), function(.x) y[1:n]-o$fitted.values[,.x])
      
      #all.rho <-sapply(1:length(taus), function(.x) sum(Rho(o$residuals[,.x], taus[.x])))
      id.coef.spline<-vector("list", length=H+1)
      id.coef.spline[[1]]<- if(ncol(X)>0) 1:ncol(X) else 0
      for(j in 1:H) id.coef.spline[[j+1]] <- max(id.coef.spline[[j]]) + seq_len(all.p[j])
      names(id.coef.spline)<-c("Xlin", names(lambda))
      id.df<-rep(1:length(id.coef.spline), sapply(id.coef.spline,length))
      attr(id.df, "nomi")<-c("Xlin", names(lambda))
      if(length(id.coef.spline[[1]])==1 && id.coef.spline[[1]]==0) {
        id.df<-id.df[-1]
        attr(id.df, "nomi")<-names(lambda)
      }
      
      all.rho <- taus
      all.sigma2u<-matrix(,H, length(taus))
      epss<-if(it.j==1) .01 else .001
      
      #browser()
      
      for(i in 1:length(taus)){
        o<-rq.fit.fnb(y=y, x=XB, tau=taus[i], eps=epss) 
        all.rho[i] <- sum(Rho(drop(o$residuals)[1:n], taus[i]))
        o$n<-n   #serve a edf.rqXB() 
        
        pesiL1 <- abs((D.matrix%*%o$coefficients))
        lambdasTutti<- rep(lambdaM[,i], sapply(D.list, nrow)) #controlla se serve -1*dropcList
        if(!is.null(pLin) && pLin>0) lambdasTutti<-c(rep(0, pLin),lambdasTutti)
        pesiL1<-sqrt(lambdasTutti/(pesiL1+.00001)) #(pesiL1[,k]+.00001)
        df.start.tau <- edf.rqXB(o, tau=taus[i], id.coef=id.df, vMonot = monotone, vConc=concave, pesiL1=pesiL1,
                               output=list(...)$df.option, df.nc=FALSE, DFvalues=NULL, D.matrix, Bconstr, D1, D2, 
                               Y=y[1:n], X=XB.orig)
        
        all.edf[,i]<- df.start.tau$edf.all
        group.edf[,i]<-df.start.tau$edf.j
        for(j in 1:H) all.sigma2u[j,i] <- sum(abs(D.list[[j]]%*%o$coefficients[id.coef.spline[[j+1]]]))
      }
      
      r<-list()
      r$edf.all<- all.edf
      r$edf.j <-group.edf
      r$id.coef <- id.df
      r$rho <- all.rho
      r$all.dev2u=all.sigma2u
      r$all.dev2e=all.rho #ll.sigma2e
      return(r)
    }
    
    #browser()

    o.start<-if(any(monotone!=0 | concave!=0)) rq.fit.fnc(x=XB, y=y, tau=start.tau, R=R.monot, r=r.monot, eps=myeps)
      else rq.fit(x=XB, y=y, tau=start.tau, method="fn") # eps=myeps) #rq.fit.pfn()?
    if(is.null(o.start$fitted.values)) o.start$fitted.values<- y - o.start$residuals
    #browser()
    #Attenzione.. guarda gli oggetti che restitruisce e che deve usare edf.rqXB(). Altrimenti non va..
    
    if(!all(is.finite(o.start$coefficients))) warning(paste("Some NA estimate in the quantile tau =", start.tau, sep="" ))
    
    o.start$n<-n
    o.start$rho<-sum(Rho(o.start$residuals[1:n], start.tau))
    id.coef.spline<-vector("list", length=H+1)
    id.coef.spline[[1]]<- if(ncol(X)>0) 1:ncol(X) else 0
    for(j in 1:H) id.coef.spline[[j+1]] <- max(id.coef.spline[[j]]) + seq_len(all.p[j])
    names(id.coef.spline)<-c("Xlin", names(lambda))
    #nomi.ok<-paste(names(lambda)[k],"ps",1:all.p[k],sep=".")
    sigma2u<-vector(length=H)
    for(j in 1:H) sigma2u[j]<- sum(abs(D.list[[j]]%*%o.start$coefficients[id.coef.spline[[j+1]]]))
    sigma2e<-if(err.rho) o.start$rho else sum(o.start$residuals[1:n]^2*abs(start.tau-I(o.start$residuals[1:n]<0))) #in realta' e' la dev..
    id.df<-rep(1:length(id.coef.spline), sapply(id.coef.spline,length))
    attr(id.df, "nomi")<-c("Xlin", names(lambda))
    if(length(id.coef.spline[[1]])==1 && id.coef.spline[[1]]==0) {
      id.df<-id.df[-1]
      attr(id.df, "nomi")<-names(lambda)
    }

    #--------calcolo edf
    pesiL1 <- abs((D.matrix%*%o.start$coefficients))
    #lambdaM ha tante colonne quanti i tau...
    lambdasTutti<- rep(lambdaM[,paste(start.tau)], sapply(D.list, nrow)) #controlla se serve -1*dropcList
    if(!is.null(pLin) && pLin>0) lambdasTutti<-c(rep(0, pLin),lambdasTutti)
    pesiL1<-sqrt(lambdasTutti/(pesiL1+.00001)) #(pesiL1[,k]+.00001)
    df.start.tau <- edf.rqXB(o.start, tau=start.tau, id.coef=id.df, vMonot = monotone, vConc=concave, pesiL1=pesiL1,
                             output=list(...)$df.option, df.nc=FALSE, DFvalues=NULL, D.matrix, 
                             Bconstr, D1, D2, Y=y[1:n], X=XB.orig)
    all.edf[,paste(start.tau)]<- df.start.tau$edf.all
    group.edf[,paste(start.tau)]<-df.start.tau$edf.j
    #-------------------

    #UNO o PIU' QUANTILI?
    if(length(taus)<=1){ #se length(taus)==1
      all.COEF<-as.matrix(o.start$coef) #era all.COEF<- o.start$coef 
      colnames(all.COEF)<-paste(taus)
      all.rho<-sum(Rho(o.start$residuals[1:n], start.tau)) 
      r<-list(coefficients=all.COEF, x=XB, rho=all.rho, #tolto df=all.df,
              fitted.values=o.start$fitted.values[1:n],residuals=o.start$residuals[1:n],
              dev2u=sigma2u, dev2e=sigma2e)
      all.dev2u=matrix(sigma2u, ncol=1, dimnames=list(names(lambda), paste(taus)))
      all.dev2e=matrix(sigma2e, ncol=1, dimnames=list(NULL, paste(taus)))
      r$all.dev2u<-all.dev2u
      r$all.dev2e<-all.dev2e
      
    } else { #se length(taus)>1
      DF.NEG <- DF.POS<- NULL
      Ident<-diag(p)
      COEF.POS<-COEF.NEG<-FIT.POS<-FIT.NEG<-RES.POS<-RES.NEG<-NULL
      df.pos.tau<-df.neg.tau<-rho.pos.tau<-rho.neg.tau<-NULL
      sigma2u.pos.tau<-sigma2u.neg.tau<-NULL
      sigma2e.pos.tau<-sigma2e.neg.tau<-NULL
      
      
      id.interc.Constr <- id.interc>0
      if(id.interc.Constr) valueNcrosInterc <- unlist(m.colmeansB)
      #indici dei coeff delle diverse basi, anche se VC
      idsmoothInter <- setdiff(unlist(id.coef.spline),id.coef.spline$Xlin)
      
      #===========================================================================

      #===========================================================================
      if(n.pos.taus>0){
        rho.pos.tau <-df.pos.tau <- vector(length=n.pos.taus)
        COEF.POS<-matrix(,ncol(XB),n.pos.taus)
        colnames(COEF.POS)<-paste(pos.taus)
        b.start<-o.start$coef
        if(any(monotone!=0 | concave!=0)){
          RR<-rbind(Ident, R.monot)
          rr<-c(b.start + eps, r.monot)
        } else {
          RR<- Ident
          rr<- b.start + eps
        }
        
        #browser()
        
        #DF per NONcrossing
        DF.POS <- matrix(,nrow(RR),n.pos.taus)
        colnames(DF.POS) <- paste(pos.taus)
        
        if(id.interc.Constr){ #se c'e' interc con smooth (incompleti, cioe' chiamati con dropc=FALSE)
          RR[1, idsmoothInter] <- valueNcrosInterc #valueNcrosInterc sono le "-colmeansB"
          rr[1] <- b.start[1] + sum(valueNcrosInterc * b.start[idsmoothInter]) #b.start[1]+sum((-all.means.B)*b.start[id.smooth])
        }
        #n.righe = p
        #id.Matr: p x "n. di termini lineari + n.smooth". Se non ci sono termini lineari la 1st colonna e' una colonna di FALSE..

        if((pLin - id.interc)>0){ #se ci sono variabili lineari (oltre all'eventuale intercetta)
#          colnamesX <- colnames(X)
#          all.max<-apply(X,2,max)
          if(id.interc>0){ #se c'e' intercetta
            RR[c(FALSE,id.Matr[-1,1]),1]<-1 #aggiungi gli 1 nella prima *colonna* in corrispondenza dei termini lineari
            rr[2:pLin] <- b.start[1] + b.start[2:pLin]*all.max[2:pLin] + eps
            
          } else {
            rr[1:pLin] <- b.start[1:pLin]*all.max[1:pLin] + eps
#            minX <- apply(X[,-1, drop=FALSE] ,2, min)
#            names(minX)<-colnames(X)
#            X <- apply(X, 2, function(.x) .x- min(.x))
          }
          diag(RR)[1:pLin]<- all.max  #aggiungi i max sulla diag principale
 #         colnames(X) <- colnamesX
        }
        
        #browser()
        
        for(j in 1:H) {
          if(ps.matrix.list[j]) {
            RR[id.Matr[,j+1], 1] <- 1
            diag(RR)[id.Matr[,j+1]]<-all.maxB[[j]]
            
            if(id.interc>0) {
              rr[id.Matr[,j+1]] <- b.start[1] + b.start[id.Matr[,j+1]]*all.maxB[[j]] + eps
            } else {
              rr[id.Matr[,j+1]] <- b.start[id.Matr[,j+1]]*all.maxB[[j]] + eps
            }
          }
        }
        
        FIT.POS<-RES.POS<-matrix(,n,n.pos.taus)
        sigma2u.pos.tau<-matrix(,H,n.pos.taus)
        sigma2e.pos.tau<-vector(length=n.pos.taus)
        
        for(i in 1:n.pos.taus){
          #AGGIORNA XB PER INCLUDERE i tau-specific lambda
          if(lambda.tau.spec){
            lambda<-lambda.matrix[,paste(pos.taus[i])] #lambda.matrix[,i]
            for(j in 1:H) D.list.lambda[[j]]<- D.list[[j]]*lambda[j]
            D.matrix.lambda<-if(length(D.list.lambda)<=1) D.list.lambda[[1]] else do.call("blockdiag",D.list.lambda) 
            D.matrix.lambda<-blockdiag(diag(rep(0,pLin), ncol=pLin),D.matrix.lambda)
            if(any(lambda>0)) XB<-rbind(XB.orig, D.matrix.lambda)
            if(lambda.ridge>0) XB<-rbind(XB, lambda.ridge*diag(ncol(XB))) #a small ridge penalty
          }
          o<-rq.fit(x=XB, y=y, tau=pos.taus[i], method="fnc", R=RR, r=rr)
          #o<-rq.fit.fnc(x=XB, y=y, tau=pos.taus[i], method="fnc", R=RR, r=rr, eps=myeps)
          if(!all(is.finite(o$coefficients))) warning(paste("Some NA estimate in the quantile tau =", pos.taus[i], sep="" ))
          
          o$rho<-sum(Rho(o$residuals[1:n], pos.taus[i]))
          #estrai fitted e residuals
          FIT.POS[,i]<-o$fitted.values[1:n]
          RES.POS[,i]<-o$residuals[1:n]
          #estrai la f. obiettivo
          #df.pos.tau[i] <- sum(abs(o$residuals[1:n])<=.000001)
          rho.pos.tau[i] <- o$rho #sarebbe sum(Rho(o$residuals[1:n], pos.taus[i]))
          #estrai sigma2u
          for(j in 1:H) sigma2u.pos.tau[j,i]<- sum(abs(drop(D.list[[j]]%*%o$coefficients[id.coef.spline[[j+1]]])))

          #estrai sigma2e
          sigma2e.pos.tau[i]<-o$rho #if(err.rho) o$rho else sum(o$residuals[1:n]^2*abs(pos.taus[i]-I(o$residuals[1:n]<0)))
          
          #estrai i coeff per "aggiornare" i vincoli rr
          b.start<-o$coef
          COEF.POS[,i]<-b.start
          #NON ANCORA IMPLEMENTATO: confrontare dove i vincoli NONCROSSING (MA ANCHE MONOT/CONC) sono attivi
          #vedi come considerare l'epsilon.. Inoltre la matrice RR include anche monot/conc
          #in particolare i primi RR[1:p,]b.start rr[1:p] sono per il NONCROSSING
          DF.POS[,i] <- 1*(drop(RR%*%b.start-rr)<=1e-8)

          rr<- if(any(monotone!=0 | concave!=0)) c(b.start+eps, r.monot) else b.start + eps
          
          if(id.interc.Constr){ #se ci sono smooth con interc
            #  rr[1] <- sum(valueNcrosInterc * b.start[idsmoothInter]) #b.start[1]-sum(all.means.B*b.start[id.smooth])
            rr[1] <- b.start[1] + sum(valueNcrosInterc * b.start[idsmoothInter])
          }
          if((pLin - id.interc)>0) { #se ci sono termini lineari oltre l'intercetta
            if(id.interc>0){
              rr[2:pLin] <- b.start[1] + b.start[2:pLin]*all.max[2:pLin] + eps
            } else {
              rr[1:pLin] <- b.start[1:pLin]*all.max[1:pLin] + eps
            }
          }

          for(j in 1:H) {
            if(ps.matrix.list[j]) {
              if(id.interc>0) {
                rr[id.Matr[,j+1]] <- b.start[1] + b.start[id.Matr[,j+1]]*all.maxB[[j]] + eps
              } else {
                rr[id.Matr[,j+1]] <- b.start[id.Matr[,j+1]]*all.maxB[[j]] + eps
              }
            }
          }

          o$n<-n
          #--------calcolo edf
          pesiL1<-abs((D.matrix%*%o$coefficients))

          #lambdaM ha tante colonne quanti i tau...
          lambdasTutti<- rep(lambdaM[,paste(pos.taus[i])], sapply(D.list, nrow)) #controlla se serve -1*dropcList
          if(!is.null(pLin) && pLin>0) lambdasTutti<-c(rep(0, pLin),lambdasTutti)
          pesiL1<-sqrt(lambdasTutti/(pesiL1+.00001)) #(pesiL1[,k]+.00001)
          df.pos.tau <- edf.rqXB(o, tau=pos.taus[i], id.coef=id.df, vMonot = monotone, vConc=concave, pesiL1=pesiL1,
                                   output=list(...)$df.option, df.nc=list(...)$df.nc, DFvalues=DF.POS[,i], 
                                 D.matrix, Bconstr, D1, D2, Y=y[1:n], X=XB.orig)
          all.edf[,paste(pos.taus[i])]<- df.pos.tau$edf.all
          group.edf[,paste(pos.taus[i])]<-df.pos.tau$edf.j
          #-------------------
          }#end for
      }#end if(n.pos.taus>0)
      ### PER I tau "negativi"..
      
      #browser()
      
      if(n.neg.taus>0){
        rho.neg.tau <-df.neg.tau <- vector(length=n.pos.taus)
        COEF.NEG<- matrix(,ncol(XB),n.neg.taus)
        colnames(COEF.NEG)<-paste(neg.taus)
        b.start<-o.start$coef
        neg.taus<-sort(neg.taus,TRUE)
        if(any(monotone!=0 | concave!=0)){
          Ident<-diag(p)
          RR<-rbind(-Ident, R.monot)
          rr <-c(-b.start - eps, r.monot)
        } else {
          RR<- -Ident
          rr<- -b.start - eps
        }
        #DF per NONcrossing
        DF.NEG<- matrix(,nrow(RR),n.neg.taus)
        colnames(DF.NEG)<-paste(neg.taus)
        
        if(id.interc.Constr){ #se ci sono smooth con interc
          RR[1, idsmoothInter] <- -valueNcrosInterc #valueNcrosInterc sono le "-colmeansB"
          rr[1] <- -(b.start[1] + sum(valueNcrosInterc * b.start[idsmoothInter])) #b.start[1]+sum((-all.means.B)*b.start[id.smooth])
        }
        
        if((pLin - id.interc)>0){ #se ci sono variabili lineari (oltre all'eventuale intercetta)
          #colnamesX <- colnames(X)
          #all.max<-apply(X,2,max)
          if(id.interc>0){ #se c'e' intercetta
            RR[c(FALSE,id.Matr[-1,1]),1]<- -1 #aggiungi i -1 nella prima *colonna* in corrispondenza dei termini lineari
            rr[2:pLin] <- -(b.start[1] + b.start[2:pLin]*all.max[2:pLin] + eps)
#            minX <- apply(X[,-1, drop=FALSE] ,2, min)
#            names(minX)<-colnames(X)[-1]
#            X <- cbind(X[,1], apply(X[,-1, drop=FALSE], 2, function(.x) .x- min(.x)))
          } else {
            rr[1:pLin] <- -(b.start[1:pLin]*all.max[1:pLin] + eps)
#            minX <- apply(X[,-1, drop=FALSE] ,2, min)
#            names(minX)<-colnames(X)
#            X <- apply(X, 2, function(.x) .x- min(.x))
          }
          diag(RR)[1:pLin]<- -all.max  #aggiungi i max sulla diag principale
#          colnames(X) <- colnamesX
        }
        
        
        for(j in 1:H) {
          if(ps.matrix.list[j]) {
            RR[id.Matr[,j+1], 1] <- -1
            diag(RR)[id.Matr[,j+1]]<- -all.maxB[[j]]
            
            if(id.interc>0) {
              rr[id.Matr[,j+1]] <- -(b.start[1] + b.start[id.Matr[,j+1]]*all.maxB[[j]] + eps)
            } else {
              rr[id.Matr[,j+1]] <- -(b.start[id.Matr[,j+1]]*all.maxB[[j]] + eps)
            }
          }
        }
        

        FIT.NEG<-RES.NEG<-matrix(,n,n.neg.taus)
        sigma2u.neg.tau<-matrix(,H,n.neg.taus)
        sigma2e.neg.tau<-vector(length=n.neg.taus)
        
        #browser()
        for(i in 1:n.neg.taus){
          #AGGIORNA XB PER INCLUDERE i tau-specific lambda
          
          #if(i==2) browser()
          
          if(lambda.tau.spec){
            lambda <- lambda.matrix[,paste(neg.taus[i])]
            for(j in 1:H) D.list.lambda[[j]]<- D.list[[j]]*lambda[j]
            D.matrix.lambda<-if(length(D.list.lambda)<=1) D.list.lambda[[1]] else do.call("blockdiag",D.list.lambda) 
            D.matrix.lambda<-blockdiag(diag(rep(0,pLin), ncol=pLin),D.matrix.lambda)
            if(any(lambda>0)) XB<-rbind(XB.orig, D.matrix.lambda)
            if(lambda.ridge>0) XB<-rbind(XB, lambda.ridge*diag(ncol(XB))) #a small ridge penalty
          }
          o<-rq.fit(x=XB,y=y,tau=neg.taus[i],method="fnc",R=RR,r=rr)
          #o<-rq.fit.fnc(x=XB,y=y,tau=neg.taus[i], R=RR, r=rr, eps=myeps)
          ##########==============================================================
          
          if(!all(is.finite(o$coefficients)))  warning(paste("At least one estimate is NA in the quantile curve tau = ", neg.taus[i], sep="" ))
          o$rho<-sum(Rho(o$residuals[1:n], neg.taus[i]))
          FIT.NEG[,i]<-o$fitted.values[1:n]
          RES.NEG[,i]<-o$residuals[1:n]
          #df.neg.tau[i] <- sum(abs(o$residuals[1:n])<=.000001)
          rho.neg.tau[i] <- o$rho #sum(Rho(o$residuals[1:n], neg.taus[i]))
          #estrai sigma2u
          for(j in 1:H){
            f.coef.pen <- abs(drop(D.list[[j]]%*%o$coefficients[id.coef.spline[[j+1]]]))
            sigma2u.neg.tau[j,i]<- if(u.rho) sum(f.coef.pen)
                else sum(drop(D.list[[j]]%*%o$coefficients[id.coef.spline[[j+1]]])^2)
          }
          #estrai sigma2e
          
          
          
          sigma2e.neg.tau[i]<-if(err.rho) o$rho else sum(o$residuals[1:n]^2*abs(neg.taus[i]-I(o$residuals[1:n]<0)))
          #estrai i coeff per aggiornare i vincoli
          b.start<-o$coef
          COEF.NEG[,i]<-b.start
          DF.NEG[,i] <- 1*(drop(RR%*%b.start-rr)<=1e-8)
          # rr<- if(monotone!=0) c(-b.start+eps, rep(0,p1-1) ) else -b.start+eps
          rr<- if(any(monotone!=0 | concave!=0)) c(-b.start+eps, r.monot) else -b.start + eps
          
          if(id.interc.Constr){ #se ci sono smooth con interc
            #  rr[1] <- sum(valueNcrosInterc * b.start[idsmoothInter]) #b.start[1]-sum(all.means.B*b.start[id.smooth])
            rr[1] <- -(b.start[1] + sum(valueNcrosInterc * b.start[idsmoothInter]))
          }
          
          if((pLin - id.interc)>0) { #se ci sono termini lineari oltre l'intercetta
            if(id.interc>0){
              rr[2:pLin] <- -(b.start[1] + b.start[2:pLin]*all.max[2:pLin] + eps)
            } else {
              rr[1:pLin] <- -(b.start[1:pLin]*all.max[1:pLin] + eps)
            }
          }
          
          for(j in 1:H) {
            if(ps.matrix.list[j]) {
              if(id.interc>0) {
                rr[id.Matr[,j+1]] <- -(b.start[1] + b.start[id.Matr[,j+1]]*all.maxB[[j]] + eps)
              } else {
                rr[id.Matr[,j+1]] <- -(b.start[id.Matr[,j+1]]*all.maxB[[j]] + eps)
              }
            }
          }
          
          o$n<-n
          #--------calcolo edf
          pesiL1<-abs((D.matrix%*%o$coefficients))
          #lambdaM ha tante colonne quanti i tau...
          lambdasTutti<- rep(lambdaM[,paste(neg.taus[i])], sapply(D.list, nrow)) #controlla se serve -1*dropcList
          if(!is.null(pLin) && pLin>0) lambdasTutti<-c(rep(0, pLin),lambdasTutti)
          pesiL1<-sqrt(lambdasTutti/(pesiL1+.00001)) #(pesiL1[,k]+.00001)
          df.neg.tau <- edf.rqXB(o, tau=neg.taus[i], id.coef=id.df, vMonot = monotone, vConc=concave, pesiL1=pesiL1,
                                 output=list(...)$df.option, df.nc=list(...)$df.nc, DFvalues=DF.NEG[,i], 
                                 D.matrix, Bconstr, D1, D2, Y=y[1:n], X=XB.orig)
          all.edf[,paste(neg.taus[i])]<- df.neg.tau$edf.all
          group.edf[,paste(neg.taus[i])]<-df.neg.tau$edf.j
          #-------------------

        }#end for
      }#end if(n.neg.taus>0)
      #-------------------------------
      #browser()
      all.COEF<-cbind(COEF.NEG[,n.neg.taus:1], o.start$coef, COEF.POS)
      colnames(all.COEF)<-paste(taus)
      all.FIT<-cbind(FIT.NEG[,n.neg.taus:1,drop=FALSE], o.start$fitted.values[1:n], FIT.POS)
      colnames(all.FIT)<-paste(taus)
      all.RES<-cbind(RES.NEG[,n.neg.taus:1,drop=FALSE], o.start$residuals[1:n], RES.POS)
      colnames(all.RES)<-paste(taus)
      
      
      all.rho<-c(rho.neg.tau[n.neg.taus:1], o.start$rho , rho.pos.tau)
      all.sigma2e<-c(sigma2e.neg.tau[n.neg.taus:1], sigma2e, sigma2e.pos.tau)
      all.sigma2u <-cbind(sigma2u.neg.tau[,n.neg.taus:1,drop=FALSE], sigma2u, sigma2u.pos.tau)
      colnames(all.sigma2u)<-paste(taus)      
      rownames(all.sigma2u)<-names(lambda)
      r<-list(coefficients=all.COEF,x=XB, rho=all.rho, fitted.values=all.FIT, residuals=all.RES,
              dev2u=sigma2u, dev2e=sigma2e, all.dev2u=all.sigma2u, all.dev2e=all.sigma2e, DF.NEG=DF.NEG, DF.POS=DF.POS) #
    } #end se length(taus)>1
    #browser()
    #Attenzione XB include nella 2nd parte la penalizzazione relativa all'*ultima* curva stimata 
    #    (quindi puo' essere un "problema" se i lambda sono diversi per quantili)
    r$edf.all<- all.edf
    r$edf.j <-group.edf
      
    r$id.coef <- id.df
    r$D.matrix<-D.matrix
    #cosa mettere?:  D.list, D.list.lambda, D.matrix, D.matrix.lambda
    r$nrowDlist<- sapply(D.list, nrow)
    r$D.matrix.lambda<-D.matrix.lambda
    r$pLin<-pLin      
    r$R.monot<-R.monot
    r$lambda<- if(lambda.tau.spec) lambda.matrix else lambda
    if(any(monotone!=0)) r$D1 <-D1
    if(any(concave!=0)) r$D2 <-D2
    if(!is.null(minX)) r$minX<- minX
    pesiPen <- NULL
    #browser()
    if(any(adList>0)) for(j in 1:H) pesiPen[[j]]<- abs(D.list[[j]]%*%all.COEF[id.coef.spline[[j+1]],,drop=FALSE])
    r$pesiPen <- pesiPen
    return(r)
  }

