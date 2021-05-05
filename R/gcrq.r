gcrq <-
function(formula, tau=c(.1,.25,.5,.75,.9), data, subset, weights, na.action, transf=NULL,
    y=TRUE, n.boot=0, eps=0.001,  #ho rimosso cv=FALSE, interc=TRUE,
    display=FALSE, method=c("REML","ML"), df.opt=2, 
    lambda0=.1, h=0.8, lambda.max=2000, tol=0.01, it.max=20, single.lambda=TRUE, 
    foldid=NULL, nfolds=10, lambda.ridge=0, sgn.constr=NULL, adjX.constr=TRUE,
    contrasts=NULL,
    ...){
#single.lambda: if TRUE (and multiple quantile curves are being estimated) a common spar across the tau values is computed
#sett 19: aggiunto df.opt=4 (tutto approx quadratica..)
  
#Ripresa il 16/11/16.. SPERIAMO CHE SIA LA VOLTA BUONA!!! :-)).. ottimismo ingiustificato... :-|

#keep in mind lambda=A/B, where A <-> err.rho   and   B<-> u.rho
#if err.rho =TRUE A=sum(|e|w)/Redf, where e=residual, w=weight depending on tau. Redf=n if method="ML" or Redf=n-edf if method="REML" 
#  edf are computed by 'df.option' (see below)
#if u.rho =TRUE A=sum(|Db|)/edf, where D=penalty, b=pen coeff and edf depend on 'df.option' (see below)
#if laplace=TRUE 'err.rho' is ignored, and A = Laplace dev.st which is based on the scale estimate: "checkLoss/n" (if method="ML") or 
#  "checkLoss/(n-edf) (if method="REML") 
#NB usa laplace=FALSE, err.rho=TRUE e u.rho=TRUE per rapporto dei parametri di scala con df regolati da method="ML"/"REML" per i df del 
#   numeratore e 'df.option' for the denominator df. Use method="ML" and df.option=1 for a *true* ML approach
#
#laplace.df="ML" quando laplace=TRUE come ottenere la stima del parametro di scala? dividendo per "n" o "n-edf"?
#   df.option
#   4 = come 3 ma tutto basato sull'approx dei valori assoluti
#   3 = the trace of the pseudo-hat matrix based on the smooth approximation (for REML-like estimate)
#   2 = the number of non-null (penalized) coefficients (for REML-like estimate)
#   1 = the number of all (penalized) coefficients (for ML-like estimate)
# Se 1 o 2, edf.j (e quindi il metodo print) vengono restituiti i df dei termini soltanto;
#l'intercetta non viene conteggiata,  sebbene inserita nel modello. Se 3 invece viene conteggiata.
#lambda.ridge: se si vuole mettere l'intercetta o ci sono multiple termini smooth deve essere >0 (.001) altrimenti la matrice non e' invertibile..
#lambda0: the starting value for the lambdas to be estimated
#g: per il calcolo dell'approx smooth
#h: il passo per il calcolo dei lambda
#it.max per il calcolo di lambda
#tol: la tol per il calcolo di lambda (.001 provides stabler estimates, but .01 suffices (in my experience))
#use.rho: quale stima di sigma usare?
#display: se visualizzare (meaningfull only if any spar is being estimated)
#interc: it is set TRUE if there is any ps() term with decompose=TRUE
#=================================
#all.perc fa riferimento alla stima di piu' percentili... 
#==>  Usare un lambda diverso per ogni tau o un unico tau? DA VEDERE 
#all.perc: quando tau e' un vettore, mean(fit$all.dev2e)/((n-sum(edf.j))^2) se all.perc=TRUE 
#(ovvero usa le info di tutti i percentili) oppure usa solo la dev.u di tau uguale a, o piu vicino a, 0.5.
#
#Growth Charts via QR
#**************weights??
#se cv=TRUE restituisce anche una componente 'cv' che e' una matrice di n.righe=n.valori di lambda e colonne nfolds
#eps in control??
#foldid, nfold usati se lambda in ps() e' un vettore
#--------------------
#  if(!single.lambda) {
#    single.lambda=TRUE
#    warning(" 'single.lambda=TRUE' is used")
#  }
dal<-function (x, mu = 0, sigma = 1, tau = 0.5, log = FALSE) {
#preso da  package lqmm 1.5
    ind <- ifelse(x < mu, 1, 0)
    val <- tau * (1 - tau)/sigma * exp(-(x - mu)/sigma * (tau - ind))
    if(log) 
      log(val) 
        else val
}
#--------------------
varAL<-function (sigma, tau){
#preso da  package lqmm 1.5
    sigma^2 * (1 - 2 * tau + 2 * tau^2)/((1 - tau)^2 * tau^2)
}
#--------------------
#-------------------------------------
blockdiag <- function(...) {
  args <- list(...)
  nc <- sapply(args,ncol)
  cumnc <- cumsum(nc)
  ##  nr <- sapply(args,nrow)
  ## NR <- sum(nr)
  NC <- sum(nc)
  rowfun <- function(m,zbefore,zafter) {
    cbind(matrix(0,ncol=zbefore,nrow=nrow(m)), #as.matrix.csr(
          m, 
          matrix(0,ncol=zafter,nrow=nrow(m)))
  }
  ret <- rowfun(args[[1]],0,NC-ncol(args[[1]]))
  for (i in 2:length(args)) {
    ret <- rbind(ret,rowfun(args[[i]],cumnc[i-1],NC-cumnc[i]))
  }
  ret
}
#-------------------------------------
edf.rq<-function(obj, tau, g=.2, id.coef, return.all=FALSE, vMonot, vConc, pesiL1, output, return.matrix=FALSE){
  #output: which df should be returned?
  #   4 = come 3, ma tutto (anche pen) attraverso su approx quadratica basata sui valori assoluti
  #   3 = the trace of the pseudo-hat matrix based on the smooth approximation (for REML-like estimate)
  #   2 = the number of non-null (penalized) coefficients (for REML-like estimate)
  #   1 = the number of all (penalized) coefficients (for ML estimate)
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
  #if(length(tau)>1) stop("multiple tau not allowed in edf.rq")
      interc<-"(Intercept)"%in% rownames(as.matrix(obj$coef))
      if(missing(id.coef)) id.coef <- obj$id.coef
      nomiGruppi<-attr(id.coef,"nomi") #"Xlin"(eventuale) + nomi smooth
      b<-if(is.matrix(obj$coef)) obj$coef[,paste(tau)] else obj$coef
      n.par<-tapply(id.coef, id.coef, length) #num. dei parametri di ciascun "gruppo" lineari, 1?smooth, 2?smooth,..
      if(output<=2){
          df.j<-n.coef.pen<- vNdx +vDeg - vDiff # un vettore - riferito ai termini smooth
          #if("Xlin" %in% nomiGruppi) df.j<-c(n.coef.pen)
          if(output>1) { #output==2 #QUA DEVI FARE LA PRE-MOLTIPLICAZIONE PER LE MATRICI D
            coef.pen<- obj$D.matrix%*%b #NB comprende termini param (con valori che sono 0)
            if(!is.null(obj$pLin) && obj$pLin>0) coef.pen<-coef.pen[-(1:obj$pLin)]
            id.coef.pen<-rep(1:length(n.coef.pen), n.coef.pen)
            df.j<-tapply(coef.pen,id.coef.pen, function(.x)sum(abs(.x)>.000001) ) #df of penalized coeffs 
          }
          #if(obj$pLin>0) df.j<-c(obj$pLin, df.j)
          df.j<- df.j + (vDiff-1) #*total* df of each smooth term - DIVIDERE PER QUESTO? NB (vDiff-1) ? il n. di coef unpen!
          if(!is.null(obj$pLin) && obj$pLin>0) df.j<-c(obj$pLin, df.j)          
          names(df.j)<-attr(id.coef, "nomi")
          return(df.j)
        }   
      
      n.Gruppi<-length(n.par)#length(n.par) ? il numero dei gruppi= gruppoLin + n.termini smooth
      delta.b<-tapply(b, id.coef, diff, simplify = FALSE)
      delta2.b<-tapply(b, id.coef, diff, simplify = FALSE, differences=2)
      names(delta.b)<- names(delta2.b)<-nomiGruppi 
      wConc<-wMon<-vector("list",n.Gruppi) #n. dei gruppi
      names(wMon) <- names(wConc) <-nomiGruppi
      istart<- if(nomiGruppi[1]=="Xlin") 2 else 1
      if(nomiGruppi[1]=="Xlin") {
        vMonot<-c(0,vMonot)
        vConc<-c(0,vConc)
        }
      if(length(vMonot)!=n.Gruppi || length(vConc)!=n.Gruppi) stop("Errore nella dim")
      for(i in istart:n.Gruppi){
            #per la monot
            Db.i<-delta.b[[i]]
            wMon[[i]]<-rep(0,length(Db.i))
            #if(vMonot[i]!=0) wMon[[i]] <- ifelse(rep(vMonot[i], length(Db.i))==1,1*(Db.i<=10e-8),1*(Db.i>=10e-8))
            if(vMonot[i]!=0) wMon[[i]] <- 1*(abs(Db.i)<=10e-8)
            D1<-diff(diag(length(Db.i)+1), diff=1)
            wMon[[i]]<-(10^6)*crossprod(wMon[[i]]*D1) #NON ? diagonale!!!
            #per la conc
            D2b.i<-delta2.b[[i]]
            wConc[[i]]<-rep(0,length(D2b.i))
            #if(vConc[i]!=0) wConc[[i]] <- ifelse(rep(vConc[i], length(D2b.i))==1,1*(D2b.i<=10e-8),1*(D2b.i>=10e-8))
            if(vConc[i]!=0) wConc[[i]] <- 1*(abs(D2b.i)<=10e-8)
            D2<-diff(diag(length(D2b.i)+2), diff=2)
            wConc[[i]]<-(10^6)* crossprod(wConc[[i]]*D2)
      }
      if(nomiGruppi[1]!="Xlin") {
        wMon <-c(1, wMon)
        wConc<-c(1, wConc)
      }
      wMon[[1]]<- matrix(0,0,0)
      wConc[[1]]<-matrix(0,0,0)
      
      wMon<-do.call("blockdiag",wMon)
      wConc<-do.call("blockdiag",wConc)
      if(nomiGruppi[1]=="Xlin") {
          #wMon <-cbind(as.matrix.csr(matrix(0,nrow=nrow(wMon),  ncol=obj$pLin)), as.matrix.csr(wMon) )
          #wConc<-cbind(as.matrix.csr(matrix(0,nrow=nrow(wConc), ncol=obj$pLin)), as.matrix.csr(wConc))
        wMon <-cbind(matrix(0,nrow=nrow(wMon),  ncol=obj$pLin), wMon ) 
        wConc<-cbind(matrix(0,nrow=nrow(wConc), ncol=obj$pLin), wConc)
              }          
      if(return.matrix){
        r<-list(mon=wMon, conc=wConc)
        r
      }
      
      #total.D<- obj$x[-(1:n),,drop=FALSE] #comprende *lambda* (ed eventualmente anche la matrice identita' della pen ridge)
      ## naturalmente P ottenuta sotto con crossprod contiene lambda^2 e quindi si deve correggere..
      ##P<- if(nrow(total.D)==0) crossprod(total.D) else crossprod(pesiL1*total.D[(1:length(pesiL1)),,drop=FALSE])
      ##matrice di penalita L1 : P<-crossprod(diag(1/sqrt(abs(diff(b, diff=1))))%*%D)
      
      total.D<-obj$D.matrix #NON comprende lambda
      P<- if(nrow(total.D)==0) crossprod(total.D) else crossprod(pesiL1*total.D[(1:length(pesiL1)),,drop=FALSE])
      
      #----
      #Se ci sono vincoli di monoton e concav
      P<- P +crossprod(wMon) + crossprod(wConc) #? megl
      
      X<-obj$x[1:n,,drop=FALSE] #design matrix
      e<-if(is.matrix(obj$residuals)) obj$residuals[,paste(tau)] else obj$residuals
      n<-length(e)
      
      if(output==3){
        g <- quantile(abs(e)[abs(e)>1e-10], probs=g, names=FALSE)
        g2<- 0
        g1<- -g*tau
        g3<- g*(1-tau)
        #v1<- I(e<g1)
        v2<- 1*I(e>=g1 & e<=g2)
        v3<- 1*I(e>=g2 & e<=g3)
        #v4<- I(e>g3)
        #w<- v2*(1-tau)/(g*tau) + v3*tau/(g*(1-tau))
        w<- v2*(1-tau)/(tau) + v3*tau/((1-tau))
      } else { #se output=4
        fittedvalues<-if(is.matrix(obj$fitted.values)) obj$fitted.values[,paste(tau)] else obj$fitted.values
        w<-ifelse(Y > fittedvalues, tau, 1-tau)
        w<-w/(abs(e)+.00001)
      }
      XtWX<-crossprod(X*sqrt(w))
      H<-try(solve(XtWX+P, XtWX), silent=TRUE)
      if(class(H)[1]=="try-error") H<-solve(make.positive.definite(XtWX+P), XtWX)
      df.all<-diag(H)
      if(return.all) return(df.all)
      df.j<-tapply(df.all, id.coef, sum)
      if(is.matrix(df.j)) rownames(df.j)<-attr(id.coef, "nomi") else names(df.j)<-attr(id.coef, "nomi")
      df.j
    }
#--------------------

#------------
bspline <- function(x, ndx, xlr = NULL, knots=NULL, deg = 3, deriv = 0, outer.ok=FALSE) {
  # x: vettore di dati
  # xlr: il vettore di c(xl,xr)
  # ndx: n.intervalli in cui dividere il range
  # deg: il grado della spline
  # Restituisce ndx+deg basis functions per ndx-1 inner nodi
  #ci sono "ndx+1" nodi interni + "2*deg" nodi esterni
  #require(splines)
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
  ##### ====================================================================================
  ##### INIZIO FUNZIONE
    #settiamo qualche opzione che forse non ? molto utile..
    
    all.perc=FALSE
    g=0.20 #usato solo se df.opt=3 (comunque sconsigliato)
    
    if(!(df.opt %in% 1:2)) {
      df.opt<-2
      warning("Meaningless value in 'df.opt'.. set to 2", call.=FALSE)
    }
    df.option<-2*df.opt
    laplace<-sd.L<-FALSE 
    err.rho<-TRUE #per il numeratore
    u.rho<-TRUE #per il denominatore
    penL1.df<-TRUE #usare sempre pen L1 per calcolare i df nella hat.matrix..  (se df.opt>=3)
    method<-match.arg(method)
    if(method=="ML") df.option<-1
    #check if the formula includes 'invalid' interactions with ps()
    s1<-strsplit(as.character(formula)[3],"\\+")[[1]] #separa i termini "additivi"..
    idC<-sapply(sapply(lapply(s1, function(x) grep("ps\\(",x)), function(x) (x>=1)), isTRUE)
    stringa<-s1[idC]  #solo i termini con ps
    
    if(any(sapply(stringa, function(.x) grepl("\\* ps\\(", .x)))) stop("invalid usage of symbol '*' in conjunction with ps()")
    if(any(sapply(stringa, function(.x) grepl("\\:ps\\(", .x)))) stop("invalid usage of symbol ':' in conjunction with ps()")
    
    if(any(sapply(stringa, function(.x) grepl("\\):", .x)))) stop("invalid usage of symbol ':' in conjunction with ps()")
    if(any(sapply(stringa, function(.x) grepl("\\) \\*", .x)))) stop("invalid usage of symbol '*' in conjunction with ps()")

    #fun <- function(s) sub("(\\().*(\\))", "\\1\\2", s) 
    #ss1 <- "z:f(5, a=3, b=4, c='1:4', d=2)"
    #ss2 <- "f(5, a=3, b=4, c=\"1:4\", d=2)*z"
    
    #fun(ss1) --> "z:f()"
    #fun(ss2) --> "f()*z" 

    
    if(any(tau>=1 || tau<=0)) stop("tau should be in (0,1)")
    call <- match.call()
    if (missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    intercMt<-attr(mt,"intercept")
    interc<-intercMt==1
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) names(Y) <- nm
      }
    if(!is.null(transf)) {
        Y.orig <- Y
        Y <- eval(parse(text=transf), list(y=Y))
        transf.inv<-splinefun(Y, Y.orig)
      } 
    #----------- vedi dove mettere i seguenti. Va bene qui?
    # colnames(mf)[id.ps]<-testo.ps #devi sostituire i nome altrimenti .getXlevels() non funziona
    # .xlivelli<-.getXlevels(mt, mf) 
    #-----------
    .xlivelli<-.getXlevels(mt, mf) 
    weights <- as.vector(model.weights(mf))
    if(!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
    if(!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
    tf <- terms(formula, specials = c("ps","ridge"))
    id.ps<-attr(tf,"specials")$ps #posizione nel modelframe; pu? essere un vettore se ci sono pi? termini..include y ma non dipende se c'? intercetta o meno..
    testo.ps<-names(mf)[id.ps]
    nomiCoefUNPEN<-names(mf)[-c(1,id.ps)]
    
    #eliminato id.ridge, testo.ridge.. 
    # id.ridge<-attr(tf,"specials")$ridge
    # testo.ridge<-names(mf)[id.ridge]
    #---
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
          else stop("error in the design matrix")#matrix(, NROW(Y), 0L)
    attrContr<-attr(X, "contrasts")
    n<-nrow(X)

    #===========================================================================
    #se NON ci sono termini ps
    #===========================================================================
    if(length(testo.ps)<=0){
      #if(length(testo.ps)<=0 && length(testo.ridge)<=0){
      #fitter<- if(length(tau)>1) get("arq.fitMulti") else get("arq.fit")
      #fit<-fitter(y=Y,X=X,tau=tau,g=g,beta0=b.start,control=control)
      cv<-FALSE #NOn so se serve.. l'ho preso dalla f. del pacchetto..
      fit<-ncross.rq.fitX(y=Y, X=X, taus=tau, eps=eps, sgn.constr=sgn.constr, adjX.constr=adjX.constr) #X gi? include l'intercetta
      if(n.boot>0){
        coef.boot<-array(NA, dim=c(nrow(as.matrix(fit$coef)), length(tau), n.boot))
        id.NA<-NULL
        for(i in 1:n.boot) {
          id<-sample(1:n, size=n, replace=TRUE)
          .o.b<-try(ncross.rq.fitX(y=Y[id], X=X[id,,drop=FALSE], taus=tau, eps=eps, sgn.constr=sgn.constr,adjX.constr=adjX.constr), silent=TRUE)
          if(class(.o.b)!="try-error") coef.boot[,,i]<-.o.b$coef else id.NA[length(id.NA)+1]<-i
        } #end i=1,..,n.boot
        if(!is.null(id.NA)) {
          coef.boot<-coef.boot[,,-id.NA]
          warning(paste(length(id.NA), "boot resamples without convergence"), call.=FALSE)
        }
      } #end if n.boot
      #se NON ci sono termini ps() e' inutile calcolare i gdl con ..
      fit$edf.all<-matrix(1, ncol=length(tau), nrow=ncol(X), dimnames=list(colnames(X),paste(tau)))
      fit$edf.j<- matrix(apply(fit$edf.all, 2,sum), nrow=1)
      #nuovo maggio 2021.. salva cio' che serve per disegnare i termini lineari
      
      X<- X[,-match("(Intercept)", colnames(X),0),drop=FALSE]
      BB<-vector("list", ncol(X))
      names(BB)<-colnames(X)
      for(i in 1:ncol(X)){
        xdisegno<-seq(min(X[,i]), max(X[,i]), l=150)
        BB1<- matrix(xdisegno, ncol=1)
        attr(BB1,"covariate.35")<- xdisegno
        attr(BB1, "coef.names")<- colnames(X)[i]
        attr(BB1,"covariate.n")<- X[,i]
        attr(BB1,"vc")<- FALSE
        attr(BB1,"smoothName") <- colnames(X)[i]
        BB[[i]] <- BB1
      }
      fit$BB<-BB  
    } else {
      #===========================================================================
      #se ci sono termini ps.
      #===========================================================================
      drop.id<-lambda<-S<-B<-BB<-nomiCoefPEN<-nomiVariabPEN<-NULL
      l<-lapply(testo.ps,function(xx)with(mf,eval(parse(text=xx),envir=data)))
      vDeg<-sapply(l,function(xx)attr(xx,"deg"))
      vNdx<-sapply(l,function(xx){if(is.null(attr(xx,"ndx"))) round(min(9,n/4)) else attr(xx,"ndx")})
      vMonot<-sapply(l,function(xx)attr(xx,"monot"))
      vConc<-sapply(l,function(xx)attr(xx,"conc"))
      vDiff<-sapply(l,function(xx)attr(xx,"pdiff"))
      lambda<-unlist(sapply(l,function(xx)attr(xx,"lambda")))
      nomiPS<-unlist(sapply(l,function(xx)attr(xx,"nomeX")))
      var.pen<-unlist(sapply(l,function(xx)attr(xx,"var.pen"))) 
      K<-unlist(sapply(l,function(xx)attr(xx,"K")))
      ridgeList<-unlist(sapply(l,function(xx)attr(xx,"ridge")))
      nomiBy<-unlist(sapply(l,function(xx)attr(xx,"nomeBy")))
      levelsBy<-lapply(l,function(xx)attr(xx,"levelsBy"))
      dimSmooth<-unlist(sapply(l,function(xx)attr(xx,"dimSmooth")))
      decomList<-unlist(sapply(l,function(xx)attr(xx,"decom")))
      knotsList<- lapply(l,function(xx)attr(xx,"nodi"))
      penMatrixList<- lapply(l,function(xx)attr(xx,"penMatrix"))
		  dropcList<- unlist(lapply(l,function(xx)attr(xx,"dropc")))
		  centerList<- unlist(lapply(l,function(xx)attr(xx,"center")))
		  dropvcList<-rep(FALSE, length(dropcList))
		  
		  
		  if(any(centerList & !dropcList)) stop(paste("center=TRUE and dropc=FALSE incompatible for the smooth term #: ", which(centerList & !dropcList)))
		  #per il msg di sotto, devi guardare bene... se ci sono VC terms?
		  #if(any(!dropcList) && "(Intercept)"%in%colnames(X)) stop() 
		  
		  mVariabili<-NULL
		  rangeSmooth<-NULL
		  #se ci sono termini ps()+ps(..,by) il nome delle variabili smooth vengono cambiati per aggiungere la variabile by
		  nomiPS.orig <- nomiPS
		  #browser()
		  nomiPS.By <-paste(nomiPS, nomiBy, sep=":")
		  nomiPS <- unlist(lapply(nomiPS.By, function(.x) sub(":NULL", "", .x)))
		  #nomiPS.By <-paste(nomiBy, nomiPS, sep=":")
		  #nomiPS <- unlist(lapply(nomiPS.By, function(.x) sub("NULL:", "", .x)))
		  
		  ####se la stessa variabile e' specificata come ps o termine lineare..
		  if(length(intersect(nomiPS,nomiCoefUNPEN))>=1) stop("The same variable specified in ps() and as linear term")

		  #ATTENZIONE.. se vuoi subito costruire i nomi ps(x), ps(x):z, ecc...usa i:
		  nomiPS.ps<-sapply(nomiPS.orig, function(x)paste("ps(",x,")",sep=""))
		  nomiPS.ps<-unlist(lapply(paste(nomiPS.ps, nomiBy, sep=":"), function(.x) sub(":NULL", "", .x)))
		  nomiPS.ps.ok<-as.list(nomiPS.ps) #serve lista
		  if(!is.null(unlist(penMatrixList))) stop(" 'pen.matrix' in ps() not yet allowed.")
		  if(any(decomList) || any(ridgeList)) stop(" 'ridge' or 'decompose' in ps() not (yet) implemented.")
		  
		  for(j in id.ps) mVariabili[length(mVariabili)+1]<-mf[j]
      BFixed<-nomiCoefPEN<-B<-BB<-Bderiv<- vector(length=length(mVariabili) , "list")
      smoothVariabili<-byVariabili<-matrix(NA,n,length(mVariabili))
      #
      #
      colMeansBorig.list <- NULL 
      for(j in 1:length(mVariabili)) {
        variabileSmooth<-mVariabili[[j]][,-ncol(mVariabili[[j]]),drop=TRUE]
        if(is.matrix(variabileSmooth)) stop("only univariate smoothing implemented so far")
        smoothVariabili[,j]<-variabileSmooth
        byVariabili[,j]<-mVariabili[[j]][, ncol(mVariabili[[j]]),drop=TRUE]
        nomiBy.j<- if(nomiBy[j]=="NULL") NULL else nomiBy[j]
        if(ridgeList[j]){
          vMonot[j]<-0 #per evitare sorprese..
          vConc[j] <-0 #per evitare sorprese..
          vDiff[j] <-0 #serve per ottenere una penalit? ridge
          vDeg[j]  <-0
          vNdx[j]  <-length(unique(variabileSmooth))
          B[[j]]<- bspline(variabileSmooth,ndx=vNdx[j],deg=vDeg[j])
          if(!is.null(nomiBy.j)) B[[j]]<- byVariabili[,j]*B[[j]]
          BB1<- bspline(unique(variabileSmooth),ndx=vNdx[j],deg=vDeg[j])
          rownames(BB1)<-unique(variabileSmooth)
          attr(BB1,"covariate.n")<- variabileSmooth #NB mVariabili[,j] (che viene assegnato a attr(,"covariate.n")) contiene altri attributi "ndx", "deg", "pdiff", "monot", "lambda","nomeX"
          attr(BB1,"covariate.35")<- "ridge"
          BB[[j]]<-BB1
          Bderiv[[j]]<- "ridge"
        } else {
          rangeSmooth[[j]] <- range(variabileSmooth)
          B[[j]]<- bspline(variabileSmooth,ndx=vNdx[j],deg=vDeg[j], knots=knotsList[[j]]) 
          if(decomList[j]) {
            Dj<-diff(diag(ncol(B[[j]])), diff = vDiff[j])
            B[[j]]<- B[[j]]%*%t(solve(tcrossprod(Dj),Dj)) #sono ncol(B[[j]])-vDiff[j]
            vMonot[j]<-0 #per evitare sorprese..
            BFixed[[j]]<-poly(variabileSmooth, degree=vDiff[j]-1, raw=TRUE)
            colnames(BFixed[[j]])<- paste(nomiPS[j],"fix",1:(vDiff[j]-1),sep=".")
          }
#browser()          
          if(!is.null(nomiBy.j)){
            if(!dropcList[j] && centerList[j]) {
              warning("in vc terms with 'dropc=FALSE', 'center=FALSE' is set.", call.=FALSE)
              centerList[j]<-FALSE
            }
          }
          if(dropcList[j]) B[[j]]<-B[[j]][,-1] #C=contr.treatment(ncol(B))); Bnew<-B%*%C; Pnew<- t(C)%*%Penalty%*%C
          if(is.na(centerList[j])){ 
            centerList[j]<-if(length(tau)>1) FALSE else TRUE
          }
          colMeansBorig.list[[j]] <-colMeans(B[[j]])
          if(centerList[j]) {
            B[[j]]<-sweep(B[[j]], 2, colMeansBorig.list[[j]])
          } else {
            colMeansBorig.list[[j]] <-rep(0, length(colMeansBorig.list[[j]]))
          }
          #-------se VC terms..
          if(!is.null(nomiBy.j)) {
            if(dropcList[j]) dropvcList[j]<-TRUE
            if(is.null(levelsBy[[j]])){
              B[[j]]<- if(dropvcList[j]) cbind(byVariabili[,j], byVariabili[,j]*B[[j]]) else byVariabili[,j]*B[[j]]
              nomiCoefPEN[[j]]<- sapply(1:ncol(B[[j]]), function(x) gsub(":", paste(".",x, ":", sep="") , nomiPS.ps[j]))
              } else {
              M<-model.matrix(~0+factor(byVariabili[,j]))
              B[[j]]<- if(dropvcList[j]) lapply(1:ncol(M), function(.x) cbind(M[,.x],M[,.x]*B[[j]])) else lapply(1:ncol(M), function(.x) M[,.x]*B[[j]])
              nomiCoefPEN[[j]] <- lapply(1:length(B[[j]]), 
                     function(.y)
                        {paste(sapply(1:ncol(B[[j]][[.y]]), 
                                   function(x){gsub(":",  paste(".",x,":",sep=""), nomiPS.ps[j])}), levelsBy[[j]][.y], sep=".")}
                     )
              nomiPS.ps.ok[[j]]<-paste(nomiPS.ps[j], levelsBy[[j]], sep=".")
            }
          } else {
            nomiCoefPEN[[j]]<-paste(nomiPS.ps[j],1:ncol(B[[j]]),sep=".")
          } 
        } #end else 
      } #end for(j in 1:length(mVariabili))
      
      #se B include liste (se ci sono vc terms con by factor), bisogna eliminarle e riportarle nella lista principale B..
      #rep() funziona anche con le liste!!!!
      repl<-pmax(sapply(B,length)*sapply(B,is.list),1)
      colMeansBorig.list <- rep(colMeansBorig.list, repl)
      rangeSmooth<-rep(rangeSmooth, repl)
      nomiBy<-rep(nomiBy, repl)
      knotsList<-rep(knotsList, repl)
      levelsBy<-rep(levelsBy, repl)
      nomiPS.orig <- rep(nomiPS.orig, repl)
      nomiPS.ps <- rep(nomiPS.ps, repl)
      mVariabili<-rep(mVariabili,repl)
      vMonot <- rep(vMonot, repl)
      vDeg <- rep(vDeg, repl)
      vNdx <- rep(vNdx, repl)
      vConc <- rep(vConc, repl)
      vDiff <- rep(vDiff, repl)
      lambda <- rep(lambda, repl)
      dropcList <- rep(dropcList, repl)
      decomList <- rep(decomList, repl)
      dropvcList <- rep(dropvcList, repl)
      centerList <- rep(centerList, repl)
      ridgeList<- rep(ridgeList, repl)
      K<-rep(K, repl)
      id.ps <- rep(id.ps, repl) #il valore restituito non e' indicativo..
      
      #vcList <- sapply(nomiBy, is.null, USE.NAMES =FALSE)
      vcList <- ifelse(nomiBy=="NULL", FALSE, TRUE) 
      
      while(any(sapply(B,is.list))){
        id.vc<-which((sapply(B, is.list)))[1]
        nc<-length(B[[id.vc]])
        B<-append(B, B[[id.vc]], after = id.vc-1)
        for(i in 1:length(B[[id.vc+nc]])) BB<-append(BB, BB[id.vc], id.vc-1)
        B[[id.vc+nc]]<-NULL
        BB[[id.vc+nc]]<-NULL
        nomiCoefPEN <- append(nomiCoefPEN, nomiCoefPEN[[id.vc]], after = id.vc-1)
        nomiCoefPEN[[id.vc+nc]]<-NULL
        #se la lista contiene solo NULL, non funziona...
        penMatrixList <- append(penMatrixList, penMatrixList[[id.vc]], after = id.vc-1)
        penMatrixList[[id.vc+nc]]<-NULL
      }
      #Attenzione se penMatrixList e' NULL, non puo' aumentare di dimensione... quindi aggiungi sotto questo artificio.. 
      if(length(nomiCoefPEN)!=length(penMatrixList)) penMatrixList<- vector(mode = "list", length =length(nomiCoefPEN))
      if(length(knotsList)!=length(vDeg)) knotsList<- vector(mode = "list", length =length(vDeg))
     
      #browser()    
      colMeansB.list <- colMeansBorig.list  
      ########################################################################
      ### un altro ciclo for.. per le colMeans e per le info per i disegni..
      ########################################################################
      for(j in 1:length(B)) {
          if(ridgeList[j]) stop("ridge terms .. non so che fare")
          colmeansB <- colMeansBorig.list[[j]] #sono le medie della B, prima dell'eventuale moltiplic per la variabile VC
          
          #variabileSmooth<-mVariabili[[j]][,-ncol(mVariabili[[j]]),drop=TRUE]
          nomiBy.j<- if(nomiBy[j]=="NULL") NULL else nomiBy[j]
          
          #smoothVariabili[,j]<-variabileSmooth
          #byVariabili[,j]<-mVariabili[[j]][, ncol(mVariabili[[j]]),drop=TRUE]
          #xdisegno<-seq(min(variabileSmooth),max(variabileSmooth), l=150)
          
          xdisegno<-seq(min(rangeSmooth[[j]]), max(rangeSmooth[[j]]), l=150)
          
          BB1<-bspline(xdisegno, ndx=vNdx[j], deg=vDeg[j], knots=knotsList[[j]])
          if(dropcList[j]) BB1<-BB1[,-1]

          #if(!is.null(nomiBy.j)) colmeansB <- colmeansB[-1] #if VC, colmeansB[1] refers to B[[j]][,1] = byVariabili[,j]
          if(centerList[j]) BB1<-sweep(BB1, 2, colMeansBorig.list[[j]])
          if(decomList[j]) {
            BB1<-BB1%*%t(solve(tcrossprod(Dj),Dj))
            attr(BB1,"name.fixed.params")<- colnames(BFixed[[j]])
          }
          attr(BB1,"covariate.n")<- variabileSmooth #NB mVariabili[,j] (che viene assegnato a attr(,"covariate.n")) contiene altri attributi "ndx", "deg", "pdiff", "monot", "lambda","nomeX"
          attr(BB1,"covariate.35")<- xdisegno
          attr(BB1,"colmeansB")<- colmeansB
          attr(BB1,"ndx")<- vNdx[j]
          attr(BB1,"deg")<- vDeg[j]
          attr(BB1,"knots")<- knotsList[[j]]
          attr(BB1,"center")<- centerList[j]
          attr(BB1,"drop")<- dropcList[j]
          #da provare se hai piu VC terms..
          #attr(BB1,"coef.names")<- if(is.list(nomiCoefPEN[[1]])) nomiCoefPEN[[1]][[j]] else nomiCoefPEN[[j]]
          #tolto.. lo fa sotto...
          attr(BB1,"vc")<- if(!is.null(nomiBy.j)) TRUE else FALSE #e' un vc terms?
          attr(BB1,"vcName")<- nomiBy.j
          attr(BB1,"vcLevels")<-levelsBy[[j]]
          attr(BB1,"smoothName")<-as.character(nomiPS.orig[[j]]) #sempre "x", anche nel caso di VC..
          attr(BB1,"smoothName1")<-as.character(nomiPS.ps[[j]]) # "ps(x)", "ps(x):z"
          BB[[j]]<-BB1
          Bderiv[[j]]<-bspline(xdisegno, ndx=vNdx[j], deg=vDeg[j], knots=knotsList[[j]], deriv=1)
          if(dropcList[j]) Bderiv[[j]]<-Bderiv[[j]][,-1]
      } #end for(j in 1:length(mVariabili))
      
    if(length(BB)!=length(nomiCoefPEN)) { #oppure length(B)? dovrebbe essere lo stesso..
      stop("error in dimensions of BB and nomiCoefPEN") 
    } else {
      for(j in 1:length(BB)) attr(BB[[j]], "coef.names") <- nomiCoefPEN[[j]]
    }
    #salvare il singolo livello a cui la base si riferisce ?(in caso di VC)  
    nomiVariabPEN <- unlist(nomiPS.ps.ok)
    #end length(testo.ps)>0
    #    if(length(testo.ridge)>0){.......[SNIP].......}
    if(any(decomList)){
        XfixedB<-do.call(cbind, BFixed)
        X<-cbind(X, XfixedB)
        #interc<-TRUE
    }
    X<- X[,-grep( "ps[(]" , colnames(X)), drop=FALSE]
    #if(!interc && ("(Intercept)"%in%colnames(X))) X<- X[,-match("(Intercept)",colnames(X)), drop=FALSE]
    #if(interc && !("(Intercept)"%in%colnames(X))) X<-cbind("(Intercept)"=1,X)
    if(length(id.ps)>1 && (length(lambda)!=length(id.ps))) stop("with several ps() terms, a single 'lambda' is allowed in each ps()")
    if(length(lambda)>1 && length(id.ps)==1){
      cv<-TRUE
      lambdas<-lambda
      r.cv<- gcrq.rq.cv(Y, B[[1]], X, tau, monotone=vMonot, concave=vConc, ndx=vNdx, lambda=lambda, deg=vDeg, dif=vDiff, 
                   var.pen=var.pen, penMatrix=penMatrixList, lambda.ridge=lambda.ridge, dropcList=dropcList, decomList=decomList, 
                   vcList=vcList, dropvcList=dropvcList, nfolds=nfolds, foldid=foldid, eps=eps, ...)
      lambda<-r.cv[[1]]
    } else {
        cv<-FALSE
    }
    ##ATTENZIONE Ai NOMI.. AGGIUNGER IN FUNZIONE DEI LIVELLI
    #nomiPS.ps
    #[1] "ps(x)"    "ps(x):z1"
    names(lambda)<-names(vMonot)<-names(vConc)<- unlist(nomiPS.ps.ok) #nomiVariabPEN #meglio nomiPS.ps
    #-----------------------------------------------
    #se qualche lambda deve essere stimato..
    #-----------------------------------------------
    if(any(lambda<0)){
        id.lambda.est<-(lambda<0) #which lambdas should be estimated?
        lambda[id.lambda.est]<-lambda0 #.5
        lambdaFixed<-lambda[!id.lambda.est]
        new.rho<-10^3
        K.matrix<-matrix(K, nrow=length(K), ncol=length(tau))
        colnames(K.matrix)<-paste(tau)
        h.ok<-h
        for(j in 1:it.max){
          fit<-suppressWarnings(ncross.rq.fitXB(y=Y, B=B, X=X, taus=tau, monotone=vMonot, concave=vConc, ndx=vNdx,
                  lambda=lambda, deg=vDeg, dif=vDiff, var.pen=var.pen, eps=eps, penMatrix=penMatrixList,
                  lambda.ridge=lambda.ridge, dropcList=dropcList, decomList=decomList, vcList=vcList, dropvcList=dropvcList, 
                  centerList=centerList, colmeansB=colMeansB.list))
          #NB colMeansBorig.list ha un valore in meno per VC terms!
          ##inizio calcolo edf
          if(df.option<=2){
              edf.j<-matrix(sapply(tau, function(xx){edf.rq(fit, xx, vMonot = vMonot, vConc=vConc, output=df.option)}), 
                            ncol=length(tau)) #matrice
          } else {
              #if(penL1.df){
              D.matrix<-fit$D.matrix
				      pesiL1<-abs((D.matrix%*%fit$coefficients)) #matrice dei "pesi". 1colonna=1tau (i primi relativi a termini lineari sono zero..)
              lambdaMatrix<-matrix(lambda, nrow=nrow(as.matrix(lambda)), ncol=length(tau))
              for(k in 1:ncol(pesiL1)){
                  lambdasTutti<- rep(lambdaMatrix[,k], c(vNdx+vDeg-vDiff)) #controlla se serve -1*dropcList
                  if(!is.null(fit$pLin) && fit$pLin>0) lambdasTutti<-c(rep(0, fit$pLin),lambdasTutti)
                  pesiL1[,k]<-sqrt(lambdasTutti/(pesiL1[,k]+.00001))
              }
              edf.j<-matrix(sapply(tau, function(xx){
                  edf.rq(fit, xx, vMonot = vMonot, vConc=vConc, g = g, 
                            pesiL1=as.matrix(pesiL1)[,paste(xx)], output=df.option)
                            }), ncol=length(tau)) #matrice
            }
          colnames(edf.j)<-paste(tau)
          edf.jM<-edf.j
          #==fine calcolo edf
          if(single.lambda){ #un unico lambda per tutti i percentili..
            dev2u<-apply(fit$all.dev2u,1,sum) #all.dev2u e' matrice . con un singolo tau era "fit$dev2u"
            if(attr(fit$id.coef, "nomi")[1]!="Xlin") edf.jM<-rbind(0, edf.jM) #edf.j<-c(0,edf.j)
            edf.j<-apply(edf.jM,1,sum)
            sigma2u <-psi.k<- dev2u/(K*edf.j[-1])  #e' un vettore
            denom.e<-if(method=="REML") max(1,n*length(tau)-sum(rbind(1,K.matrix)*edf.jM)) else (n*length(tau))
            sigma2e<-phi<-sum(fit$all.dev2e)/denom.e  #sopra, se non vuoi pesare usa "sum(edf.j)"  
          } else {
            if(attr(fit$id.coef, "nomi")[1]!="Xlin") edf.j<-edf.jM<-rbind(0,edf.j)   
            sigma2u <-psi.k.tau<-fit$all.dev2u/(K.matrix*edf.j[-1,]) ##e' una matrice
            denom.e<-if(method=="REML") (n-apply(edf.j,2,sum)) else n
            sigma2e<-phi.tau<-fit$all.dev2e/denom.e 
            sigma2e<-matrix(sigma2e, ncol=length(tau), nrow=nrow(sigma2u), byrow = TRUE)
          }
          
          #browser()
          #fit$all.dev2e #corrisponde a fit$rho, vettore (o scalare) delle loss
          #apply(edf.j,2,sum) #un vettore.. sono i edf (NON pesati per il fattore K) per ogni curva Q(tau).. comprendono anche i df parametrici
           # 
          #if(laplace) {
          #  est.scale<- if(method=="ML") fit$rho/n else fit$rho/(n-sum(edf.j)) #n  #n-sum(edf.j) #(se REML)
          #  sigma2e<- varAL(est.scale, tau) #? vettore se tau ? vettore
          #  sigma2e<-sum(sigma2e*tau/sum(tau))
          #} else {
          #  if(all.perc && length(tau)>1){
          #      sigma2e<- if(err.rho) mean(fit$all.dev2e)/((n-sum(edf.j*c(1,K)))^2) else mean(fit$all.dev2e)/(n-sum(edf.j*c(1,K))) #oppure: sum(fit$all.dev2e)/sum(n-sum(edf.j)) ?
          #  } else {
          #      den<-if(method=="ML") n else (n-sum(edf.j*c(1,K)))  
          #      sigma2e<- if(err.rho) fit$dev2e/(den^2) else fit$dev2e/den 
          #    }
          #}
          lambda.old<-lambda
          ###########################################################
          lambda<- sigma2e/(sigma2u+.00001*sigma2e) #scalare, vettore o MATRICE!!! CONTROLLA!!! naturalemnete vedi anche fitXB.r()
          ###########################################################
          #attenzione ai nomi... magari names(lambda)<-names(lambda.old) ??
          
          if(!single.lambda) lambda.old<-matrix(lambda.old, ncol=ncol(lambda), nrow=nrow(lambda))
          
          d<-lambda-lambda.old
          if(j>=it.max) h.ok<-h^(j-it.max)
          lambda<-pmin(lambda.old+h.ok*d,lambda.max)
          lambda[!id.lambda.est]<-lambdaFixed #funziona anche con matrici !!! :-))) che culo!!! :-)))
          #---------------------------------------------
          #se vuoi valutare diff in rho..
          old.rho<-new.rho
          new.rho<-fit$rho
          #---------------------------------------------
          if(all(abs((lambda.old-lambda)/(lambda.old+1e-8))< tol)) {
            break
          }
          #if() Fellner: change in Variance/sum(AllVariances)
          ###### non convergenza: continua!
          if(display) {
            
            flush.console()
            cat("iteration:", if(j<10) paste("",j) else j,
                " lambda:", formatC(ifelse(lambda==lambda.max,Inf,lambda), digits=3, width=8, format="f"), "\n") #prima era format="f"/"g"
            #writeLines(strwrap(paste("iteration", j, " lambda", formatC(ifelse(lambda==lambda.max,Inf,lambda), digits=3, width=7, format="f")),"\n"))
            }
          #---------------------------------------------
          #if(all(abs((old.rho-new.rho)/(old.rho+1e-8))< tol)) break
          #if(all(abs((lambda.old-lambda)/(lambda.old+1e-8))< tol)) break
          #if() Fellner: change in Variance/sum(AllVariances)
        } #end for(j in 1.. it.max)
        lambda<-lambda.old
        names(edf.j)<-rownames(edf.jM)<-c("Xlin",nomiVariabPEN)        
        if(attr(fit$id.coef, "nomi")[1]!="Xlin") {
          edf.j<-edf.j[-1] #rimuovi il primo lo zero eventualmente aggiunto prima.. #e' importante altrimenti i nomi non corrispondono..
          edf.jM<-edf.jM[-1,,drop=FALSE]
        }
        #----------------------------------------------------------
        #se tutti i lambda sono *assegnati*
        #----------------------------------------------------------
      } else {  #stima il modello *assegnati* lambda
        fit<-ncross.rq.fitXB(y=Y, B=B, X=X, taus=tau, monotone=vMonot, concave=vConc, ndx=vNdx, 
          lambda=lambda, deg=vDeg, dif=vDiff, var.pen=var.pen, eps=eps, penMatrix=penMatrixList, 
          lambda.ridge=lambda.ridge,dropcList=dropcList,decomList=decomList, vcList=vcList, dropvcList=dropvcList, 
          centerList=centerList, colmeansB=colMeansB.list)
        
        #colnames(edf.j)<-paste(tau)
        #rownames(edf.j)<- if(nrow(edf.j)==length(nomiVariabPEN)) nomiVariabPEN else c("Xlin",nomiVariabPEN)
        #edf.jM<-edf.j
        if(df.option<=2){
              edf.j<-matrix(sapply(tau, function(xx){edf.rq(fit, xx, vMonot = vMonot, vConc=vConc, output=df.option)}), ncol=length(tau)) #matrice
        } else {
              #if(penL1.df){
              D.matrix<-fit$D.matrix
				      pesiL1<-abs((D.matrix%*%fit$coefficients)) #matrice dei "pesi". 1colonna=1tau (i primi relativi a termini lineari sono zero..)
              lambda<-matrix(lambda, nrow=nrow(as.matrix(lambda)), ncol=length(tau))
              for(k in 1:ncol(pesiL1)){
                lambdasTutti<- rep(lambda[,k], c(vNdx+vDeg-vDiff)) #controlla se serve -1*dropcList
                if(!is.null(fit$pLin) && fit$pLin>0) lambdasTutti<-c(rep(0, fit$pLin),lambdasTutti)
                pesiL1[,k]<-sqrt(lambdasTutti/(pesiL1[,k]+.00001))
              }
              edf.j<-edf.jM<-matrix(sapply(tau, function(xx){
                  edf.rq(fit, xx, vMonot = vMonot, vConc=vConc, g = g, 
                        pesiL1=as.matrix(pesiL1)[,paste(xx)], output=df.option)}), ncol=length(tau)) #matrice
        }
        #browser()
        #nomiPS.ps.ok oppure nomiVariabPEN?
        colnames(edf.j)<-colnames(edf.jM)<-paste(tau)
        rownames(edf.j)<- if(nrow(edf.j)==length(nomiVariabPEN)) nomiVariabPEN else c("Xlin",nomiVariabPEN)
        #==fine calcolo edf
        edf.jM<-edf.j
      } #end di stima il modello *assegnati* lambda
      ########################################################
    
      edf.all<- matrix(sapply(tau, function(xx){
          edf.rq(fit, xx, return.all=TRUE, output=df.option, vMonot=vMonot, vConc=vConc, g=g, 
                  pesiL1=as.matrix(pesiL1)[,paste(xx)])}), ncol=length(tau))
      colnames(edf.all)<-paste(tau)
      #-----------------------------
      #attenzione a questo punto edf.j dovrebbe gia' esserci.. METTERLO comunque? 
      #da verificare se c'e un termine smooth e se c'e' o no l'intercetta..
      #con lambda fissati, tutto ok.. non ci vuole questo pezzo di codice..
      #---------------------------
      if(n.boot>0){
        B.orig<-B
        coef.boot<-array(NA, dim=c(nrow(as.matrix(fit$coef)), length(tau), n.boot))
        id.NA<-NULL
        for(i in 1:n.boot) {
          id<-sample(1:n, size=n, replace=TRUE)
          #la linea di sotto non va bene perche' quando si passano a ncross.rq.fitXB() gli argomenti  x, nomiBy e byVariabili (invece che B)
          #poi nella costruzione della base con vc elimina una colonna e quindi un parametro non viene stimato..
          #
          #.o.b<- try(ncross.rq.fitXB(y=Y[id], B=NULL, x=smoothVariabili[id,,drop=FALSE], nomiBy=nomiBy, 
          #                byVariabili=byVariabili[id,,drop=FALSE], X=X[id,, drop=FALSE], 
          #                taus=tau, monotone=vMonot, concave=vConc, ndx=vNdx, lambda=lambda, deg=vDeg, dif=vDiff, 
          #                var.pen=var.pen, eps=eps, penMatrix=penMatrixList,lambda.ridge=lambda.ridge,
          #                dropcList=dropcList,decomList=decomList,dropvcList=dropvcList)
          #                , silent=TRUE)
          
          #suppressWarnings?
          for(j in 1:length(B.orig)) B[[j]]<-B.orig[[j]][id,,drop=FALSE]
          .o.b<- try(ncross.rq.fitXB(y=Y[id], B=B, X=X[id,, drop=FALSE], taus=tau, monotone=vMonot, concave=vConc, ndx=vNdx,
                            lambda=lambda, deg=vDeg, dif=vDiff, var.pen=var.pen, eps=eps, penMatrix=penMatrixList,
                            lambda.ridge=lambda.ridge, dropcList=dropcList, decomList=decomList, vcList=vcList, dropvcList=dropvcList,
                            centerList=centerList, colmeansB=colMeansB.list),
                                silent=TRUE)

          if(class(.o.b)!="try-error") coef.boot[,,i]<-.o.b$coef else id.NA[length(id.NA)+1]<-i
        } #end i=1,..,n.boot
        if(!is.null(id.NA)) {
          coef.boot<-coef.boot[,,-id.NA]
          warning(paste(length(id.NA), "boot resamples without convergence"), call.=FALSE)
        }
      } #end if n.boot
      nomiCoefUNPEN<-rownames(as.matrix(fit$coefficients))[1:ncol(X)]
      nn<-c(nomiCoefUNPEN,unlist(nomiCoefPEN))
      #if(ncol(X)>0) nn<-c("(Intercept)",nn)
      if(is.matrix(fit$coefficients)) rownames(fit$coefficients)<-nn else names(fit$coefficients)<-nn
      #if(length(tau)>1) rownames(fit$coefficients)<-nn else names(fit$coefficients)<-nn
      
      names(BB)<-nomiVariabPEN
      
      
      #per i nomi? names(colMeansB.list) <-
      names(dropvcList)<-names(centerList)<-names(decomList)<-names(dropcList)<- names(vNdx)<- names(vDeg)<- names(vDiff)<- nomiVariabPEN
	    fit$info.smooth<-list(monotone=vMonot, ndx=vNdx, lambda=lambda, deg=vDeg, dif=vDiff, concave=vConc, 
                    dropcList=dropcList, decomList=decomList, centerList=centerList, vcList=vcList, dropvcList=dropvcList, testo.ps=nomiVariabPEN) #che e' nomiPS.ps.ok. Prima era testo.ps

	    ###########################
      ######### SERVE??????????????
      #if(is.matrix(edf.j)) rownames(edf.j)<-attr(fit$id.coef, "nomi") else names(edf.j)<-attr(fit$id.coef, "nomi")
      if(is.matrix(lambda)) {
        rownames(lambda)<-nomiVariabPEN
        colnames(lambda)<-paste(tau)
      } else {
        names(lambda)<-nomiVariabPEN
      }
      #########
      if(df.option<=2){
        edf.j<-as.matrix(edf.j)
        colnames(edf.j)<-paste(tau)
        fit$edf.j<-if(!interc) rbind(1,edf.j)  else edf.j
      } else {
        #if(df.option>=3){ #salva edf.all solo se df.option=3; infatti edf.all non ? definito con df.option<=2
        if(is.matrix(edf.all)) rownames(edf.all)<-nn else names(edf.all)<-nn
        fit$edf.all<-edf.all
        fit$edf.j<-edf.jM
      }
      #names(BB)<-names(fit$lambda)<-names(fit$edf.smooth)<-nomiVariabPEN
      #fit$B<-B  #n righe..
      fit$BB<-BB #35 (in realt? 100) righe e attr "covariate.n" "covariate.35". #NB attr(,"covariate.n") contiene altri attributi "ndx", "deg", "pdiff", "monot", "lambda","nomeX"
      fit$Bderiv<-Bderiv
      fit$lambda<-lambda
    } #fine del "se ci sono termini smooth"
    #========================================
    #========================================
    if(y) fit$y<-Y
    if(!is.null(transf)) {
      #if(!is.null(transf)) fit$transf.inv<-transf.inv
      attr(fit$fitted.values, "transf.inv") <- transf.inv
      attr(fit$fitted.values, "transf") <- eval(parse(text=paste("function(y){", transf, "}"))) 
    }
    fit$contrasts <- attrContr  
    fit$xlevels <- .xlivelli 
    fit$taus<-tau
    fit$call<-call
    if(n.boot){ 
      fit$boot.coef<- coef.boot
    }
    if(cv) {
      fit$cv <- cbind(lambdas,r.cv[[2]])
      colnames(fit$cv)<-c("lambdas",paste("fold",1:ncol(r.cv[[2]]), sep=""))
      #fit$foldid<-attr(r.cv[[2]], "foldid")
      attr(fit$cv, "foldid")<-attr(r.cv[[2]], "foldid")
      #warning se lambda selezionato e' al limite..
      valoriL<-fit$cv[,1]
      valoriCV<-apply(fit$cv[,-1],1,mean)
      lambda.ott<- valoriL[which.min(valoriCV)]
      boundary<-FALSE
      if(lambda.ott<=min(valoriL)) {
        boundary<-TRUE
        etic<- "left"
      }
      if(lambda.ott>=max(valoriL)) {
        boundary<-TRUE
        etic<- "right"
      }
      if(boundary) warning("The lambda selected by CV is on the", paste(" (",etic,") ",sep=""),"boundary of the range set in ps(..,lambda=..)", call. = FALSE)
    }
    class(fit)<-c("gcrq")
    fit
  }

