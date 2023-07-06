gcrq <-
function(formula, tau=c(.1,.25,.5,.75,.9), data, subset, weights, na.action, transf=NULL,
    y=TRUE, n.boot=0, eps=0.001,  #ho rimosso cv=FALSE, interc=TRUE,
    display=FALSE, method=c("REML","ML"), df.opt=2, df.nc=FALSE,
    lambda0=.1, h=0.8, lambda.max=1000, tol=0.01, it.max=15, single.lambda=TRUE, 
    foldid=NULL, nfolds=10, lambda.ridge=0, adjX.constr=TRUE, #sgn.constr=NULL, ,
    contrasts=NULL, sparse=FALSE){
#df.nc: se TRUE, i df calcolati tengono conto dei vincoli noncross (ignorato se la stima e' di una singolo tau)
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
#-----------------------
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
    B <- splineDesign(knots, x, ord = deg + 1, derivs = rep(deriv, length(x)), outer.ok=outer.ok) #,sparse=TRUE
    B
}
  ##### ====================================================================================
  ##### INIZIO FUNZIONE
    #settiamo qualche opzione che forse non ? molto utile..
    #sparse=FALSE
    n.points=200 #questo serve per i disegni e per i vincoli di monot e conc (se ps(..,constr.fit=TRUE))
    
    #all.perc=FALSE
    #g=0.20 #usato solo se df.opt=3 (comunque sconsigliato)
    
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
    
    if(length(tau)<=1) df.nc<-FALSE
    if(any(tau>=1 | tau<=0)) stop("tau should be in (0,1)")
    tau<-sort(tau)
    call <- match.call()
    if (missing(data)) data <- environment(formula)
    tf <- terms(formula, specials = "ps")
    id.ps<-attr(tf,"specials")$ps #posizione nel modelframe; vettore se ci sono piu' termini..include y ma non da interc

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
   
    #browser()
    
    mf <- eval(mf, parent.frame())
    
    mt <- attr(mf, "terms")
    intercMt<-attr(mt,"intercept")
    interc<-intercMt==1
    Y <- stats::model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) names(Y) <- nm
      }
    if(!is.null(transf)) {
        Y.orig <- Y
        Y <- eval(parse(text=transf), list(y=Y))
        transf.inv<-splinefun(Y, Y.orig, ties=min, method="monoH.FC")
    }
    
    
    #browser()
    
    #----------- vedi dove mettere i seguenti. Va bene qui?
    # colnames(mf)[id.ps]<-testo.ps #devi sostituire i nome altrimenti .getXlevels() non funziona
    # .xlivelli<-.getXlevels(mt, mf) 
    #-----------
    .xlivelli<-.getXlevels(mt, mf) 
    weights <- as.vector(model.weights(mf))
    if(!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
    if(!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
    testo.ps<-names(mf)[id.ps]
    nomiCoefUNPEN<-names(mf)[-c(1,id.ps)]
    
    #eliminato id.ridge, testo.ridge.. 
    # id.ridge<-attr(tf,"specials")$ridge
    # testo.ridge<-names(mf)[id.ridge]
    #---

    X <- if(!is.empty.model(mt)){
      stats::model.matrix(mt, mf, contrasts) 
      } else {
        stop("error in the design matrix")
        }#matrix(, NROW(Y), 0L)
    attrContr<-attr(X, "contrasts")
    n<-nrow(X)
    #browser()
    
    if(ncol(X)>0 && any(apply(X,2,sd)==0)){
      colnames(X)[which(apply(X,2,sd)==0)] <- "(Intercept)"
    }
    #===========================================================================
    #se NON ci sono termini ps
    #===========================================================================
    if(length(testo.ps)<=0){
      #if(length(testo.ps)<=0 && length(testo.ridge)<=0){
      #fitter<- if(length(tau)>1) get("arq.fitMulti") else get("arq.fit")
      #fit<-fitter(y=Y,X=X,tau=tau,g=g,beta0=b.start,control=control)
      cv<-FALSE #NOn so se serve.. l'ho preso dalla f. del pacchetto..
      fit <-ncross.rq.fitX(y=Y, X=X, taus=tau, eps=eps, adjX.constr=adjX.constr) #X gi? include l'intercetta
      if(n.boot>0){
        coef.boot<-array(NA, dim=c(nrow(as.matrix(fit$coef)), length(tau), n.boot))
        id.NA<-NULL
        for(i in 1:n.boot) {
          id<-sample(1:n, size=n, replace=TRUE)
          .o.b<-try(ncross.rq.fitX(y=Y[id], X=X[id,,drop=FALSE], taus=tau, eps=eps, adjX.constr=adjX.constr), silent=TRUE)
          #if(class(.o.b)!="try-error") 
          if(!inherits(.o.b,"try-error"))  coef.boot[,,i]<-.o.b$coef else id.NA[length(id.NA)+1]<-i
        } #end i=1,..,n.boot
        if(!is.null(id.NA)) {
          coef.boot<-coef.boot[,,-id.NA]
          warning(paste(length(id.NA), "boot resamples without convergence"), call.=FALSE)
        }
      } #end if n.boot
      #se NON ci sono termini ps() e' inutile calcolare i gdl con ..
      fit$edf.all<-matrix(1, ncol=length(tau), nrow=ncol(X), dimnames=list(colnames(X),paste(tau)))
      if(df.nc) {
        DF.all<-cbind(fit$DF.NEG[,ncol(fit$DF.NEG):1,drop=FALSE], 0, fit$DF.POS)
        fit$edf.all <- fit$edf.all - DF.all
      }
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
      K<- sapply(l,function(xx)attr(xx,"K"))
      ridgeList<-unlist(sapply(l,function(xx)attr(xx,"ridge"))) #TRUE se la variabile e' una matrice o un fattore
      nomiBy<-unlist(sapply(l,function(xx)attr(xx,"nomeBy")))
      levelsBy<-lapply(l,function(xx)attr(xx,"levelsBy"))
      dimSmooth<-unlist(sapply(l,function(xx)attr(xx,"dimSmooth")))
      decomList<-unlist(sapply(l,function(xx)attr(xx,"decom")))
      knotsList<- lapply(l,function(xx)attr(xx,"nodi"))
      penMatrixList <- lapply(l,function(xx)attr(xx,"penMatrix"))
		  dropcList<- unlist(lapply(l,function(xx)attr(xx,"dropc")))
		  centerList<- unlist(lapply(l,function(xx)attr(xx,"center")))
		  constr.fitList<- unlist(lapply(l,function(xx)attr(xx,"constr.fit")))
		  dropvcList<-rep(FALSE, length(dropcList))
		  shared.penList<- unlist(lapply(l,function(xx)attr(xx,"shared.pen")))
		  adList<- unlist(sapply(l,function(xx)attr(xx,"adapt")))
		  scList<- unlist(sapply(l,function(xx)attr(xx,"sc")))
		  
		  #controlla se ci soo vincoli monot e concav con constr.fit=FALSE e base centrata
		  #vMonot, vConc, constr.fitList,dropcList,centerList
		  
		  #browser()
		  if(any(centerList & !dropcList)) stop(paste("center=TRUE and dropc=FALSE incompatible for the smooth term #: ", which(centerList & !dropcList)))
		  #per il msg di sotto, devi guardare bene... se ci sono VC terms?
		  #if(any(!dropcList) && "(Intercept)"%in%colnames(X)) stop() 
		  
		  if(any(shared.penList) && !single.lambda) stop("single.lambda=FALSE is not (yet) allowed with shared penalty")
		  
		  
		  if(any(decomList)) stop("spline decomposition not yet allowed")
		  
		  if(any(decomList) && length(tau)>1) stop("spline decomposition is not allowed with multiple quantile")
		  
		  mVariabili<-NULL
		  rangeSmooth<-NULL
		  #se ci sono termini ps()+ps(..,by) il nome delle variabili smooth vengono cambiati per aggiungere la variabile by
		  nomiPS.orig <- nomiPS
		  
		  nomiPS.By <- paste(nomiPS, nomiBy, sep=":")
		  nomiPS <- unlist(lapply(nomiPS.By, function(.x) sub(":NULL", "", .x)))
		  #nomiPS.By <-paste(nomiBy, nomiPS, sep=":")
		  #nomiPS <- unlist(lapply(nomiPS.By, function(.x) sub("NULL:", "", .x)))
		  
		  ####se la stessa variabile e' specificata come ps o termine lineare..
		  if(length(intersect(nomiPS,nomiCoefUNPEN))>=1) stop("The same variable specified in ps() and as linear term")

		  #ATTENZIONE.. se vuoi subito costruire i nomi ps(x), ps(x):z, ecc...usa i:
		  nomiPS.ps<- sapply(nomiPS.orig, function(.x)paste("ps(",list(.x),")",sep=""))
		  nomiPS.ps<-unlist(lapply(paste(nomiPS.ps, nomiBy, sep=":"), function(.x) sub(":NULL", "", .x)))
		  nomiPS.ps.ok<-as.list(nomiPS.ps) #serve lista

		  for(j in id.ps) mVariabili[length(mVariabili)+1]<-mf[j]
      BFixed<-nomiCoefPEN<-B<-BB<-Bderiv<- vector(length=length(mVariabili) , "list")
      #
      #
      ps.matrix.list<-colMeansBorig.list <- NULL 
      #ps.matrix.list include informazioni se il termine ps si riferisce ad un fattore (per random intercepts)
      
      #browser()
      
      for(j in 1:length(mVariabili)) {
        if(nomiBy[j]=="NULL"){
          nomiBy.j <- NULL
          #variabileSmooth<- if(is.matrix(mVariabili[[j]]))  else c(mVariabili[[j]]) #c() rimuove gli attributi! Prima era drop()
          variabileSmooth<- drop(mVariabili[[j]]) #c() non va bene con matrici! le converte in vettori!
        } else { #se ci sono termini VC
          nomiBy.j <- nomiBy[j]
          variabileSmooth<- mVariabili[[j]][,-ncol(mVariabili[[j]]),drop=TRUE] 
          variabileBy<- mVariabili[[j]][, ncol(mVariabili[[j]]),drop=TRUE]
          #variabileSmooth<- c(mVariabili[[j]][,-ncol(mVariabili[[j]]),drop=TRUE]) 
          #variabileBy<- c(mVariabili[[j]][, ncol(mVariabili[[j]]),drop=TRUE])
        } 
        #browser()
        
        #if(adList[j]>0) stop(" adaptive smoothing not (yet) implemented")
        if(adList[j]<0 || adList[j]>1 ) stop(" 'ad' should be in [0,1]")
        
        for(jj in c("center", "dropc", "shared.pen", "decom", "dimSmooth", "nomeBy", 
          "ridge", "K", "nomeX", "lambda", "constr.fit", "conc", "monot", 
          "pdiff", "deg","adapt","sc")) attr(variabileSmooth,jj)<-NULL
        
        #browser()
        if(sparse) decomList[j]<-centerList[j] <- dropcList[j]<- ridgeList[j] <- constr.fitList[j]<- FALSE

        if(ridgeList[j]){
          dropcList[j] <- FALSE
          vMonot[j]<-0 #per evitare sorprese..
          vConc[j] <-0 #per evitare sorprese..
          vDiff[j] <-0 #serve per ottenere una penalita' ridge
          vDeg[j]  <-0
          vNdx[j]  <-length(unique(variabileSmooth))
          centerList[j] <- FALSE
          #browser()
          if(is.matrix(variabileSmooth)) {
            #if(length(tau)>1) stop("Variable selection is not allowed with multiple quantile curves")
            ps.matrix.list[j] <- TRUE
            B[[j]]<-variabileSmooth
            if(scList[j]) B[[j]] <- scale(variabileSmooth)
            colnames(B[[j]])<- if(is.null(colnames(variabileSmooth))) paste(nomiPS[j], 1:ncol(variabileSmooth), sep="") else colnames(variabileSmooth)
            if(is.null(K[[j]])) K[[j]] <- log(nrow(B[[j]])/(ncol(B[[j]])^(2/3))) #log(n/p^(1/2)) e' troppo severo! (lambda tende ad essere grande)
            } else { #se random intercepts
              ps.matrix.list[j] <- FALSE
              B[[j]]<-stats::model.matrix(~0+factor(variabileSmooth))
              colnames(B[[j]])<- paste(nomiPS[j], sort(unique(variabileSmooth)), sep="")
              if(is.null(K[[j]])) K[[j]] <- log(nrow(B[[j]])/(ncol(B[[j]])^(2/3))) #log(n/p^(1/2))
            }
          #browser()
          #se NON e' specificata, la matrice di Pen e' Ident
          if(is.null(penMatrixList[[j]])) penMatrixList[[j]] <- diag(ncol(B[[j]]))
          
          if(!is.null(nomiBy.j)) B[[j]]<- variabileBy*B[[j]]
          colMeansBorig.list[[j]]<-rep(0, ncol(B[[j]]))
          nomiCoefPEN[[j]] <- colnames(B[[j]])
          rangeSmooth[[j]] <- range(variabileSmooth)
        } else {
          if(is.null(K[[j]])) K[[j]] <- 2
          ps.matrix.list[j] <- FALSE
          rangeSmooth[[j]] <- range(variabileSmooth)
          B[[j]]<- bspline(variabileSmooth,ndx=vNdx[j],deg=vDeg[j], knots=knotsList[[j]])

          if(decomList[j]) {
            Dj<-diff(diag(ncol(B[[j]])), diff = vDiff[j])
            B[[j]]<- B[[j]]%*%t(solve(tcrossprod(Dj),Dj)) #sono ncol(B[[j]])-vDiff[j]
            BFixed[[j]]<-poly(variabileSmooth, degree=vDiff[j]-1, raw=TRUE)
            colnames(BFixed[[j]])<- paste(nomiPS[j],"fix",1:(vDiff[j]-1),sep=".")
            vDiff[j]<-0
            vMonot[j]<-0 #per evitare sorprese.. Pero' se uno vuole imporre la monot? Dovrebbe essere possibile...
          }
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
            B[[j]]<-sweep(B[[j]], 2, colMeansBorig.list[[j]]) ## subtract the column means
          } else {
            colMeansBorig.list[[j]] <-rep(0, length(colMeansBorig.list[[j]]))
          }
          #-------se VC terms..
          
          if(!is.null(nomiBy.j)) {
            if(dropcList[j]) dropvcList[j]<-TRUE
            if(is.null(levelsBy[[j]])){ #se e' vc con variabile continua
              #CAMBIO 25/10 non aggiungere MAI variabileBy alla base..
              #B[[j]]<- if(dropvcList[j]) cbind(variabileBy, variabileBy*B[[j]]) else variabileBy*B[[j]]
              B[[j]] <-variabileBy*B[[j]]
              nomiCoefPEN[[j]]<- sapply(1:ncol(B[[j]]), function(x) gsub(":", paste(".",x, ":", sep="") , nomiPS.ps[j]))
            } else {
              names(levelsBy)[j] <- nomiBy[j]
              M<-stats::model.matrix(~0+factor(variabileBy))
              #B[[j]]<- if(dropvcList[j]) lapply(1:ncol(M), function(.x) cbind(M[,.x],M[,.x]*B[[j]])) else lapply(1:ncol(M), function(.x) M[,.x]*B[[j]])
              B[[j]]<- lapply(1:ncol(M), function(.x) M[,.x]*B[[j]])
              nomiCoefPEN[[j]] <- lapply(1:length(B[[j]]), 
                     function(.y)
                        {paste(sapply(1:ncol(B[[j]][[.y]]), 
                                   function(x){gsub(":",  paste(".",x,":",sep=""), nomiPS.ps[j])}), levelsBy[[j]][.y], sep=".")}
                     )
              nomiPS.ps.ok[[j]]<-paste(nomiPS.ps[j], levelsBy[[j]], sep=".")
            }
          } else {
            nomiCoefPEN[[j]]<- paste(nomiPS.ps[j],1:ncol(B[[j]]),sep=".")
          } 
        } #end else se non e' ridge..
      } #end for(j in 1:length(mVariabili))
      
      nomiCateg <-unlist(sapply(levelsBy, function(.x) if(is.null(.x)) NA else .x))
      
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
      constr.fitList<-rep(constr.fitList, repl)
      shared.penList<-rep(shared.penList, repl)
      K<-rep(K, repl)
      adList <-rep(adList, repl)
      scList <-rep(scList, repl)
      id.ps <- rep(id.ps, repl) #il valore restituito non e' indicativo..
      ps.matrix.list <- rep(ps.matrix.list, repl)
      #browser()
      
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
     
      colMeansB.list <- colMeansBorig.list  
      ########################################################################
      ### un altro ciclo for.. per le colMeans e per le info per i disegni..
      ########################################################################
      #browser()
      for(j in 1:length(B)) {
          nomiBy.j<- if(nomiBy[j]=="NULL") NULL else nomiBy[j]
          variabileSmooth <- if(is.null(nomiBy.j)) drop(mVariabili[[j]]) else mVariabili[[j]][,-ncol(mVariabili[[j]]),drop=TRUE]
          if(ridgeList[j]){
            #variabileSmooth <- if(is.null(nomiBy.j)) drop(mVariabili[[j]]) else mVariabili[[j]][,-ncol(mVariabili[[j]]),drop=TRUE]
            BB1 <- if(is.matrix(variabileSmooth)) variabileSmooth else stats::model.matrix(~0+ factor(unique(variabileSmooth)))
            colnames(BB1)<- colnames(B[[j]])
            attr(BB1,"covariate.n")<- variabileSmooth #NB mVariabili[,j] (che viene assegnato a attr(,"covariate.n")) contiene altri attributi "ndx", "deg", "pdiff", "monot", "lambda","nomeX"
            attr(BB1,"covariate.35")<- "ridge"
            attr(BB1,"vc")<- if(!is.null(nomiBy.j)) TRUE else FALSE #e' un vc terms?
            attr(BB1,"vcName")<- nomiBy.j
            #attr(BB1,"vcLevels")<-levelsBy[[j]]
            attr(BB1,"smoothName")<-  deparse(nomiPS.orig[[j]]) #sempre "x", anche nel caso di VC.. deparse() o as.character()?
            attr(BB1,"smoothName1")<- deparse(nomiPS.ps[[j]]) # "ps(x)", "ps(x):z"
            BB[[j]]<-BB1
            Bderiv[[j]]<- "ridge"
          } else { #per termini smooth usuali..
            colmeansB <- colMeansBorig.list[[j]] #sono le medie della B, prima dell'eventuale moltiplic per la variabile VC
            xdisegno <- seq(min(rangeSmooth[[j]]), max(rangeSmooth[[j]]), l= n.points)
            BB1<-bspline(xdisegno, ndx=vNdx[j], deg=vDeg[j], knots=knotsList[[j]])
            if(dropcList[j]) BB1<-BB1[,-1]
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
            attr(BB1,"diff")<- vDiff[j]
            attr(BB1,"knots")<- knotsList[[j]]
            attr(BB1,"center")<- centerList[j]
            attr(BB1,"drop") <- dropcList[j]
            attr(BB1,"vc") <- if(!is.null(nomiBy.j)) TRUE else FALSE #e' un vc terms?
            attr(BB1,"vcName") <- nomiBy.j
            attr(BB1,"vcLevels") <- levelsBy[[j]]
            attr(BB1,"vcCategory") <- if(is.na(nomiCateg[j])) NULL else nomiCateg[j] 
            #dopo la 1.2-1 alla versione 1.3.0 ho eliminato as.character() e ho messo deparse()perche' se c'era log(x) lo separava in "log" "x" 
            #     e plot.gcrq() che usa "smoothName" per xlab metteva l'xlab non corretto..
            #attr(BB1,"smoothName")<-as.character(nomiPS.orig[[j]]) #sempre "x", anche nel caso di VC..
            #attr(BB1,"smoothName1")<-as.character(nomiPS.ps[[j]]) # "ps(x)", "ps(x):z"
            
            attr(BB1,"smoothName") <- deparse(nomiPS.orig[[j]]) #sempre "x", anche nel caso di VC.. metti nomiPS se vuoi i nomi con le interaz g:x
            attr(BB1,"smoothName1") <- deparse(nomiPS.ps[[j]]) # "ps(x)", "ps(x):z"
            attr(BB1,"constr.fit") <- constr.fitList[j]
            BB[[j]]<-BB1
            Bderiv[[j]]<-bspline(xdisegno, ndx=vNdx[j], deg=vDeg[j], knots=knotsList[[j]], deriv=1)
            if(dropcList[j]) Bderiv[[j]]<-Bderiv[[j]][,-1]
            
            if(is.null(penMatrixList[[j]])) {
              penMatrixList[[j]] <- if(dropcList[j]) diff(diag(ncol(BB[[j]])+1), diff=vDiff[j])[,-1]
                                          else diff(diag(ncol(BB[[j]])), diff=vDiff[j])
              
              #penMatrixList[[j]] <- diff(diag(vNdx[j]+vDeg[j]), diff=vDiff[j])
              #if(dropcList[j]) penMatrixList[[j]] <- penMatrixList[[j]][,-1]
              
            }
          }
      } #end for(j in 1:length(mVariabili))
      
    #browser()
    
    
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
    
   if(sparse && ncol(X)>0) {
     X<-X[,-(1:ncol(X)),drop=FALSE]
     warning("removing linear terms (including the intercept)' not allowed with 'sparse=TRUE' ")
   }
    
    BBlin<-vector("list", ncol(X[,setdiff(colnames(X),"(Intercept)"),drop=FALSE]))
    if(length(setdiff(colnames(X),"(Intercept)"))>=1){ #se ci sono termini lineari (oltre l'intercetta)
      names(BBlin)<-colnames(X[,setdiff(colnames(X),"(Intercept)"),drop=FALSE])
      for(i in 1:ncol(X[,setdiff(colnames(X),"(Intercept)"),drop=FALSE])){
        xdisegno<-seq(min(X[,setdiff(colnames(X),"(Intercept)"),drop=FALSE][,i]), max(X[,setdiff(colnames(X),"(Intercept)"),drop=FALSE][,i]), l=150)
        BB1<- matrix(xdisegno, ncol=1)
        attr(BB1,"covariate.35")<- xdisegno
        attr(BB1,"covariate.n")<- X[,setdiff(colnames(X),"(Intercept)"),drop=FALSE][,i]
        attr(BB1,"vc")<- FALSE
        attr(BB1, "coef.names")<- attr(BB1,"smoothName") <- colnames(X[,setdiff(colnames(X),"(Intercept)"),drop=FALSE])[i]
        BBlin[[i]] <- BB1
      }
    }
    
    if(length(id.ps)>1 && (length(lambda)!=length(id.ps))) stop("with several ps() terms, a single 'lambda' is allowed in each ps()")
    if(length(lambda)>1 && length(id.ps)==1){
      cv<-TRUE
      lambdas<-lambda
      r.cv<- gcrq.rq.cv(Y, B[[1]], X, tau, monotone=vMonot, concave=vConc, ndx=vNdx, lambda=lambdas, deg=vDeg, dif=vDiff, 
                   var.pen=var.pen, penMatrix=penMatrixList, lambda.ridge=lambda.ridge, dropcList=dropcList, decomList=decomList, 
                   vcList=vcList, dropvcList=dropvcList, centerList=centerList, colmeansB=colMeansB.list, 
                   nfolds=nfolds, foldid=foldid, eps=eps, Bconstr=BB, df.option=df.option, df.nc=df.nc, sparse=sparse)
      lambda<-r.cv[[1]]
    } else {
        cv<-FALSE
    }
    ##ATTENZIONE Ai NOMI.. AGGIUNGER IN FUNZIONE DEI LIVELLI
    #nomiPS.ps
    #[1] "ps(x)"    "ps(x):z1"
    
    names(lambda)<-names(vMonot)<-names(vConc)<- unlist(nomiPS.ps.ok) #nomiVariabPEN #meglio nomiPS.ps
    #-------------------------------------------
    shared.penList1 <- unlist(tapply(shared.penList, unlist(nomiPS.ps.ok), function(.x) if(!.x[1]) seq_len(length(.x)) else rep(100,length(.x))))
    id.shared.pen   <- as.numeric(interaction(shared.penList, unlist(nomiPS.ps.ok), drop=TRUE))+shared.penList1 #era nomiPS.ps
    id.shared.pen   <- id.shared.pen[names(lambda)]
    #-----------------------------------------------
    idCoefGroup <- rep(0:length(B), c(ncol(X),sapply(B, ncol)))
    
    #metti il n.degli smooth piuttosto che 1.
    
    #browser()
    #length(B)= n.termini con penalizz
    #pesiPen deve essere una lista di lunghezza pari a length(B)
    
    
    #fit<-list(pesiPen = as.list(rep(1, length(B))))
    #fit<-list(pesiPen=matrix(1, length(idCoefGroup), length(B))) #serve solo se c'e' lo smoothing adattivo..
    pesiPen <- lapply(1:length(B), function(.x) matrix(1, nrow=length(idCoefGroup[idCoefGroup==.x]), ncol=1))
    fit <- list(pesiPen = pesiPen )
    #browser()
    #-----------------------------------------------
    #se qualche lambda deve essere stimato..
    #-----------------------------------------------
    penMatrixList0 <- penMatrixList
    
    
    # fitter<- if(sparse) get("ncross.rq.fitXBsparse", mode="function") else get("ncross.rq.fitXB", mode="function") 
    
    if(any(lambda<0)){
        sigma2uValues<-matrix(,it.max,length(B))
        lambda<- tapply(lambda, id.shared.pen, mean) #riduci i lambda se c'e' qualcuno condiviso..
        names(lambda) <- names(vMonot) #il tapply() perde i nomi.. quindi recuperali..
        id.lambda.est<-(lambda<0) #which lambdas should be estimated?
        lambda[id.lambda.est]<-lambda0 #.5
        lambdaFixed<-lambda[!id.lambda.est]
        new.rho<-10^3
        K <- unlist(K)
        K.matrix<-matrix(K, nrow=length(K), ncol=length(tau))
        colnames(K.matrix)<-paste(tau)
        h.ok<-h
        myepsV <- rep(1e-06,it.max) #seq(1e-01,1e-06,l=it.max)
        ############################################################################################
        #==================================================================
        #browser()
        sigma2u.Rel<-rep(1, length(B))
        sigma2u<- rep(1, length(B))
        All.lambda<-matrix(NA,it.max,length(B))
        if(!single.lambda && length(B)>1) stop(" tau-specific lambdas allowed only with one smooth ")
        for(j in 1:it.max){ #
          
          expand.lambda<- rep(lambda, tapply(id.shared.pen, id.shared.pen, length))
          if(!single.lambda) expand.lambda <- matrix(expand.lambda, ncol=length(tau))
          #==se adattivo
          #browser()
          if(any(adList>0)){
            for(i in 1:length(B)) {
              if(adList[i]>0){
                #pesij <- rowMeans(abs(fit$pesiPen[idCoefGroup==i, ,drop=FALSE])+.0001)^adList[i]
                pesij <- rowMeans(abs(fit$pesiPen[[i]])+.0001)^adList[i]
                #pesij<- drop(dscad(abs(fit$coefficients[idCoefGroup==i])+.00001, expand.lambda))
                #diag(penMatrixList[[i]]) <- diag(penMatrixList0[[i]])/pesij
                penMatrixList[[i]] <- penMatrixList0[[i]]/pesij
                #penMatrixList[[i]] <- if(group) pesij/sqrt(sum(pesij^2))*penMatrixList0[[i]] else penMatrixList0[[i]]/pesij
              }
            }
          }
          
          if(sparse){
            fit <-suppressWarnings(ncross.rq.fitXBsparse(y=Y, B=B, X=X, taus=tau, monotone=vMonot, concave=vConc, ndx=vNdx,
                                                   lambda=expand.lambda, deg=vDeg, dif=vDiff, var.pen=var.pen, eps=eps, penMatrix=penMatrixList,
                                                   lambda.ridge=lambda.ridge, dropcList=dropcList, decomList=decomList, vcList=vcList, dropvcList=dropvcList,
                                                   centerList=centerList, colmeansB=colMeansB.list, ridgeList=ridgeList, ps.matrix.list=ps.matrix.list,
                                                   Bconstr=BB, df.option=df.option, df.nc=df.nc, adjX.constr=adjX.constr, myeps=myepsV[j], it.j=10,
                                                   adList=adList))
          } else {
            fit <-suppressWarnings(ncross.rq.fitXB(y=Y, B=B, X=X, taus=tau, monotone=vMonot, concave=vConc, ndx=vNdx,
                    lambda=expand.lambda, deg=vDeg, dif=vDiff, var.pen=var.pen, eps=eps, penMatrix=penMatrixList,
                    lambda.ridge=lambda.ridge, dropcList=dropcList, decomList=decomList, vcList=vcList, dropvcList=dropvcList,
                    centerList=centerList, colmeansB=colMeansB.list, ridgeList=ridgeList, ps.matrix.list=ps.matrix.list,
                    Bconstr=BB, df.option=df.option, df.nc=df.nc, adjX.constr=adjX.constr, myeps=myepsV[j], it.j=10,
                    adList=adList))
                  #NB colMeansBorig.list ha un valore in meno per VC terms!
          }

          edf.jM<-edf.j<-fit$edf.j
          sigma2u.old <- sigma2u
          
          if(single.lambda){ #un unico lambda per tutti i percentili..
            dev2u<-apply(fit$all.dev2u,1,sum) #all.dev2u e' matrice . con un singolo tau era "fit$dev2u"
            if(attr(fit$id.coef, "nomi")[1]!="Xlin") edf.jM<-rbind(0, edf.jM) #edf.j<-c(0,edf.j)
            edf.j<-apply(edf.jM,1,sum)
            #sigma2u <-psi.k<- dev2u/(K*edf.j[-1])  #e' un vettore
            sigma2u <-psi.k<- tapply(dev2u, id.shared.pen, sum)/(tapply(K, id.shared.pen, mean)* tapply(edf.j[-1], id.shared.pen, sum))  #e' un vettore
            denom.e<-if(method=="REML") max(1,n*length(tau)-sum(rbind(1,K.matrix)*edf.jM)) else (n*length(tau))
            sigma2e<-phi<-sum(fit$all.dev2e)/denom.e  #sopra, se non vuoi pesare usa "sum(edf.j)"  
          } else {
            if(attr(fit$id.coef, "nomi")[1]!="Xlin") edf.j<-edf.jM<-rbind(0,edf.j)   
            sigma2u <-psi.k.tau<-fit$all.dev2u/(K.matrix*edf.j[-1,]) ##e' una matrice
            denom.e<-if(method=="REML") (n-apply(edf.j,2,sum)) else n
            sigma2e<-phi.tau<-fit$all.dev2e/denom.e 
            sigma2e<-matrix(sigma2e, ncol=length(tau), nrow=nrow(sigma2u), byrow = TRUE)
          }
          
          lambda.old<-lambda
          ###########################################################
          lambda<- sigma2e/(sigma2u+.00001*sigma2e) #scalare, vettore o MATRICE!!! CONTROLLA!!! naturalemnete vedi anche fitXB.r()
          ###########################################################
          #attenzione ai nomi... magari names(lambda)<-names(lambda.old) ??
          
          #if(j==8) browser()
          
          if(!single.lambda) lambda.old<-matrix(lambda.old, ncol=ncol(lambda), nrow=nrow(lambda))
          
          d<-lambda-lambda.old
          if(j>=it.max) h.ok<-h^(j-it.max)
          lambda <- pmin(lambda.old+h.ok*d,lambda.max)
          lambda[!id.lambda.est] <- lambdaFixed #funziona anche con matrici !!! :-))) che culo!!! :-)))
          if(single.lambda) All.lambda[j,]<- lambda
          #---------------------------------------------
          #se vuoi valutare diff in rho..
          old.rho<-new.rho
          new.rho<-fit$rho
          #---------------------------------------------
          if(all(abs((lambda.old-lambda)/(lambda.old+1e-8))< tol)) break
          
          if(single.lambda && j>=5){
            #controlla e fermati se tra gli ultimi 4 valori dei lambda stimati quelli diversi sono <=2 (il che significa 
            #   che si stanno ripetendo sempre gli stessi valori)
            #browser()
            n.diversi.lambda <- apply(round(All.lambda[j:(j-3),, drop=FALSE], 3), 2, function(.x) length(unique(.x)))
            if(all(n.diversi.lambda <=2)) break
          }
          
          #if(all((abs(sigma2u-sigma2u.old)/sum(sigma2u)) < tol)) break
           
          #if() Fellner: change in Variance/sum(AllVariances)
          #sigma2u.Rel.old <- sigma2u.Rel
          #sigma2uValues[j,] <- sigma2u.Rel<- sigma2u /sum(sigma2u)
          #if(all(abs(sigma2u.Rel.old - sigma2u.Rel)< tol)) break
          
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
        lambda<- rep(lambda, tapply(id.shared.pen, id.shared.pen,length)) #espandi lambda

        names(edf.j)<-rownames(edf.jM)<-c("Xlin",nomiVariabPEN)        
        if(attr(fit$id.coef, "nomi")[1]!="Xlin") {
          edf.j<-edf.j[-1] #rimuovi il primo lo zero eventualmente aggiunto prima.. #e' importante altrimenti i nomi non corrispondono..
          edf.jM<-edf.jM[-1,,drop=FALSE]
        }
      } else {
    
        #browser()
        if(any(adList>0)) warning("The adaptive penalty is not implemented with fixed lambda")
        # if(any(adList>0)){
        #   for(i in 1:length(B)) {
        #     if(adList[i]>0){
        #       pesij <- drop(abs(fit$coefficients[idCoefGroup==i])+.0001)^adList[i]
        #       diag(penMatrixList[[i]]) <- diag(penMatrixList0[[i]])/pesij #e' penMatrixList0, giusto? o senza 0?
        #       #pesij<- drop(dscad(abs(fit$coefficients[idCoefGroup==i])+.00001, expand.lambda))
        #       #diag(penMatrixList[[i]]) <- diag(penMatrixList0[[i]])*pesij
        #     }
        #   }
        # }
        if(sparse){
          fit <-suppressWarnings(ncross.rq.fitXBsparse(y=Y, B=B, X=X, taus=tau, monotone=vMonot, concave=vConc, ndx=vNdx,
                                                       lambda=lambda, deg=vDeg, dif=vDiff, var.pen=var.pen, eps=eps, penMatrix=penMatrixList,
                                                       lambda.ridge=lambda.ridge, dropcList=dropcList, decomList=decomList, vcList=vcList, dropvcList=dropvcList,
                                                       centerList=centerList, colmeansB=colMeansB.list, ridgeList=ridgeList, ps.matrix.list=ps.matrix.list,
                                                       Bconstr=BB, df.option=df.option, df.nc=df.nc, adjX.constr=adjX.constr, it.j=10))
        } else {
          fit <-suppressWarnings(ncross.rq.fitXB(y=Y, B=B, X=X, taus=tau, monotone=vMonot, concave=vConc, ndx=vNdx,
                                                 lambda=lambda, deg=vDeg, dif=vDiff, var.pen=var.pen, eps=eps, penMatrix=penMatrixList,
                                                 lambda.ridge=lambda.ridge, dropcList=dropcList, decomList=decomList, vcList=vcList, dropvcList=dropvcList,
                                                 centerList=centerList, colmeansB=colMeansB.list, ridgeList=ridgeList, ps.matrix.list=ps.matrix.list,
                                                 Bconstr=BB, df.option=df.option, df.nc=df.nc, adjX.constr=adjX.constr, it.j=10))
          #NB colMeansBorig.list ha un valore in meno per VC terms!
        }
        
        #browser()
            # fit<-ncross.rq.fitXB(y=Y, B=B, X=X, taus=tau, monotone=vMonot, concave=vConc, ndx=vNdx, 
            #   lambda=lambda, deg=vDeg, dif=vDiff, var.pen=var.pen, eps=eps, penMatrix=penMatrixList, 
            #   lambda.ridge=lambda.ridge,dropcList=dropcList,decomList=decomList, vcList=vcList, dropvcList=dropvcList, 
            #   centerList=centerList, colmeansB=colMeansB.list,ridgeList=ridgeList, ps.matrix.list=ps.matrix.list, Bconstr=BB, 
            #   df.option=df.option, df.nc=df.nc, adjX.constr=adjX.constr, adList=adList, it.j=10)
        
        K.matrix<-NULL
        #K.matrix<-matrix(K, nrow=length(K), ncol=length(tau))
        #colnames(K.matrix)<-paste(tau)
        
      }
        
      #---------------------------
      if(n.boot>0){
        B.orig<-B
        coef.boot<-array(NA, dim=c(nrow(as.matrix(fit$coef)), length(tau), n.boot))
        id.NA<-NULL
        for(i in 1:n.boot) {
          id<-sample(1:n, size=n, replace=TRUE)
          for(j in 1:length(B.orig)) B[[j]]<-B.orig[[j]][id,,drop=FALSE]
          .o.b<- try(ncross.rq.fitXB(y=Y[id], B=B, X=X[id,, drop=FALSE], taus=tau, monotone=vMonot, concave=vConc, ndx=vNdx,
                            lambda=lambda, deg=vDeg, dif=vDiff, var.pen=var.pen, eps=eps, penMatrix=penMatrixList,
                            lambda.ridge=lambda.ridge, dropcList=dropcList, decomList=decomList, vcList=vcList, dropvcList=dropvcList,
                            centerList=centerList, colmeansB=colMeansB.list,ridgeList=ridgeList, ps.matrix.list=ps.matrix.list, 
                            Bconstr=BB, df.option=df.option, df.nc=df.nc, adjX.constr=adjX.constr),
                                silent=TRUE)

          if(!inherits(.o.b,"try-error"))  coef.boot[,,i]<-.o.b$coef else id.NA[length(id.NA)+1]<-i
        } #end i=1,..,n.boot
        if(!is.null(id.NA)) {
          coef.boot<-coef.boot[,,-id.NA]
          warning(paste(length(id.NA), "boot resamples without convergence"), call.=FALSE)
        }
      } #end if n.boot
    
    nomiCoefUNPEN<-rownames(as.matrix(fit$coefficients))[1:ncol(X)]
      nn<-c(nomiCoefUNPEN,unlist(nomiCoefPEN))

      if(is.matrix(fit$coefficients)) rownames(fit$coefficients)<-nn else names(fit$coefficients)<-nn

      names(BB)<-names(dropvcList)<-names(centerList)<-names(decomList)<-names(dropcList)<- names(vNdx)<- names(vDeg)<- 
        names(vDiff)<- names(knotsList)<-nomiVariabPEN
	    fit$info.smooth<-list(monotone=vMonot, ndx=vNdx, lambda=lambda, deg=vDeg, dif=vDiff, concave=vConc, 
                    dropcList=dropcList, decomList=decomList, centerList=centerList, vcList=vcList, dropvcList=dropvcList, 
                    testo.ps=nomiVariabPEN, knots=knotsList, K=K.matrix) 

	    ###########################
      if(is.matrix(lambda)) {
        rownames(lambda)<-nomiVariabPEN
        colnames(lambda)<-paste(tau)
      } else {
        names(lambda)<-nomiVariabPEN
      }
      #########nomina le righe di fit$edf.all e fit$edf.j
	    #NON e' chiaro cosa succede se df.option<=2. fit$edf.all esiste? non dovrebbe essere definito..
	    if(is.matrix(fit$edf.all)) rownames(fit$edf.all)<-nn else names(fit$edf.all)<-nn
	    rownames(fit$edf.j)<- if(nrow(fit$edf.j)==length(nomiVariabPEN)) nomiVariabPEN else c("Xlin",nomiVariabPEN)
      ##############
	    
	    names(Bderiv)<-names(BB)
	    fit$Bderiv<-Bderiv
	    fit$BB<-c(BB, BBlin) #35 (in realt? 100) righe e attr "covariate.n" "covariate.35". #NB attr(,"covariate.n") contiene altri attributi "ndx", "deg", "pdiff", "monot", "lambda","nomeX"
      
      fit$lambda<-lambda
    } #fine del "se ci sono termini smooth"
    #========================================
    #========================================
    attr(fit$edf.j,"df.nc") <- df.nc
    if(y) fit$y<-Y
    if(!is.null(transf)) {
      #if(!is.null(transf)) fit$transf.inv<-transf.inv
      attr(fit$fitted.values, "transf.inv") <- transf.inv
      attr(fit$fitted.values, "transf") <- eval(parse(text=paste("function(y){", transf, "}"))) 
    }
    fit$contrasts <- attrContr
    
    #se il modello e' stato chiamato con  y ~ 0 + ps(x, by = g, dropc = FALSE, center = FALSE), .xlivelli non ha informazioni
    #sulla categoriale nomiBy con categoria levelsBy . Questo puo' essere un problema per predict.gcrq() che deve usare le info
    #di .xlivelli. Quindi "aggiorniamo" .xlivelli
    #browser()
    
    if(length(testo.ps)>0 && any(!sapply(levelsBy, is.null))){
      .xlivelli<-c(.xlivelli, levelsBy)
      .xlivelli <- .xlivelli[unique(names(.xlivelli))]
    }
    
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

