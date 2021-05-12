ncross.rq.fitXB <-
  function(y, x, B=NULL, X=NULL, taus, monotone=FALSE, concave=FALSE, #tolto argomento interc=FALSE,
           nomiBy=NULL, byVariabili=NULL,
           ndx=10, deg=3, dif=3, lambda=0, eps=.0001, var.pen=NULL, penMatrix=NULL,
           lambda.ridge=0,  dropcList=FALSE, decomList=FALSE, vcList=FALSE, dropvcList=FALSE, centerList=FALSE, colmeansB=NULL, ...){
    #Stima dei non-crossing rq, possibly monotone
    #B: la base di spline o una lista di Bspline. Se NULL viene costruita attraverso la variabile x
    #x: la variabile o la matrice di variabili rispetto a cui viene (o vengono) costruita la base, ammesso che B sia NULL.
    #nomiBy, byVariabili:  nomi e variabili (matrice)  rispetto a cui viene (o vengono) costruita la base. Inutili se B e' fornita.  
    #err.rho: se TRUE usa la sigma2 basata su rho, altrimenti quella basata sui residui al quadrato pesati
    #u.rho: if TRUE, returns the (squared) sum of absolute values of the penalized coeffs
    #
    #plott {0,1,2} se 0 non disegna, se 1 aggiunge se 2 apre un nuovo device. Se x non ? fornita
    #   plott viene posto a 0.
    plott=0 #prima era un argomento
    adj.middle=FALSE
    err.rho=TRUE #prima era un argomento
    u.rho=TRUE ##prima era un argomento
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
    build.D<-function(var.pen.ok, p.ok, dif.ok, lambda.ok, penMatrix.ok=NULL, dropc=FALSE, decomp=FALSE, dropvc=FALSE){
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
        
        D.ok<-if(decomp) xx.var.pen*diag(p.ok) else xx.var.pen*diff(diag(p.ok), diff=dif.ok)
        #D.ok<-if(dif.ok<=0) xx.var.pen*diag(p.ok) else xx.var.pen*diff(diag(p.ok), diff=dif.ok)
      }
      if(dif.ok>0 && dropc) D.ok<-D.ok[,-1] # D%*%C infatti deve risultare t(C)%*%Penalty%*%C
      if(dropvc) D.ok<- cbind(0,D.ok) #rbind(0,cbind(0,D.ok))
      D.ok<-lambda.ok*D.ok
      D.ok
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
    pSmooth<- sum(all.p)
    H<-length(B) #no. of smooth terms
    B<-matrix(unlist(B), n, pSmooth) #matrix(unlist(B), nrow=3,byrow = FALSE)
    XB<-cbind(X,B)
    p<-ncol(XB) #all coeffs, lineari + smooth
    pLin<-p-pSmooth #pSmooth=p1 #n. termini lineari
    #id.interc<- "(Intercept)" %in% colnames(XB)
    #id.smooth<-c(rep(FALSE, pLin),rep(TRUE,pSmooth))
    #if(!is.null(colmeansB)) 
    #listaMedieB<- lapply(1:length(colmeansB), function(.x){if(is.null(colmeansB[[.x]])) rep(0, all.p[.x]) else colmeansB[[.x]]})
    #all.means.B <- unlist(listaMedieB) #con basi noncentrate e' 0,0...
    
    #browser()
    
    #fallo per tutte le basi.. tanto quando non e' centrato la media e' zero..
    for(j in 1:H){
      colmeansB[[j]] <- -colmeansB[[j]]
      if(vcList[j] && dropvcList[j]) colmeansB[[j]] <- c(1, colmeansB[[j]])
    }
    
    
    

    if(length(dropvcList)!=length(all.p)) stop("errore..") #non dovrebbe servire..
    id.intercVC <- tapply(rep(1:length(all.p), all.p), rep(1:length(all.p), all.p), function(.x)c(TRUE, rep(FALSE, length(.x)-1)))
    id.intercVC <- unlist(id.intercVC) & rep(dropvcList, all.p)
    
    id.intc.VC<- which((id.intercVC)) #le positzioni delle intercette dei VC (quando ci sono gruppi)
    exist.VC <-sum(vcList) #se >0 ci sono VC
    #which.VC <-which(vcList)
    
    id.interc <- match("(Intercept)", colnames(XB), 0)
    if(("(Intercept)" %in% colnames(XB)) && length(id.interc)<=0) stop("errore") #un piccolo controllo..
    
    id.interc.Constr <- (id.interc>0 && any(!vcList))
    id.VC.Constr <- (exist.VC>0 && any(centerList))
      
    #n.righe = p
    #n.colonne = n. di termini lineari + n.smooth. Se non ci sono termini lineari la 1st colonna e' una colonna di FALSE..
    id.Matr <- sapply(0:length(all.p), function(.x) I(c(rep(0, pLin), rep(1:length(all.p), all.p))==.x))
    
#    if(any(!vcList)){ #se ci sono termini smooth non VC (indipendentemente dalla intercetta!)
#      idsmoothInter <-rowSums(id.Matr[,-1, drop=FALSE][, !dropvcList, drop=FALSE])==1
#      idsmoothInter[id.interc] <- TRUE #questa riga non fa niente se id.interc = 0 (cioe se non c'e' interc) 
#      id.Matr <- cbind(idsmoothInter, id.Matr[,dropvcList]) 
#    } else {
#      id.Matr <- id.Matr[,-1]
#    }
    
    
    if(id.interc.Constr){
      idsmoothInter <-rowSums(id.Matr[,-1, drop=FALSE][, !dropvcList, drop=FALSE])==1
      idsmoothInter[id.interc] <- TRUE #questa riga non fa niente se id.interc = 0 (cioe se non c'e' interc)
      valueNcrosInterc <- c(1, unlist(colmeansB[!vcList])) #intercetta si assume nella posiz 1.. colmeansB e' 0 per basi noncentrate
    }
    if(exist.VC) {
      id.Matr <- id.Matr[,-1][, vcList, drop=FALSE]
      valueNcrosVC <- colmeansB[vcList]
    }

    
    
    #QUANTI VINCOLI CI SONO???? quelli di sotto sono uguali?
    #ncol(id.Matr)
    #sum(sapply(valueNcrosList, length)>0)
    #NB ncol(id.Matr)==length(colmeansB)
    
    id.start.tau<-which.min(abs(taus-0.5))
    lambda.tau.spec<-FALSE
    if(is.matrix(lambda)) {
      lambda.tau.spec<-TRUE
      lambda.matrix<-lambda
      colnames(lambda.matrix)<-paste(taus)
      lambda<-lambda.matrix[,id.start.tau]
    }
    
    
    D.list<-D.list.lambda<-D1<-D2<-vector("list", length=H)
    for(j in 1:H){
      
      if(monotone[j]!=0) {
        D1[[j]]<-if(dropcList[j]) sign(monotone[j])*(diff(diag(all.p[j]+1), diff=1)[,-1]) else sign(monotone[j])*diff(diag(all.p[j]), diff=1) 
      } else { D1[[j]]<-matrix(0, all.p[j]-1, all.p[j]) }
      
      if(concave[j]!=0) {
        D2[[j]]<-if(dropcList[j]) sign(-concave[j])*(diff(diag(all.p[j]+1), diff=2)[,-1]) else sign(-concave[j])*diff(diag(all.p[j]), diff=2) 
      } else { D2[[j]]<-matrix(0, all.p[j]-2, all.p[j]) }
      
      #D1[[j]]<-if(monotone[j]!=0) sign(monotone[j])*diff(diag(all.p[j]), diff=1) else matrix(0,all.p[j]-1,all.p[j])
      #D2[[j]]<-if(concave[j]!=0) sign(-concave[j])*diff(diag(all.p[j]), diff=2) else matrix(0,all.p[j]-2,all.p[j])
      
      D.list[[j]] <-build.D(var.pen[j], all.p[j]+1*dropcList[j], dif[j], lambda.ok=1, penMatrix[[j]], 
                            dropc=dropcList[j], decomp=decomList[j], dropvc = dropvcList[j])
      D.list.lambda[[j]]<- D.list[[j]]*lambda[j]
    }      
    
    R.monot<-if(length(D1)<=1) D1[[1]] else do.call("blockdiag",D1)
    R.monot<-cbind(matrix(0,nrow=nrow(R.monot), ncol=pLin), R.monot)
    if(any(concave!=0)) {
      R.conc<-if(length(D2)<=1) D2[[1]] else do.call("blockdiag",D2)
      R.conc<-cbind(matrix(0,nrow=nrow(R.conc), ncol=pLin), R.conc)
      R.monot<-rbind(R.monot,R.conc)
    }
    #r.monot<-rep(0, sum(all.p-1))
    r.monot<-rep(0, nrow(R.monot))
    #matrice D *NON* comprende lambda
    D.matrix<-if(length(D.list)<=1) D.list[[1]] else do.call("blockdiag",D.list) 
    D.matrix<- blockdiag(diag(rep(0,pLin), ncol=pLin),D.matrix)
    #matrice D che comprende lambda. D.list.lambda prima era DD
    D.matrix.lambda<-if(length(D.list.lambda)<=1) D.list.lambda[[1]] else do.call("blockdiag",D.list.lambda) 
    D.matrix.lambda<- blockdiag(diag(rep(0,pLin), ncol=pLin),D.matrix.lambda)
    
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
    o.start<-if(any(monotone!=0 | concave!=0)) rq.fit(x=XB,y=y,tau=start.tau,method="fnc",R=R.monot,r=r.monot)
      else rq.fit(x=XB, y=y, tau=start.tau)
    o.start$rho<-sum(Rho(o.start$residuals[1:n], start.tau))
    id.coef.spline<-vector("list", length=H+1)
    id.coef.spline[[1]]<- if(ncol(X)>0) 1:ncol(X) else 0
    for(j in 1:H) id.coef.spline[[j+1]] <- max(id.coef.spline[[j]]) + seq_len(all.p[j])
    names(id.coef.spline)<-c("Xlin", names(lambda))
    #nomi.ok<-paste(names(lambda)[k],"ps",1:all.p[k],sep=".")
    sigma2u<-vector(length=H)
    for(j in 1:H){
      #sigma2u[j]<-drop(crossprod((DD[[j]]/max(lambda[j],1e-15))%*%o.start$coefficients[id.coef.spline[[j+1]]]))
      sigma2u[j]<- if(u.rho) (sum(abs(D.list[[j]]%*%o.start$coefficients[id.coef.spline[[j+1]]])))
      else sum((D.list[[j]]%*%o.start$coefficients[id.coef.spline[[j+1]]])^2)
      #sigma2u[j]<- sum((D.list[[j]]%*%o.start$coefficients[id.coef.spline[[j+1]]])^2)
      #sigma2u[j]<- sum(abs(D.list[[j]]%*%o.start$coefficients[id.coef.spline[[j+1]]]))
    }
    sigma2e<-if(err.rho) o.start$rho else sum(o.start$residuals[1:n]^2*abs(start.tau-I(o.start$residuals[1:n]<0))) #in realt? ? la dev..
    #--------calcola i df
    id.df<-rep(1:length(id.coef.spline), sapply(id.coef.spline,length))
    attr(id.df, "nomi")<-c("Xlin", names(lambda))
    if(length(id.coef.spline[[1]])==1 && id.coef.spline[[1]]==0) {
      id.df<-id.df[-1]
      attr(id.df, "nomi")<-names(lambda)
    }
    #UNO o PIU' QUANTILI?
    if(length(taus)<=1){ #se length(taus)==1
      all.COEF<-as.matrix(o.start$coef) #era all.COEF<- o.start$coef 
      colnames(all.COEF)<-paste(taus)
      #all.df<- sum(round(o.start$residuals[1:n],2)==0)
      all.rho<-sum(Rho(o.start$residuals[1:n], start.tau)) 
      r<-list(coefficients=all.COEF, x=XB, rho=all.rho, #tolto df=all.df,
              fitted.values=o.start$fitted.values[1:n],residuals=o.start$residuals[1:n],
              dev2u=sigma2u, dev2e=sigma2e)
      all.dev2u=matrix(sigma2u, ncol=1, dimnames=list(names(lambda), paste(taus)))
      all.dev2e=matrix(sigma2e, ncol=1, dimnames=list(NULL, paste(taus)))
      r$all.dev2u<-all.dev2u
      r$all.dev2e<-all.dev2e
      #r$DF.NEG <- r$DF.POS<- NULL
    } else { #se length(taus)>1
      DF.NEG <- DF.POS<- NULL
      Ident<-diag(p)
      COEF.POS<-COEF.NEG<-FIT.POS<-FIT.NEG<-RES.POS<-RES.NEG<-NULL
      df.pos.tau<-df.neg.tau<-rho.pos.tau<-rho.neg.tau<-NULL
      sigma2u.pos.tau<-sigma2u.neg.tau<-NULL
      sigma2e.pos.tau<-sigma2e.neg.tau<-NULL
      #===========================================================================
      if(n.pos.taus>0){
        rho.pos.tau <-df.pos.tau <- vector(length=n.pos.taus)
        COEF.POS<-matrix(,ncol(XB),n.pos.taus)
        colnames(COEF.POS)<-paste(pos.taus)
        b.start<-o.start$coef
        if(any(monotone!=0 | concave!=0)){
          RR<-rbind(Ident,R.monot)
          rr<-c(b.start + eps, r.monot)
        } else {
          RR<- Ident
          rr<- b.start + eps
        }
        #DF per NONcrossing
        DF.POS <- matrix(,nrow(RR),n.pos.taus)
        colnames(DF.POS) <- paste(pos.taus)
        
        #NUOVO APRILE 2021: se c'e' intercetta e spline centrate devi cambiare il vincolo
        if(id.interc.Constr){ #se ci sono smooth con interc
          RR[1,idsmoothInter] <- valueNcrosInterc
          rr[1] <- sum(valueNcrosInterc * b.start[idsmoothInter]) #b.start[1]-sum(all.means.B*b.start[id.smooth])
        }
        if(id.VC.Constr) {#se ci sono VC 
          for(i in 1:length(id.intc.VC)){
            RR[id.intc.VC[i], id.Matr[,i]] <- valueNcrosVC[[i]]
            rr[id.intc.VC[i]] <-  sum( valueNcrosVC[[i]] * b.start[id.Matr[,i]] )
          }
        }

        #if(id.interc>0){
        #  RR[1,idsmoothInter] <- valueNcrosList[[1]]
        #  rr[1] <- sum(valueNcrosList[[1]]*b.start[idsmoothInter]) #b.start[1]-sum(all.means.B*b.start[id.smooth])
        #}

        #if(id.interc){
        #  rr[1] <- b.start[1]-sum(all.means.B*b.start[id.smooth])
        #  RR[1, id.smooth]<- -all.means.B
        #  #rr[1] <- b.start[1]-sum(all.means.B*b.start[-1])
        #  #RR[1,]<-c(1, -all.means.B)
        #}
        #se global penalty
        #if(monotone!=0){
        #RR<-rbind(Ident,R)
        # rr<-c(b.start + eps, rep(0,p1-1))
        # } else {
        # RR<- Ident
        # rr<- b.start + eps
        # }
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
          o<-rq.fit(x=XB,y=y,tau=pos.taus[i],method="fnc",R=RR,r=rr)
          
          o$rho<-sum(Rho(o$residuals[1:n], pos.taus[i]))
          #estrai fitted e residuals
          FIT.POS[,i]<-o$fitted.values[1:n]
          RES.POS[,i]<-o$residuals[1:n]
          #estrai la f. obiettivo
          #df.pos.tau[i] <- sum(abs(o$residuals[1:n])<=.000001)
          rho.pos.tau[i] <- o$rho #sarebbe sum(Rho(o$residuals[1:n], pos.taus[i]))
          #estrai sigma2u
          for(j in 1:H){
            sigma2u.pos.tau[j,i]<-if(u.rho) (sum(abs(drop(D.list[[j]]%*%o$coefficients[id.coef.spline[[j+1]]]))))
            else sum(drop(D.list[[j]]%*%o$coefficients[id.coef.spline[[j+1]]])^2)
          }
          
          #estrai sigma2e
          sigma2e.pos.tau[i]<-if(err.rho) o$rho else sum(o$residuals[1:n]^2*abs(pos.taus[i]-I(o$residuals[1:n]<0)))
          
          #estrai i coeff per "aggiornare" i vincoli rr
          b.start<-o$coef
          COEF.POS[,i]<-b.start
          #NON ANCORA IMPLEMENTATO: confrontare dove i vincoli NONCROSSING (MA ANCHE MONOT/CONC) sono attivi
          #vedi come considerare l'epsilon.. Inoltre la matrice RR include anche monot/conc
          #in particolare i primi RR[1:p,]b.start rr[1:p] sono per il NONCROSSING
          DF.POS[,i] <- 1*(drop(RR%*%b.start-rr)<=1e-8)

          rr<- if(any(monotone!=0 | concave!=0)) c(b.start+eps, r.monot) else b.start + eps
          
          #NUOVO APRILE 2021: se c'e' intercetta e spline centrate devi cambiare il vincolo
          if(id.interc.Constr){ #se ci sono smooth con interc
            rr[1] <- sum(valueNcrosInterc * b.start[idsmoothInter]) #b.start[1]-sum(all.means.B*b.start[id.smooth])
          }
          if(id.VC.Constr) {#se ci sono VC 
            for(i in 1:length(id.intc.VC)){
              rr[id.intc.VC[i]] <-  sum( valueNcrosVC[[i]] * b.start[id.Matr[,i]] )
            }
          }
          
          #if(id.interc) rr[1] <- b.start[1]-sum(all.means.B*b.start[id.smooth])
          #if(id.interc>0) rr[1] <- sum(valueNcrosList[[1]]*b.start[idsmoothInter]) 
        }#end for
      }#end if(n.pos.taus>0)
      ### PER I tau "negativi"..
      if(n.neg.taus>0){
        rho.neg.tau <-df.neg.tau <- vector(length=n.pos.taus)
        COEF.NEG<- matrix(,ncol(XB),n.neg.taus)
        colnames(COEF.NEG)<-paste(neg.taus)
        b.start<-o.start$coef
        neg.taus<-sort(neg.taus,TRUE)
        if(any(monotone!=0 | concave!=0)){
          Ident<-diag(p)
          RR<-rbind(-Ident,R.monot)
          rr<-c(-b.start + eps, r.monot)
        } else {
          RR<- -Ident
          rr<- -b.start + eps
        }
        #DF per NONcrossing
        DF.NEG<- matrix(,nrow(RR),n.neg.taus)
        colnames(DF.NEG)<-paste(neg.taus)
        
        #NUOVO APRILE 2021: se c'e' intercetta e spline centrate devi cambiare il vincolo
        if(id.interc.Constr){ #se ci sono smooth con interc
          RR[1,idsmoothInter] <- -valueNcrosInterc
          rr[1] <- -sum(valueNcrosInterc * b.start[idsmoothInter]) #b.start[1]-sum(all.means.B*b.start[id.smooth])
        }
        if(id.VC.Constr) {#se ci sono VC 
          for(i in 1:length(id.intc.VC)){
            RR[id.intc.VC[i], id.Matr[,i]] <- -valueNcrosVC[[i]]
            rr[id.intc.VC[i]] <-  -sum( valueNcrosVC[[i]] * b.start[id.Matr[,i]] )
          }
        }
        
        
        #if(id.interc>0){
        #  RR[1,idsmoothInter] <- -valueNcrosList[[1]]
        #  rr[1] <- -sum(valueNcrosList[[1]]*b.start[idsmoothInter]) #b.start[1]-sum(all.means.B*b.start[id.smooth])
        #}
        
        #if(id.interc>0){
        #  rr[1] <- -(b.start[1]-sum(all.means.B*b.start[id.smooth]))
        #  RR[1, id.smooth]<- all.means.B
        #  #rr[1] <- -(b.start[1]-sum(all.means.B*b.start[-1]))
        #  #RR[1,]<- -c(1, -all.means.B)
        #} 
        
        FIT.NEG<-RES.NEG<-matrix(,n,n.neg.taus)
        sigma2u.neg.tau<-matrix(,H,n.neg.taus)
        sigma2e.neg.tau<-vector(length=n.neg.taus)
        
        for(i in 1:n.neg.taus){
          #AGGIORNA XB PER INCLUDERE i tau-specific lambda
          if(lambda.tau.spec){
            lambda <- lambda.matrix[,paste(neg.taus[i])]
            for(j in 1:H) D.list.lambda[[j]]<- D.list[[j]]*lambda[j]
            D.matrix.lambda<-if(length(D.list.lambda)<=1) D.list.lambda[[1]] else do.call("blockdiag",D.list.lambda) 
            D.matrix.lambda<-blockdiag(diag(rep(0,pLin), ncol=pLin),D.matrix.lambda)
            if(any(lambda>0)) XB<-rbind(XB.orig, D.matrix.lambda)
            if(lambda.ridge>0) XB<-rbind(XB, lambda.ridge*diag(ncol(XB))) #a small ridge penalty
          }
          o<-rq.fit(x=XB,y=y,tau=neg.taus[i],method="fnc",R=RR,r=rr)
          o$rho<-sum(Rho(o$residuals[1:n], neg.taus[i]))
          FIT.NEG[,i]<-o$fitted.values[1:n]
          RES.NEG[,i]<-o$residuals[1:n]
          #df.neg.tau[i] <- sum(abs(o$residuals[1:n])<=.000001)
          rho.neg.tau[i] <- o$rho #sum(Rho(o$residuals[1:n], neg.taus[i]))
          #estrai sigma2u
          for(j in 1:H){
            sigma2u.neg.tau[j,i]<- if(u.rho) (sum(abs(drop(D.list[[j]]%*%o$coefficients[id.coef.spline[[j+1]]]))))
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
          
          #NUOVO APRILE 2021: se c'e' intercetta e spline centrate devi cambiare il vincolo
          if(id.interc.Constr){ #se ci sono smooth con interc
            rr[1] <- -sum(valueNcrosInterc * b.start[idsmoothInter]) #b.start[1]-sum(all.means.B*b.start[id.smooth])
          }
          if(id.VC.Constr) {#se ci sono VC 
            for(i in 1:length(id.intc.VC)){
              rr[id.intc.VC[i]] <-  -sum( valueNcrosVC[[i]] * b.start[id.Matr[,i]] )
            }
          }

          #rr[1] <- -(b.start[1]-sum(all.means.B*b.start[-1]))
          #if(id.interc) rr[1] <- -(b.start[1]-sum(all.means.B*b.start[id.smooth]))
          #if(id.interc.Constr) rr[1] <- sum(valueNcrosList[[1]]*b.start[idsmoothInter])
          #if(id.interc>0) rr[1] <- -sum(valueNcrosList[[1]]*b.start[idsmoothInter])
        }#end for
      }#end if(n.neg.taus>0)
      #------------------------------
      #      if(adj.middle){
      #            if(monotone!=0){
      #              RR<-rbind(Ident,-Ident,R)
      #              rr<-c(COEF.NEG[,1],-COEF.POS[,1], rep(0, p-1))
      #              } else {
      #              RR<-rbind(Ident, -Ident)
      #              rr<-c(COEF.NEG[,1],-COEF.POS[,1])
      #              }
      #      o.start<-rq.fit(x=XB,y=y,tau=start.tau,method="fnc",R=RR,r=rr)
      #        }
      #-------------------------------
      #browser()
      all.COEF<-cbind(COEF.NEG[,n.neg.taus:1], o.start$coef, COEF.POS)
      colnames(all.COEF)<-paste(taus)
      all.FIT<-cbind(FIT.NEG[,n.neg.taus:1,drop=FALSE], o.start$fitted.values[1:n], FIT.POS)
      colnames(all.FIT)<-paste(taus)
      all.RES<-cbind(RES.NEG[,n.neg.taus:1,drop=FALSE], o.start$residuals[1:n], RES.POS)
      colnames(all.RES)<-paste(taus)
      #all.df<- c(df.neg.tau[n.neg.taus:1], sum(abs(o.start$residuals[1:n])<=.000001), df.pos.tau)
      all.rho<-c(rho.neg.tau[n.neg.taus:1], sum(Rho(o.start$residuals[1:n], start.tau)) , rho.pos.tau)
      all.sigma2e<-c(sigma2e.neg.tau[n.neg.taus:1], sigma2e, sigma2e.pos.tau)
      all.sigma2u <-cbind(sigma2u.neg.tau[,n.neg.taus:1,drop=FALSE], sigma2u, sigma2u.pos.tau)
      colnames(all.sigma2u)<-paste(taus)      
      rownames(all.sigma2u)<-names(lambda)
      r<-list(coefficients=all.COEF,x=XB, rho=all.rho, fitted.values=all.FIT, residuals=all.RES,
              dev2u=sigma2u, dev2e=sigma2e, all.dev2u=all.sigma2u, all.dev2e=all.sigma2e, DF.NEG=DF.NEG, DF.POS=DF.POS) #
    } #end se length(taus)>1
    #browser()
    #Attenzione XB include nella 2nd parte la penalizzazione relativa all'*ultima* curva stimata (quindi puo' essere un "problema" se i lambda sono diversi per quantili)
    r$id.coef <- id.df
    r$D.matrix<-D.matrix
    #cosa mettere?:  D.list, D.list.lambda, D.matrix, D.matrix.lambda
    r$D.matrix.lambda<-D.matrix.lambda
    r$pLin<-pLin      
    r$R.monot<-R.monot
    r$lambda<- if(lambda.tau.spec) lambda.matrix else lambda
    return(r)
  }

