plot.gcrq <-
  function(x, term=NULL, add=FALSE, res=FALSE, conf.level=0, axis.tau=FALSE, interc=TRUE, se.interc=FALSE,
           legend=FALSE, select.tau, deriv=FALSE, cv=FALSE, transf=NULL, lambda0=FALSE, shade=FALSE, 
           overlap=NULL, rug=FALSE, overall.eff=TRUE, grid=NULL, smoos=NULL, split=FALSE, 
           shift=0, type=c("sandw","boot"), ...){ #palette="Roma",
    #===============================================================================
    bspline <- function(x, ndx, xlr = NULL, knots=NULL, deg = 3, deriv = 0, outer.ok=FALSE) {
      # x: vettore di dati
      # xlr: il vettore di c(xl,xr)
      # ndx: n.intervalli in cui dividere il range
      # deg: il grado della spline
      # Restituisce ndx+deg basis functions per ndx-1 inner nodi
      # ci sono "ndx+1" nodi interni + "2*deg" nodi esterni
      # require(splines)
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
    #===============================================================================
    type<- match.arg(type)
    if(conf.level>0) VAR <-vcov.gcrq(x, type=type)  
    b<-x$coefficients
    
    if("terms" %in% names(list(...))) stop(" 'terms' is not allowed.. Do you mean 'term'? ")
    
    if(axis.tau){
      if(!is.matrix(b) || ncol(b)<=1) stop("axis.tau=TRUE is meaningless with a single quantile fit")
      tt<-as.numeric(colnames(x$coefficients)) #i valori di tau 0.10, 0.25,...
      
      if(nrow(b)>1 && split) {
        n.plot<-c(ceiling(nrow(b)/2),2)
        oldpar<-par(mfrow=n.plot)
      }
      oldAsk<-options("device.ask.default")[[1]]
      
      if(is.null(term)){
        term<-1:nrow(b)
      } else {
        if(is.character(term)) term<-  match(term, rownames(b), NA)
      }
      if(any(is.na(term))) stop("undefined term")
      all.terms.id<- term
      for(j in all.terms.id){
        b<-b[j,]
        if(conf.level>0) {
          se<-sapply(VAR, function(.x)sqrt(.x[j,j]))
          z<- -qnorm((1-conf.level)/2)
          val<-cbind(b-z*se, b, b+z*se)
          if(!add){
            matplot(tt, val, lty=c(2,1,2), type="o", pch=20, 
              xlab="Probability values", 
              ylab=paste("coefficient of ", rownames(x$coefficients)[j]),...)
          } else {
            matpoints(tt, val, lty=c(2,1,2), type="o", pch=20, ...)
          }
        } else {
          if(!add) {
            plot(tt, b, type="o", pch=20, xlab="Probability values", 
                ylab=paste("coefficient of ", rownames(x$coefficients)[j]),...)
          } else {
          points(tt, b, type="o", pch=20, ...)
          }
        }
        devAskNewPage(ask =TRUE)
      } #end for(j in all.terms.id)
      
      devAskNewPage(ask =oldAsk)
      if(length(all.terms.id)>1 && split) par(oldpar)
      return(invisible(NULL))      
    }
    if(cv){
      if(is.null(x$cv)) stop("the object does not include the 'cv' component")
      valoriL<-x$cv[,1]
      valoriCV<-apply(x$cv[,-1],1,mean)
      #aa <- apply(x$cv[, -1], 1, quantile, probs=c(.25,.75))  #se vuoi riportare anche i quantili .25, .75
      if(!lambda0) {
        valoriL<-valoriL[-1]
        valoriCV<-valoriCV[-1]
        #aa<-aa[,-1]
      }
      plot(valoriL, valoriCV ,type="o",ylab="Cross Validation score",
           xlab="lambda values",xaxt="n",...)
      axis(1, at=x$cv[,1], labels=round(x$cv[,1],2), las=2, cex.axis=.7)
      points(valoriL[which.min(valoriCV)], min(valoriCV, na.rm=TRUE), pch=19)
      #segments(valoriL, aa[1,], valoriL, aa[2,], lty=3) #comunque devi aggiustare ylim..
      #------------
      #matlines(valoriL, x$cv[-1, -1], type="l", col=grey(.5)) #non saprei..
      #boxplot(unlist(apply(x$cv[,-1],1,function(.x)list(.x)), recursive=FALSE)) #ma non si vede NULLA!!!
      return(invisible(NULL))
    }
    if(!is.null(overlap)) legend<-TRUE
    type<-match.arg(type)
    y<-res

    if(is.null(x$info.smooth)){ #if((x$pLin-match("(Intercept)", rownames(x$coefficients),0))>0)
      if(is.null(term)){
        term<-colnames(x$x)
        term <- setdiff(term, "(Intercept)")
      } else {
        if(is.numeric(term)) term<- colnames(x$x[, term+match("(Intercept)", colnames(x$x),0), drop=FALSE])
      }
      term <- setdiff(term, "(Intercept)")
      if(!is.character(term)) stop("problems in 'term' ")
    } else {
      #per termini smooth
      if(is.null(term)) {
        term=names(x$BB)
      } else {
        if(is.numeric(term)) term <- names(x$BB)[term]
      }
      if(!is.character(term)) stop("problems in 'term' ")
      if(!all(term %in% names(x$BB))) stop(paste("Unknown term. It should be numeric or one of: ",paste(names(x$BB),collapse=" ")))
    }
    
    n.smooths<-length(term)  #include anche termini lineari..
    
    if(n.smooths>1 && split) {
      n.plot<-c(ceiling(n.smooths/2),2)
      oldpar<-par(mfrow=n.plot)
    }
    oldAsk<-options("device.ask.default")[[1]]
    #ATTENZIONE per resettare i parametri grafici, mi sa che alla fine del for devi mettere:
    #devAskNewPage(ask =old1)
    #if(split) par(mfrow=old2)
    #dove
    #   old1<-options("device.ask.default")[[1]]
    #   old2<-par()$mfrow
    if(length(x$taus)<=1) select.tau<-1
    if(missing(select.tau)) {
      select.tau<-1:ncol(x$coefficients)
    } else {
      if(length(select.tau)>length(x$taus)) stop("`length(select.tau)<=length(taus)' is requested")
      if(all(select.tau<1 & select.tau>0)) select.tau<- match(select.tau,x$taus)
      if(!all(select.tau<=length(x$taus) & select.tau>=1)) stop("'select.tau' is not correctly specified")
    }
    
        
    if(deriv) y<-interc<-FALSE
    if(add) y<-FALSE
    
    if(length(shift)>1 && length(shift)!=length(select.tau)) stop(" 'shift' should be a scalar or an appropriate vector")
    
    shift00<-shift  
    interc00<- interc
    all.term.names<-term
    
    #browser()
    
    #if(length(x$info.smooth$lambda)!=n.smooths)
    #eliminare da all.term.names i nomi dei termini lineari?
    #################################################################################################
    for(term in all.term.names){
      #if(term=="x") browser()
      #browser()
      interc <- interc00
      shift<-shift00
      BB<-x$BB[[term]] 
      if(!is.null(x$minX) && term%in%names(x$minX)) BB<- BB - x$minX[[term]] 
      nomi.ok<- attr(BB, "coef.names") ##nomi.ok<-paste(term,"ps",1:ncol(BB),sep=".")
      xvar.n <-attr(BB,"covariate.n") 
      xvar.35 <-attr(BB,"covariate.35") 
      vc.term<-attr(BB,"vc")
      smoothName<-attr(BB,"smoothName") #attr(BB,"smoothName1") ps(x) o ps(x):z
      
      
      #if(isTRUE(vc.term)) {
      #  if(isTRUE(attr(BB,"drop"))){
      #    #if(interc) {BB<-cbind(1,BB)} else {nomi.ok<-nomi.ok[-1]}
      #  } else {
      #    if(!interc) warning(" 'interc=FALSE' ignored with vc terms called with 'dropc=FALSE' ", call.=FALSE)
      #  }
      #  interc<-FALSE  
      #} else {
      #  if(is.null(shift)) shift<-0 #se NON e' un VC e shift=NULL (default), poni shift<-0
      #}
      #shift=0
      
      
      ###################################
      ### ATTENZIONE: con vc models interc e' posto FALSE il che preclude di disegnare l'effetto con la model intercept
      #  comunque da una prova fatta se metti interc=TRUE, poi i codici seguenti funzionano...
      #################################
      
      b<-if(length(x$tau)<=1) drop(x$coefficients)[nomi.ok] else x$coefficients[nomi.ok,select.tau]
      
      fit.35<-if(deriv) x$Bderiv[[term]]%*%b else BB%*%b #matrici
      ###########=================================================
      corr.df<-0
      #browser()
      if(is.character(xvar.35) && xvar.35[1]=="ridge") { #if ridge, plot the single coeffs 
        if(conf.level>0) warning(" 'conf.level>0' not (yet) implemented", call.=FALSE, immediate. = FALSE)
        #browser()
        tt <- sub("\\)","",gsub("ps\\(","" ,term, ))
        if(is.matrix(b)){
          etich <- sub(tt,"",rownames(b), fixed=TRUE)
        } else {
          etich <- sub(tt,"",names(b), fixed=TRUE)
          }
        #b<-b # +shift, nell'altra versione c'era come argomento shift=0.. non so se utile..
        #browser()
        if(is.matrix(b)) {
          xx1 <- matrix(1:nrow(b), nrow(b), ncol(b)) + matrix(seq(0,.3,l=ncol(b)), nrow(b), ncol(b), byrow = TRUE) 
        } else {
          xx1 <- 1:length(b)
        }
          
        #browser()
        
        if(add) {
          if(is.matrix(b)) matpoints(xx1, b, ...) else points(xx1, b, ...)
        } else {
          #sqrt(diag(VAR[[1]]))[attr(BB,"coef.names")] #se vuoi metterci gli IC
           opz<- list(...)
           opz$xaxt <- "n"
           if(is.null(opz$xlab)) opz$xlab<-tt
           if(is.null(opz$ylab)) opz$ylab<-"Estimate"
           opz$x<- xx1
           opz$y<- b #b[,select.tau,drop=FALSE]
           if(is.matrix(b)){
             if(is.null(opz$col)) opz$col<-(1:ncol(b))+1
             if(is.null(opz$pch)) opz$pch<-(1:ncol(b))
             do.call(matplot, opz)
           } else {
             if(is.null(opz$col)) opz$col<-2
             do.call(plot, opz)
           }
          #if(is.matrix(b)) matplot(xx1, b, xlab=tt, xaxt="n", ylab="Estimate", ... ) else plot(xx1, b, xlab=tt, xaxt="n", ylab="Estimate", ...)
          axis(1, at=1:nrow(as.matrix(b)), labels=etich, cex.axis=.5) #labels=attr(BB,"coef.names")
          abline(h=0, lty=3)
          abline(v=1:nrow(as.matrix(b))+.7, col=grey(.8))
          if(legend) {
            opz.leg<-list(pch=opz$pch, col=opz$col, bty="n", legend=formatC(x$taus[select.tau], digits=2, format="f"), cex=.5, x="topright")
            do.call("legend", opz.leg)
            #legend("topright", legend=x$taus, bty="n", cex  = .5, ...)
          }
        }
      } else { #altrimenti disegni del segnale smooth..
  #browser()
        if(!is.null(attr(BB,"name.fixed.params"))){ #se c'e' stata una decomp di spline 
        #if(!missing(n.points)) warning("argument 'n.points' ignored", call. = FALSE)
          if(overall.eff){ #e vuoi disegnare tutto l'effetto  
            nomi.okF<-attr(BB,"name.fixed.params")
            b.fix<-if(is.matrix(x$coefficients)) x$coefficients[nomi.okF,select.tau] else x$coefficients[nomi.okF]
            fit.35<- fit.35 +  drop(poly(xvar.35, degree=length(b.fix), raw=TRUE)%*%b.fix)  
            corr.df <-length(nomi.okF) #serve per correggere i df in Ylab
          } else {
            interc <- FALSE
          }
        }
        
        #blocco commentato il 25/10/21
        #if(is.null(shift) && isTRUE(vc.term)){ 
        #  if("(Intercept)"%in%rownames(as.matrix(x$coefficients))){
        #    fit.35<-fit.35 + matrix(as.matrix(x$coefficients)["(Intercept)", select.tau], ncol=ncol(fit.35), nrow=nrow(fit.35), byrow=TRUE)
        # }
        #  shift<-0
        #}
        #browser()
        
        id.model.interc <- "(Intercept)"%in%rownames(as.matrix(x$coefficients)) 
        if(interc && id.model.interc) {#NB qui interc e' sempre FALSE per termini VC 
          fit.35<-fit.35 + matrix(as.matrix(x$coefficients)["(Intercept)", select.tau], ncol=ncol(fit.35), nrow=nrow(fit.35), byrow=TRUE)
        }
        
        fit.35<-fit.35+ matrix(shift, nrow(fit.35), ncol(fit.35),byrow = TRUE) 
        m.x<- if(is.null(list(...)$xlim)) min(xvar.35) else min(list(...)$xlim)
        M.x<- if(is.null(list(...)$xlim)) max(xvar.35) else max(list(...)$xlim)
        select.n<-xvar.n>=m.x & xvar.n<= M.x
        xvar.n<-xvar.n[select.n]
        select.35<-xvar.35>=m.x & xvar.35<=M.x
        xvar.35<-xvar.35[select.35]
        fit.35<-fit.35[select.35, , drop=FALSE]
        l<-c(list(x=xvar.35, y=fit.35),list(...))
        cexL<-if(is.null(l$cex)) .65 else l$cex #sara' usato solo se legend=TRUE
        if(!is.null(l$text.col)) {
          text.col<- l$text.col
          l$text.col<-NULL
        } else {
          text.col<-1
        }
        Ylab.ok<-all.vars(formula(x))[1]
        if(is.null(transf)){
          f.transf.inv<-attr(x$fitted.values, "transf.inv")
          #if(!is.null(f.transf.inv)) Ylab.ok <- paste(Ylab.ok, " (", x$call$transf,")", sep="")
        } else {
          if(!is.character(transf)) stop(" 'transf' should be NULL or character")
          Ylab.ok <- paste(Ylab.ok, " (", transf,")", sep="")
          f.transf.inv <- eval(parse(text=paste("function(y){", transf, "}"))) #assegna la funzione a partire dal carattere..
        } 
        if(!is.null(f.transf.inv)) l$y<- apply(as.matrix(l$y), 2, f.transf.inv) #l$y <- eval(parse(text=transf), list(y=l$y))               
        if(length(select.tau)==1 && ncol(BB)!=1 ){ #ncol(BB)==1 se il termine e' lineare, altrimenti se e' una base ha molte piu' colonne..
          .df.term <- as.matrix(x$edf.j)[term, select.tau, drop=FALSE]
          #if(isTRUE(vc.term)) .df.term <- .df.term-1  #se e' un VC term, dal conteggio df togli l'intercetta (cosi come avviene per i termini smooth non-VC)
          Ylab.ok<- gsub(")", paste(", df=", round(.df.term[term,]+corr.df,2),")",sep=""), term)
        } else {
          Ylab.ok <- paste("Effect of ", term)
        }
        #se col<0
        
        #browser()
        
        if(!is.null(l$col) && !is.character(l$col) && all(l$col < 0)){ 
          Lab.palette <- colorRampPalette(c("blue", "green", "red"), space = "Lab")      #c("blue", "orange", "red")
          l$col<-Lab.palette(length(select.tau))
          #terrain.colors(length(select.tau), alpha = .6) #anche rainbow() or heat.colors(..)
        }
        ####per disegnare IC
        if(conf.level>0){
          alpha <- 1 - conf.level
          a <- alpha/2
          a <- c(a, 1 - a)
          fac <- qnorm(a)
          #se vc e interc=FALSE la BB ha una colonna in meno!! (certo perche' interc=FALSE ha portato alla rimozione ...)
          
          if(id.model.interc && se.interc) {
            nomi.ok <- c("(Intercept)", nomi.ok)
            BB <- cbind(1, BB)
          }

          V<- lapply(VAR, function(x) x[nomi.ok, nomi.ok])
          se<-sapply(V, function(x)sqrt(rowSums((BB %*% x * BB)))) 
          low.q <- fit.35 + fac[1]*se[select.35, select.tau, drop=FALSE]
          up.q <-  fit.35 + fac[2]*se[select.35, select.tau, drop=FALSE]
          if(shade){
            yy <- sapply(seq.int(ncol(fit.35)), function(.i) c(low.q[, .i], tail(up.q[, .i], 1), rev(up.q[, .i]), low.q[, .i][1]))
            xx <- c(xvar.35, tail(xvar.35, 1), rev(xvar.35), xvar.35[1])
            l2<-l
            l2$x<-xx
            l2$y<-yy
            #l2 <- c(list(x=xx, y=yy), list(...))
            l2$col <- if(is.null(l2$col)) adjustcolor(grey(.3), alpha.f = 0.25) else adjustcolor(l2$col, alpha.f = 0.25)
            if(is.null(l2$border)) l2$border <- NA
          } else {
            l2 <- c(list(x = xvar.35, y = low.q), list(...))
            l3 <- c(list(x = xvar.35, y = up.q), list(...))
            l2$type<-l3$type<-"l"
            if(is.null(l2$cex)) l2$cex <- l3$cex <- .25
          }
        }          
        #######
        #browser()
        if(y){
          id.res.ok<-which.min(abs(x$taus[select.tau]-.5))
          ff<-splinefun(xvar.35, fit.35[,id.res.ok])
          #x$y <- as.matrix(x$residuals)[, id.res.ok] + ff(xvar.n)
          x$y<- (as.matrix(x$residuals)[select.n,select.tau,drop=FALSE])[,id.res.ok] + ff(xvar.n)
          l1<-c(list(x=xvar.n, y=x$y),list(...))
          l1$cex<-NULL
          l1$text.col<-NULL
          #if(!is.null(transf)) l1$y <- eval(parse(text=transf), list(y=l1$y))              
          if(!is.null(f.transf.inv)) l1$y<- apply(as.matrix(l1$y), 2, f.transf.inv) #l1$y <- f.transf.inv(l1$y)
          if(is.null(l1$xlab)) l1$xlab<- smoothName #term
          if(is.null(l1$ylab)) l1$ylab<- Ylab.ok
          if(!is.null(l1$lwd)) l1$lwd<-NULL                            
          if(!is.null(l1$cex.p)) l1$cex<-l1$cex.p else l1$cex<-.75
          l1$cex.p<-NULL
          if(!is.null(l1$col.p)) l1$col<-l1$col.p else l1$col<-grey(.35, alpha = .15)
          l1$col.p<-NULL
          if(!is.null(l1$pch.p)) l1$pch<-l1$pch.p else l1$pch=19
          l1$pch.p<-NULL
          #if(!is.null(l1$col)) coll<-l1$col; l1$col<-NULL               
          if(is.null(l1$xlim)) l1$xlim <- c(min(xvar.35), max(xvar.35))  #era xvar.n
          l1$x[l1$x>max(l1$xlim)]<-NA #metti NA prima di incrementare il limite
          if(legend && is.null(overlap)) l1$xlim <- l1$xlim*c(1,1.03) 
          if(is.null(l1$ylim)){
            if(conf.level>0) {
              l1$ylim <- if(shade) range(c(l2$y, l1$y)) else c(min(c(l2$y,l1$y)), max(c(l3$y,l1$y)))
            } else {
              l1$ylim <- range(c(l$y,l1$y))
            }
          }
          if(is.null(smoos)) { smoos <- if(length(l1$x)>10000) TRUE else FALSE }
          if(smoos){
            l1$type<-"n"
            do.call(plot, l1)
            smoothScatter(l1$x, l1$y, add=TRUE, nrpoints = 0, colramp= colorRampPalette(c("white", grey(.4))))
          #browser()
            } else {
            do.call(plot, l1)              
          }
          if(!is.null(grid)){
            xval<- if(length(grid$x)==1) seq(par()$usr[1], par()$usr[2] , l= grid$x+2)[-c(1,grid$x+2)] else grid$x
            yval<- if(length(grid$y)==1) seq(par()$usr[3], par()$usr[4] , l= grid$y+2)[-c(1,grid$y+2)] else grid$y
            if(is.null(grid$col)) grid$col<-grey(.7)
            if(is.null(grid$lty)) grid$lty<- 3
            if(is.null(grid$lwd)) grid$lwd<- .7
            abline(h=yval, v=xval, col=grid$col, lty=grid$lty, lwd=grid$lwd)
          }
          add<-TRUE
        } #end if(res)
        if(legend){
          if(!is.null(overlap)) { #se xlim e' fornito come modificare??
            #if(!(is.numeric(overlap) && length(overlap)==1)) stop(" 'overlap' should be NULL or a numeric scalar")
            #xleg <- l$x[which.min(abs(l$x-overlap))] #
            #yleg <- l$y[which.min(abs(l$x-overlap)),]
            if(!is.numeric(overlap)) stop(" 'overlap' should be NULL or numeric (scalar/vector)")
            overlap<-rep(overlap, length.out=ncol(l$y))
            xleg <- sapply(overlap, function(.x)l$x[which.min(abs(l$x-.x))] )
            yleg <- sapply(overlap, function(.x)l$y[which.min(abs(l$x-.x)),] )  
          } else {
            #l$xlim <- if(is.null(l$xlim)) c(min(xvar.n),1.1*max(xvar.n)) else c(min(l$xlim),1.1*max(l$xlim))
            #xleg<-1.05*max(xvar.n)
            if(is.null(l$xlim)) l$xlim <- c(min(xvar.35),max(xvar.35)) 
            l$xlim<- l$xlim*c(1, 1.03) 
            xleg<-max(l$xlim) *1.01
            #yleg<- l$y[100,] #yleg<-tail(l$y,1)
            yleg<-tail(l$y[l$x<max(l$xlim), ],1)
          }
        }
        M.x<- if(is.null(l$xlim)) max(xvar.35) else max(l$xlim) 
        l$x[l$x>M.x]<-NA
        l$col.p<-NULL
        l$cex.p<-NULL
        l$pch.p<-NULL
        #browser()
        if(is.null(l$col)) l$col<- 2#hcl.colors(length(select.tau),alpha = .6, palette="Roma") #seq(2, length.out=length(select.tau))
        if(is.null(l$lwd)) l$lwd<-1.6
        if(is.null(l$lty)) l$lty<-1
        if(add){
          if(!is.null(l$col)) l$col<-adjustcolor(l$col,.65)
          do.call(matlines, l)
        } else {
          if(is.null(l$xlab)) l$xlab<- smoothName # term
          if(is.null(l$ylab)) {l$ylab<-if(deriv) paste(Ylab.ok," (first derivative)") else Ylab.ok}
          l$type<-"l"
          if(conf.level>0 & is.null(l$ylim)) {
            if(shade) l$ylim <- range(l2$y) else l$ylim <- c(min(l2$y), max(l3$y))
          }
          do.call(matplot, l)
          if(!is.null(grid)){
            xval<- if(length(grid$x)==1) seq(par()$usr[1], par()$usr[2] , l= grid$x+2)[-c(1,grid$x+2)] else grid$x
            yval<- if(length(grid$y)==1) seq(par()$usr[3], par()$usr[4] , l= grid$y+2)[-c(1,grid$y+2)] else grid$y
            if(is.null(grid$col)) grid$col<-grey(.7)
            if(is.null(grid$lty)) grid$lty<- 3
            if(is.null(grid$lwd)) grid$lwd<- .7
            abline(h=yval, v=xval, col=grid$col, lty=grid$lty, lwd=grid$lwd)
          }
          #if(y && !is.null(x$y)) points(xvar.n, x$y)
          #matplot(xvar.35, fit.35, type="l", xlab=term, ylab="", ...)
        }
        if(rug) { #rug(xvar)
          segments(xvar.n, rep(par()$usr[3],length(xvar.n)), xvar.n,
                  rep(par()$usr[3],length(xvar.n))+ abs(diff(par()$usr[3:4]))/40)
        }
        
        if(conf.level>0){
          if(shade){
            l2$cex.p<- l2$col.p<-l2$pch.p<-NULL
            if(length(l2$col)!=length(select.tau)) l2$col<-rep(l2$col, length(select.tau))
            sapply(seq.int(ncol(fit.35)), function(.i){
              l3 <- l2; l3$y <- l3$y[,.i]; l3$col <- l3$col[.i]; do.call("polygon", l3)
            })
          } else {
            l2$col<-l3$col<-l$col
            l2$lwd<-l3$lwd<-l$lwd*.5
            l2$lty<-l3$lty<-if(length(l$lty)==1) l$lty+1 else l$lty
            if(is.null(l2$lty)) l2$lty<-2
            if(is.null(l3$lty)) l3$lty<-2
            l2$x<-  l3$x<- l$x
            do.call(matlines, l2)
            do.call(matlines, l3) 
          }
        }
        #browser()
        if(legend) {
          text.col<-rep(text.col, l=ncol(yleg))
          if(is.null(overlap)) {
            text(xleg, yleg, formatC(x$taus[select.tau], digits=2, format="f"), cex=cexL, col=text.col)
            text(xleg, yleg, formatC(x$taus[select.tau], digits=2, format="f"), cex=cexL, col=text.col)
          } else {
            for(i in 1:ncol(yleg)) {
              legend(xleg[i], yleg[i,i], formatC(x$taus[select.tau][i], digits=2, format="f"), 
                    bty="o", xjust=.5, yjust=.5, bg=adjustcolor("white",.4), box.col=adjustcolor("white",.4),
                    adj=.5, cex=cexL, text.col=text.col[i])
            }
            for(i in 1:ncol(yleg)) { 
              legend(xleg[i], yleg[i,i], formatC(x$taus[select.tau][i], digits=2, format="f"), 
                    bty="n", xjust=.5, yjust=.5, adj=.5, cex=cexL, text.col=text.col[i])
            }
          }
        }
      devAskNewPage(ask =TRUE)
      }
    } #end for(term in all.terms)
    devAskNewPage(ask =oldAsk)
    if(n.smooths>1 && split) par(oldpar)
    #invisible(NULL) #necessario?          
  }  
