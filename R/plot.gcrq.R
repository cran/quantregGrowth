plot.gcrq <-
  function(x, term=NULL, add=FALSE, res=FALSE, conf.level=0, interc=TRUE, legend=FALSE, select.tau, deriv=FALSE, cv=FALSE, 
           transf=NULL, lambda0=FALSE, shade=FALSE, overlap=NULL, rug=FALSE, overall.eff=TRUE, #n.points=100,  
           grid=NULL, smoos=NULL, split=FALSE, shift=NULL, type=c("sandw","boot"), ...){
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
    type<-match.arg(type)
    y<-res
    if(is.null(x$info.smooth)){ #stop("plot for simple linear fits is not allowed")
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
    
    n.smooths<-length(term)  
    
    if(n.smooths>1 && split) {
      n.plot<-c(ceiling(n.smooths/2),2)
      oldpar<-par(mfrow=n.plot)
    }
    #in realta' oldpar sotto non serve perche' alla fine viene richiamato solo se split=TRUE
    #} else { 
    #oldpar<-par()
    #}
    oldAsk<-options("device.ask.default")[[1]]
    #ATTENZIONE per resettare i parametri grafici, mi sa che alla fine del for devi mettere:
    #devAskNewPage(ask =old1)
    #if(split) par(mfrow=old2)
    #dove
    #   old1<-options("device.ask.default")[[1]]
    #   old2<-par()$mfrow
    #  if(n.smooths>1) {if(missing(term)) stop("please provide 'term'")} else {term=names(n.smooths)}
    #  if(!term%in%names(x$BB)) stop("'term' is not a smooth variable")
    if(length(x$taus)<=1) select.tau<-1
    if(missing(select.tau)) {
      select.tau<-1:ncol(x$coefficients)
    } else {
      if(length(select.tau)>length(x$taus)) stop("`length(select.tau)<=length(taus)' is requested")
      if(all(select.tau<1 & select.tau>0)) select.tau<- match(select.tau,x$taus)
      if(!all(select.tau<=length(x$taus) & select.tau>=1)) stop("'select.tau' is not correctly specified")
    }
    
    if(conf.level>0) VAR <-vcov.gcrq(x, type=type)      
    if(deriv) y<-FALSE
    if(add) y<-FALSE
    
    shift00<-shift  
    interc00<- interc
    all.term.names<-term
    #################################################################################################
    for(term in all.term.names){
      interc <- interc00
      shift<-shift00
      BB<-x$BB[[term]] 
      if(!is.null(x$constX)) BB<- BB - x$constX[[term]] 
      nomi.ok<- attr(BB, "coef.names") ##nomi.ok<-paste(term,"ps",1:ncol(BB),sep=".")
      xvar.n <-attr(BB,"covariate.n") 
      xvar.35 <-attr(BB,"covariate.35") 
      vc.term<-attr(BB,"vc")
      smoothName<-attr(BB,"smoothName") #attr(BB,"smoothName1") ps(x) o ps(x):z
      
      if(isTRUE(vc.term)) {
        if(isTRUE(attr(BB,"drop"))){
          if(interc) {BB<-cbind(1,BB)} else {nomi.ok<-nomi.ok[-1]}
        } else {
          if(!interc) warning(" 'interc=FALSE' ignored with vc terms called with 'dropc=FALSE' ", call.=FALSE)
        }
        interc<-FALSE  
      } else {
        if(is.null(shift)) shift<-0 #se NON e' un VC e shift=NULL (default), poni shift<-0
      }
      ###################################
      ### ATTENZIONE: con vc models interc e' posto FALSE il che preclude di disegnare l'effetto con la model intercept
      #  comunque da una prova fatta se metti interc=TRUE, poi i codici seguenti funzionano...
      #################################
      
      b<-if(length(x$tau)<=1) drop(x$coefficients)[nomi.ok] else x$coefficients[nomi.ok,select.tau]
      fit.35<-if(deriv) x$Bderiv%*%b else BB%*%b #matrici
      ###########=================================================
      if(is.character(xvar.35) && xvar.35[1]=="ridge") { #if ridge, plot the single coeffs **UNTESTED**
        b<-b # +shift, nell'altra versione c'era come argomento shift=0.. non so se utile..
        if(add) {
          points(1:length(b), b,...)
          return(invisible(NULL))
        }
        plot(b, xlab=term, xaxt="n")
        axis(1, at=1:length(b), labels=rownames(BB), cex.axis=.6,...)
        return(invisible(NULL))
      }
      if(!is.null(attr(BB,"name.fixed.params"))){ #se c'e' stata una decomp di spline 
        #if(!missing(n.points)) warning("argument 'n.points' ignored", call. = FALSE)
        if(overall.eff){ #e vuoi disegnare tutto l'effetto  
          nomi.okF<-attr(BB,"name.fixed.params")
          b.fix<-if(is.matrix(x$coefficients)) x$coefficients[nomi.okF,select.tau] else x$coefficients[nomi.okF]
          fit.35<- fit.35 +  drop(poly(xvar.35, degree=length(b.fix), raw=TRUE)%*%b.fix)  
        }
      }
      if(is.null(shift) && isTRUE(vc.term)){ 
        if("(Intercept)"%in%rownames(as.matrix(x$coefficients))){
          fit.35<-fit.35 + matrix(as.matrix(x$coefficients)["(Intercept)", select.tau], ncol=ncol(fit.35), nrow=nrow(fit.35), byrow=TRUE)
        }
        shift<-0
      }
      if(interc && "(Intercept)"%in%rownames(as.matrix(x$coefficients))) {#NB qui interc e' sempre FALSE per termini VC 
        fit.35<-fit.35 + matrix(as.matrix(x$coefficients)["(Intercept)", select.tau], ncol=ncol(fit.35), nrow=nrow(fit.35), byrow=TRUE)
      }
      fit.35<-fit.35+shift
      m.x<- if(is.null(list(...)$xlim)) min(xvar.35) else min(list(...)$xlim)
      M.x<- if(is.null(list(...)$xlim)) max(xvar.35) else max(list(...)$xlim)
      select.n<-xvar.n>=m.x & xvar.n<= M.x
      xvar.n<-xvar.n[select.n]
      select.35<-xvar.35>=m.x & xvar.35<=M.x
      xvar.35<-xvar.35[select.35]
      fit.35<-fit.35[select.35, , drop=FALSE]
      l<-c(list(x=xvar.35, y=fit.35),list(...))
      cexL<-if(is.null(l$cex)) .65 else l$cex #sara' usato solo se legend=TRUE
      Ylab.ok<-all.vars(formula(x))[1]
      if(is.null(transf)){
        f.transf.inv<-attr(x$fitted.values, "transf.inv")
        #if(!is.null(f.transf.inv)) Ylab.ok <- paste(Ylab.ok, " (", x$call$transf,")", sep="")
      } else {
        if(!is.character(transf)) stop(" 'transf' should be NULL or character")
        Ylab.ok <- paste(Ylab.ok, " (", transf,")", sep="")
        f.transf.inv<- eval(parse(text=paste("function(y){", transf, "}"))) #assegna la funzione a partire dal carattere..
      } 
      if(!is.null(f.transf.inv)) l$y<- apply(as.matrix(l$y), 2, f.transf.inv) #l$y <- eval(parse(text=transf), list(y=l$y))               
      if(length(select.tau)==1 && ncol(BB)!=1 ){ #ncol(BB)==1 se il termine e' lineare, altrimenti se e' una base ha molte piu' colonne..
        .df.term <- as.matrix(x$edf.j)[term, select.tau, drop=FALSE]
        if(isTRUE(vc.term)) .df.term <- .df.term-1  #se e' un VC term, dal conteggio df togli l'intercetta (cosi come avviene per i termini smooth non-VC)
        Ylab.ok<- gsub(")", paste(", df=", round(.df.term[term,],2),")",sep=""), term)
      } else {
        Ylab.ok <- term
      }
      #se col<0
      if(!is.null(l$col) && !is.character(l$col) && l$col < 0){ 
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
      if(y){
        id.res.ok<-which.min(abs(x$taus[select.tau]-.5))
        ff<-splinefun(xvar.35, fit.35[,id.res.ok])
        #x$y <- as.matrix(x$residuals)[, id.res.ok] + ff(xvar.n)
        x$y<- (as.matrix(x$residuals)[select.n,select.tau,drop=FALSE])[,id.res.ok] + ff(xvar.n)
        l1<-c(list(x=xvar.n, y=x$y),list(...))
        l1$cex<-NULL
        #if(!is.null(transf)) l1$y <- eval(parse(text=transf), list(y=l1$y))              
        if(!is.null(f.transf.inv)) l1$y<- apply(as.matrix(l1$y), 2, f.transf.inv) #l1$y <- f.transf.inv(l1$y)
        if(is.null(l1$xlab)) l1$xlab<- smoothName #term
        if(is.null(l1$ylab)) l1$ylab<- Ylab.ok
        if(!is.null(l1$lwd)) l1$lwd<-NULL                            
        if(!is.null(l1$cex.p)) l1$cex<-l1$cex.p else l1$cex<-.5
        l1$cex.p<-NULL
        if(!is.null(l1$col.p)) l1$col<-l1$col.p else l1$col<-grey(.7)
        l1$col.p<-NULL
        if(!is.null(l1$pch.p)) l1$pch<-l1$pch.p else l1$pch=20
        l1$pch.p<-NULL
        #if(!is.null(l1$col)) coll<-l1$col; l1$col<-NULL               
        if(is.null(l1$xlim)) l1$xlim <- c(min(xvar.35), max(xvar.35))  #era xvar.n
        l1$x[l1$x>max(l1$xlim)]<-NA #metti NA prima di incrementare il limite
        if(legend && is.null(overlap)) l1$xlim <- l1$xlim*c(1,1.1) 
        if(conf.level>0 & is.null(l1$ylim)) {
          l1$ylim <- if(shade) range(c(l2$y, l1$y)) else c(min(c(l2$y,l1$y)), max(c(l3$y,l1$y)))
        }
        if(is.null(smoos)) { smoos <- if(length(l1$x)>10000) TRUE else FALSE }
        if(smoos){
          l1$type<-"n"
          do.call(plot, l1)
          smoothScatter(l1$x, l1$y, add=TRUE, nrpoints = 0, colramp= colorRampPalette(c("white", grey(.4))))
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
          if(!(is.numeric(overlap) && length(overlap)==1)) stop(" 'overlap' should be NULL or a numeric scalar")
          xleg <- l$x[which.min(abs(l$x-overlap))]
          yleg <- l$y[which.min(abs(l$x-overlap)),]
          
          #xleg<-l$x[floor(.92*n.points)]        #xleg<-rev(l$x[l$x<max(l$xlim)])[8]
          #yleg<-l$y[floor(.92*n.points),]
          #l$x[floor(.89*n.points):floor(.95*n.points)]<-NA
        } else {
          #l$xlim <- if(is.null(l$xlim)) c(min(xvar.n),1.1*max(xvar.n)) else c(min(l$xlim),1.1*max(l$xlim))
          #xleg<-1.05*max(xvar.n)
          if(is.null(l$xlim)) l$xlim <- c(min(xvar.35),max(xvar.35)) 
          l$xlim<- l$xlim*c(1, 1.1) 
          xleg<-.985*max(l$xlim)
          #yleg<- l$y[100,] #yleg<-tail(l$y,1)
          yleg<-tail(l$y[l$x<max(l$xlim), ],1)
        }
      }
      M.x<- if(is.null(l$xlim)) max(xvar.35) else max(l$xlim) 
      l$x[l$x>M.x]<-NA
      l$col.p<-NULL
      l$cex.p<-NULL
      l$pch.p<-NULL
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
      if(legend) {
        if(is.null(overlap)) {
          text(xleg, yleg, formatC(x$taus[select.tau], digits=2, format="f"), cex=cexL)
        } else {
          for(i in 1:length(yleg)) {
            legend(xleg, yleg[i], formatC(x$taus[select.tau][i], digits=2, format="f"), 
                   bty="o", xjust=.5, yjust=.5, bg=adjustcolor("white",.4), box.col=adjustcolor("white",.4),
                   adj=.5, cex=cexL)
          }
          for(i in 1:length(yleg)) { 
            legend(xleg, yleg[i], formatC(x$taus[select.tau][i], digits=2, format="f"), 
                   bty="n", xjust=.5, yjust=.5, adj=.5, cex=cexL)
          }
        }
      }
      devAskNewPage(ask =TRUE)
    }
    devAskNewPage(ask =oldAsk)
    if(n.smooths>1 && split) par(oldpar)
    #invisible(NULL) #necessario?          
  }  
