plot.gcrq <-
function(x, term, add=FALSE, res=FALSE, conf.level=0, legend=FALSE, select.tau, deriv=FALSE, cv=FALSE, 
  transf=NULL, lambda0=FALSE, shade=FALSE, overlap=FALSE, rug=FALSE, n.points=100, grid=NULL, smoos=NULL,...){ #, se=FALSE, intercept=FALSE, resid=TRUE, alpha=0.01, legend=TRUE, ...){
#x: un oggetto restituito da gcrq()
#add: se TRUE aggiunge le linee.....
#y: se TRUE e se l'oggetto x contiene y (i dati) allora li disegna. Se add=TRUE, y viene posto a FALSE
#... argomenti da passare a plot() e a matplot/matlines e text(). Quindi se c'e' col questo viene applicato sia a plot (se i dati devono essere disegnati) 
#   sia a matlines/matplot se lwd solo a matlines() se cex a plot e a text per legend.
#  alla legenda (se legend=TRUE) e sia a plot() 
#select.tau: which quantile curves should be drawn? default (missing) is all
#f.deriv: if TRUE, the first derivatives ofthe growth curves  are plotted
#cv: se TRUE disegna la cross validation versus lambdas
#===============================================================================
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
#===============================================================================
          y<-res
          if(is.null(x$info.smooth)) stop("plot for simple linear fits is not allowed")
          if(is.null(x$BB)) stop(" plot.gcrq() only works with smooth terms")          
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
               points(valoriL[which.min(valoriCV)], min(valoriCV), pch=19)
               #segments(valoriL, aa[1,], valoriL, aa[2,], lty=3) #comunque devi aggiustare ylim..
               #------------
               #matlines(valoriL, x$cv[-1, -1], type="l", col=grey(.5)) #non saprei..
               #boxplot(unlist(apply(x$cv[,-1],1,function(.x)list(.x)), recursive=FALSE)) #ma non si vede NULLA!!!
               return(invisible(NULL))
              }
          if(length(x$BB)>1) {if(missing(term)) stop("please provide 'term'")} else {term=names(x$BB)}
          if(!term%in%names(x$BB)) stop("'term' is not a smooth variable")
          if(length(x$tau)<=1) select.tau<-1
          if(missing(select.tau)) {
              select.tau<-1:ncol(x$coefficients)
              } else {
              if(length(select.tau)>length(x$taus)) stop("`length(select.tau)<=length(taus)' is requested")
              if(all(select.tau<1 & select.tau>0)) select.tau<- match(select.tau,x$taus)
              if(!all(select.tau<=length(x$taus) & select.tau>=1)) stop("'select.tau' is not correctly specified")
              }
          if(deriv) y<-FALSE
          if(add) y<-FALSE
          BB<-x$BB[[term]]
          xvar.n<-attr(BB,"covariate.n")
          xvar.35<-attr(BB,"covariate.35")
          ###########=================================================
          if(n.points!=100){
                 xvar.35<- seq(min(xvar.35),max(xvar.35), length=n.points)
                 BB<-bspline(xvar.35, ndx=attr(BB,"ndx"), deg= attr(BB,"deg"))
                 }
          ###########=================================================
          nomi.ok<-paste(term,"ps",1:ncol(BB),sep=".")
          b<-if(length(x$tau)<=1) x$coefficients[nomi.ok] else x$coefficients[nomi.ok,select.tau]
          fit.35<-if(deriv) x$Bderiv%*%b else BB%*%b #matrici
          if("(Intercept)"%in%rownames(as.matrix(x$coefficients))) {
                fit.35<-fit.35 + matrix(as.matrix(x$coefficients)["(Intercept)",], ncol=ncol(fit.35), nrow=nrow(fit.35), byrow=TRUE)
                }
          m.x<- if(is.null(list(...)$xlim)) min(xvar.35) else min(list(...)$xlim)
          M.x<- if(is.null(list(...)$xlim)) max(xvar.35) else max(list(...)$xlim)
          select.n<-xvar.n>=m.x & xvar.n<= M.x
          xvar.n<-xvar.n[select.n]
          select.35<-xvar.35>=m.x & xvar.35<=M.x
          xvar.35<-xvar.35[select.35]
          fit.35<-fit.35[select.35, , drop=FALSE]

          l<-c(list(x=xvar.35, y=fit.35),list(...))
          cexL<-if(is.null(l$cex)) .6 else l$cex #sara' usato solo se legend=TRUE
          
          
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
            
          
          #se col<0
          if(!is.null(l$col) && !is.character(l$col) && l$col < 0){ 
              Lab.palette <- colorRampPalette(c("blue", "green", "red"), space = "Lab")      #c("blue", "orange", "red")
              l$col<-Lab.palette(length(select.tau))
              }

          ####per disegnare IC
          if(conf.level>0){
              alpha <- 1 - conf.level
              a <- alpha/2
              a <- c(a, 1 - a)
              fac <- qnorm(a)
              V<-vcov.gcrq(x,term)
              se<-sapply(V, function(x)sqrt(rowSums((BB %*% x * BB))))
              low.q <- fit.35 + fac[1]*se[select.35, select.tau, drop=FALSE]
              up.q <-  fit.35 + fac[2]*se[select.35, select.tau, drop=FALSE]
              if (shade) {
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
          if(y) {
              id.res.ok<-which.min(abs(x$taus[select.tau]-.5))
              ff<-splinefun(xvar.35, fit.35[,id.res.ok])
              #x$y <- as.matrix(x$residuals)[, id.res.ok] + ff(xvar.n)
              x$y<- (as.matrix(x$residuals)[select.n,select.tau,drop=FALSE])[,id.res.ok] + ff(xvar.n)
              l1<-c(list(x=xvar.n, y=x$y),list(...))
              l1$cex<-NULL
              #if(!is.null(transf)) l1$y <- eval(parse(text=transf), list(y=l1$y))              
              if(!is.null(f.transf.inv)) l1$y<- apply(as.matrix(l1$y), 2, f.transf.inv) #l1$y <- f.transf.inv(l1$y)
              if(is.null(l1$xlab)) l1$xlab<-term
              if(is.null(l1$ylab)) l1$ylab<- Ylab.ok
              if(!is.null(l1$lwd)) l1$lwd<-NULL                            
              if(!is.null(l1$col)) l1$col<-NULL                            
              if(!is.null(l1$col.p)) l1$col<-l1$col.p;l1$col.p<-NULL
              if(!is.null(l1$cex.p)) l1$cex<-l1$cex.p;l1$cex.p<-NULL
              if(!is.null(l1$pch.p)) l1$pch<-l1$pch.p;l1$pch.p<-NULL              

              if(is.null(l1$xlim)) l1$xlim <- c(min(xvar.35), max(xvar.35))  #era xvar.n
              l1$x[l1$x>max(l1$xlim)]<-NA #metti NA prima di incrementare il limite
              if(legend && !overlap) l1$xlim <- l1$xlim*c(1,1.1) 
              if(conf.level>0 & is.null(l1$ylim)) {
                      l1$ylim <- if(shade) range(c(l2$y, l1$y)) else c(min(c(l2$y,l1$y)), max(c(l3$y,l1$y)))
                              }
#browser()
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
#smoothScatter(x, nrpoints = 0)
#smoothScatter(x)
#smoothScatter(x, colramp= colorRampPalette(c("white", "black")))
              } #end if(res)
              
              if(legend) {
                         if(overlap) { #se xlim e' fornito come modificare??
                            xleg<-l$x[floor(.92*n.points)]  #xleg<-rev(l$x[l$x<max(l$xlim)])[8]
                            yleg<-l$y[floor(.92*n.points),]
                            l$x[floor(.89*n.points):floor(.95*n.points)]<-NA
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
#browser()          
          if(add){
                do.call(matlines, l)
              } else {
                if(is.null(l$xlab)) l$xlab<-term
                if(is.null(l$ylab)) {l$ylab<-if(deriv) paste(Ylab.ok," (first derivative)") else Ylab.ok}
                l$type<-"l"
                l$col.p<-NULL
                l$cex.p<-NULL
                l$pch.p<-NULL
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
        if(rug) {
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
               l2$x<-  l3$x<- l$x
               do.call(matlines, l2)
               do.call(matlines, l3) 
               }
          }

          if(legend) {
             
             text(xleg, yleg, formatC(x$taus[select.tau], digits=2, format="f"), cex=cexL)
             }
          }
