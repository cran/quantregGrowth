plot.gcrq <-
function(x, add=FALSE, y=FALSE, legend=FALSE, select.tau, deriv=FALSE, cv=FALSE,...){ #, se=FALSE, intercept=FALSE, resid=TRUE, alpha=0.01, legend=TRUE, ...){
#x: un oggetto restituito da gcrq()
#add: se TRUE aggiunge le linee.....
#y: se TRUE e se l'oggetto x contiene y (i dati) allora li disegna. Se add=TRUE, y viene posto a FALSE
#... argomenti da passare a plot() e a matplot/matlines e text(). Quindi se c'è col questo viene applicato sia a plot (se i dati devono essere disegnati) 
#   sia a matlines/matplot se lwd solo a matlines() se cex a plot e a text per legend.
#  alla legenda (se legend=TRUE) e sia a plot() 
#select.tau: which quantile curves should be drawn? default (missing) is all
#f.deriv: if TRUE, the first derivatives ofthe growth curves  are plotted
#cv: se TRUE disegna la cross validation versus lambdas
          if(cv){
               if(is.null(x$cv)) stop("the object does not include the 'cv' component")
               plot(x$cv[,1], apply(x$cv[,-1],1,mean),type="o",ylab="Cross Validation score",
                xlab="lambda values",xaxt="n")
               axis(1, at=x$cv[,1], labels=round(x$cv[,1],2), las=2)
               #boxplot(unlist(apply(m$cv[,-1],1,function(x)list(x)), recursive=FALSE))
               return(invisible(NULL))
              }
          if(is.null(x$BB)) stop(" plot.gcrq() only works with a single smooth variable")          
          if(missing(select.tau)) {
              select.tau<-1:ncol(x$coefficients)
              } else {
              if(length(select.tau)>length(x$taus)) stop("`length(select.tau)<=length(taus)' is requested")
              if(all(select.tau<1 & select.tau>0)) select.tau<- match(select.tau,x$taus)
              if(!all(select.tau<=length(x$taus) & select.tau>=1)) stop("select.tau is not correctly specified")
              }
          if(length(x$BB)==1) term<-names(x$BB) else stop("which (smooth) term?")
          if(deriv) y<-FALSE
          if(add) y<-FALSE
          BB<-x$BB[[term]]
          xvar.n<-attr(BB,"covariate.n")
          xvar.35<-attr(BB,"covariate.35")
          nomi.ok<-paste(term,"ps",1:ncol(BB),sep=".")
          b<-x$coefficients[nomi.ok,select.tau]
          fit.35<-if(deriv) drop(x$Bderiv%*%b) else drop(BB%*%b)
          l<-c(list(x=xvar.35, y=fit.35),list(...))
          if(y && is.null(x$y)) warning("y=TRUE ignored.. the fit does not include the data", call.=FALSE)
          if(y && !is.null(x$y)) {
              l1<-c(list(x=xvar.n, y=x$y),list(...))
              if(is.null(l1$xlab)) l1$xlab<-term
              if(is.null(l1$ylab)) l1$ylab<-"Growth variable"
              if(legend) l1$xlim <- c(min(xvar.n),1.1*max(xvar.n))
              do.call(plot, l1)              
              #plot(xvar.n, x$y, xlab=term, ylab="Growth variable")
              if(legend) {
                  cexL<-if(is.null(l1$cex)) .6 else l1$cex
                  text(1.05*max(xvar.n),  l$y[nrow(l$y),], x$taus[select.tau], cex=cexL)
                  legend<-FALSE
                  }
              add<-TRUE
              }
          if(legend) {
            x.leg<-l$x[(length(l$x)-8)]
            l$x[(length(l$x)-11):(length(l$x)-5)]<-NA
          }
          if(add){
              do.call(matlines, l)
              } else {
              if(is.null(l$xlab)) l$xlab<-term
              if(is.null(l$ylab)) {l$ylab<-if(deriv) "Growth variable (first derivative)" else "Growth variable"}
              l$type<-"l"
              do.call(matplot, l)
              #if(y && !is.null(x$y)) points(xvar.n, x$y)
              #matplot(xvar.35, fit.35, type="l", xlab=term, ylab="", ...)
              }
          if(legend)  {
              cexL<-if(is.null(l$cex)) .6 else l$cex
              text(x.leg,  l$y[(length(l$x)-8),], x$taus[select.tau], cex=cexL)
            #tau<-x$tau; mtext(bquote(tau == .(tau)),line=-2)
            }
          }
