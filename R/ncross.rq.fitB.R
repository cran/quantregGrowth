ncross.rq.fitB <-
function(y, x, B=NULL, taus, monotone=FALSE, adj.middle=FALSE, ndx=10, lambda=0,
  deg=3, dif=3, plott=0, var.pen=NULL, eps=.0001, ...){
#manca l'eps
#stima in cui NON ci sono altre esplicative lineari.
#Stima dei non-crossing rq, possibly monotone
#A differenza di ncross.rq.fit1() questa usa semplici B-spline con un linear inequality constraint on the
#   B-spline coefficients
#B: la base di spline, se NULL viene costruita attraverso la variabile x
#x: la variabile rispetto a cui viene costruita la base, ammesso che B sia NULL. Quando B è fornita
#   questa viene usata per la stima, e x viene usata solo per disegnare, ammesso che plott>0
#plott {0,1,2} se 0 non disegna, se 1 aggiunge se 2 apre un nuovo device. Se x non è fornita
#   plott viene posto a 0.
#Aggiungere una matrice di esplicative lineari?
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
#      require(quantreg)
      n<-length(y)
      if(is.null(B)) B<-bspline(x, ndx=ndx, deg=deg)
      if(missing(x)) plott<-0
      p<-ncol(B)
      Ident<-diag(p)

      taus<-sort(taus)

      id.start.tau<-which.min(abs(taus-0.5))
      start.tau<-taus[id.start.tau]

      pos.taus<-taus[(taus-start.tau)>0]
      neg.taus<-taus[(taus-start.tau)<0]
      n.pos.taus<-length(pos.taus)
      n.neg.taus<-length(neg.taus)

      #if varying penalty
      if(is.null(var.pen)){
        xx.var.pen <- rep(1,(p-dif))
        } else {
      f.var.pen <- function(k) eval(parse(text = var.pen))
      xx.var.pen <- 1:(p-dif)
      xx.var.pen <- sqrt(f.var.pen(max(xx.var.pen)))
      }

      if(monotone!=0) {
          R<-sign(monotone)*diff(diag(p), diff=1)
          DD<- xx.var.pen*diff(diag(p), diff=dif) #ridge: diag(c(0,rep(1,p-1)))
          B<-rbind(B, lambda*DD)
          y<-c(y, rep(0,nrow(DD)))
          o.start<-rq.fit(x=B,y=y,tau=start.tau,method="fnc",R=R,r=rep(0,p-1))
          } else {
              DD<-xx.var.pen*diff(diag(p), diff=dif) #ridge: diag(p) #diff(diag(p), diff=1)
              B<-rbind(B, lambda*DD)
              y<-c(y, rep(0,nrow(DD)))
              o.start<-rq.fit(x=B,y=y,tau=start.tau)
              }

    if(length(taus)<=1){ #se length(taus)==1
      all.COEF<-o.start$coef
      #colnames(all.COEF)<-paste(taus)
      all.df<- sum(round(o.start$residuals[1:n],2)==0)
      all.rho<-sum(Rho(o.start$residuals[1:n], start.tau)) 
      r<-list(coefficients=all.COEF,B=B, df=all.df, rho=all.rho, 
          fitted.values=o.start$fitted.values[1:n],residuals=o.start$residuals[1:n])
      } else { #se length(taus)>1
    COEF.POS<-COEF.NEG<-FIT.POS<-FIT.NEG<-RES.POS<-RES.NEG<-NULL
    df.pos.tau<-df.neg.tau<-rho.pos.tau<-rho.neg.tau<-NULL
    if(n.pos.taus>0){
      rho.pos.tau <-df.pos.tau <- vector(length=n.pos.taus)
      COEF.POS<-matrix(,ncol(B),n.pos.taus)
      colnames(COEF.POS)<-paste(pos.taus)
      b.start<-o.start$coef
      #se global penalty
      if(monotone!=0){
          RR<-rbind(Ident,R)
          rr<-c(b.start, rep(0,p-1))
          } else {
          RR<- Ident
          rr<- b.start
          }
      FIT.POS<-RES.POS<-matrix(,n,n.pos.taus)
      for(i in 1:n.pos.taus){
            o<-rq.fit(x=B,y=y,tau=pos.taus[i],method="fnc",R=RR,r=rr)
            FIT.POS[,i]<-o$fitted.values[1:n]
            RES.POS[,i]<-o$residuals[1:n]
            #estrai la f. obiettivo
            df.pos.tau[i] <- sum(round(o$residuals[1:n],3)==0)
            rho.pos.tau[i] <- sum(Rho(o$residuals[1:n], pos.taus[i]))
            b.start<-o$coef
            COEF.POS[,i]<-b.start
            rr<- if(monotone) c(b.start+eps, rep(0,p-1) ) else b.start+eps
            }#end for
      }#end if(n.pos.taus>0)

    if(n.neg.taus>0){
      rho.neg.tau <-df.neg.tau <- vector(length=n.pos.taus)
      COEF.NEG<-matrix(,ncol(B),n.neg.taus)
      colnames(COEF.NEG)<-paste(neg.taus)
      b.start<-o.start$coef
      neg.taus<-sort(neg.taus,TRUE)
      if(monotone!=0){
          RR<-rbind(-Ident,R)
          rr<-c(-b.start, rep(0, p-1) )
          } else {
          RR<- -Ident
          rr<- -b.start
          }
      FIT.NEG<-RES.NEG<-matrix(,n,n.neg.taus)
      for(i in 1:n.neg.taus){
            o<-rq.fit(x=B,y=y,tau=neg.taus[i],method="fnc",R=RR,r=rr)
            FIT.NEG[,i]<-o$fitted.values[1:n]
            RES.NEG[,i]<-o$residuals[1:n]
            df.neg.tau[i] <- sum(round(o$residuals[1:n],3)==0)
            rho.neg.tau[i] <- sum(Rho(o$residuals[1:n], neg.taus[i]))
            b.start<-o$coef
            COEF.NEG[,i]<-b.start
            rr<- if(monotone) c(-b.start+eps, rep(0,p-1) ) else -b.start + eps
            }#end for
      }#end if(n.neg.taus>0)
#------------------------------
      if(adj.middle){
            if(monotone){
              RR<-rbind(Ident,-Ident,R)
              rr<-c(COEF.NEG[,1],-COEF.POS[,1], rep(0, p-1))
              } else {
              RR<-rbind(Ident, -Ident)
              rr<-c(COEF.NEG[,1],-COEF.POS[,1])
              }
      o.start<-rq.fit(x=B,y=y,tau=start.tau,method="fnc",R=RR,r=rr)
        }
#-------------------------------
      all.COEF<-cbind(COEF.NEG[,n.neg.taus:1], o.start$coef, COEF.POS)
      colnames(all.COEF)<-paste(taus)
      all.FIT<-cbind(FIT.NEG[,n.neg.taus:1,drop=FALSE], o.start$fitted.values[1:n], FIT.POS)
      colnames(all.FIT)<-paste(taus)
      all.RES<-cbind(RES.NEG[,n.neg.taus:1,drop=FALSE], o.start$residuals[1:n], RES.POS)
      colnames(all.RES)<-paste(taus)
      all.df<- c(df.neg.tau[n.neg.taus:1], sum(round(o.start$residuals[1:n],2)==0), df.pos.tau)
      all.rho<-c(rho.neg.tau[n.neg.taus:1], sum(Rho(o.start$residuals[1:n], start.tau)) , rho.pos.tau)
      r<-list(coefficients=all.COEF,B=B, df=all.df, rho=all.rho, fitted.values=all.FIT, residuals=all.RES)
      } ##end se length(taus)>1
      if(plott>0){
          if(plott==1) {matlines(x, B[1:n,]%*%all.COEF ,lwd=2,...)
              } else {plot(x,y[1:n]); matpoints(x, B[1:n,]%*%all.COEF ,lwd=2, type="l",...)}
          }
      return(r)
      }
