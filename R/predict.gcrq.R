predict.gcrq <-
function(object, newdata, se.fit=FALSE, transf=NULL, xreg, ...){
#xreg o newdata have to be provided. Nessun controllo (neanche sull'ordine dei coef)
#if xreg is provided, newdata is ignored
#
#Attenzione: puo' dare problemi (in realta' non funziona) se la formula contiene "factor" o "poly"
#model.matrix(as.formula(m1$call$formula), data=newdata)
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
    b<-as.matrix(object$coefficients)  #corretto in 0.3-1

    if(missing(xreg)){
      if(missing(newdata)) {#stop("please, provide `newdata' or 'xreg'")
         fit<- object$fitted.values
            } else {
      nomiCoef<-rownames(b)
      info.smooth<-object$info.smooth #estrai 
      p<- info.smooth$deg + info.smooth$ndx
      nome.smooth<-names(object$BB) #nomi variabili smooth
      n.smooth<-length(nome.smooth)
      #se p NON e' scalare:
      #if(length(nome.smooth)>1){
      nomiCoefPen.List<-mapply(function(.x, .y) paste(.x, ".ps." ,1:.y,sep=""), nome.smooth, p, SIMPLIFY=FALSE)
      #}      
      nomiCoefPen <- unlist(nomiCoefPen.List)      
      #se p e' scalare: nomiCoefPen<-paste(nome.smooth, ".ps." ,1:p,sep="")
      id.coef.smooth<-match(nomiCoefPen, nomiCoef)
      b.smooth<-b[id.coef.smooth, ]
      nomiCoefUnpen<-nomiCoef[-id.coef.smooth]
      nomiVarModello<-all.vars(object$call$formula)[-1]
      id.var<-match(nomiVarModello, names(newdata))
      if(any(is.na(id.var))) stop("`newdata' does not include all the covariates in the model")
      newdata<-newdata[,id.var,drop=FALSE] 
      x.new<-newdata[,match(nome.smooth,names(newdata)),drop=FALSE] #individua le variabili smooth
      newdata<-newdata[,-match(nome.smooth,names(newdata)),drop=FALSE] #togli le variabili smooth
      #costruire la base su x.new.. prendere il min max..
      
      B.new.list<-vector("list", length=n.smooth)
      for(i in 1:n.smooth){
          m<-min(attr(object$BB[[i]], "covariate.35"))       #as.numeric(attr(object$BB[[1]], "covariate.n"))
          M<-max(attr(object$BB[[i]], "covariate.35"))
          B.new<-bspline(c(m, x.new[,i], M), ndx=info.smooth$ndx[i], deg=info.smooth$deg[i])
          B.new<-B.new[-c(1,nrow(B.new)),,drop=FALSE] #rimuovi le righe relative al min e max aggiunte sopra!
          B.new.list[[i]]<- B.new           
            }
      B.new<- do.call(cbind, B.new.list)
      xreg <-as.matrix(cbind(newdata,B.new))
      colnames(xreg)<-c(nomiCoefUnpen, nomiCoefPen)
      if("(Intercept)" %in% nomiCoef ) xreg<-cbind("(Intercept)"=1, xreg)
      fit<-drop(xreg[,nomiCoef]%*%b)
          }
      } else {
      if(!missing(newdata)) warning("`newdata' ignored when 'xreg' is provided")      
      fit<-drop(xreg%*%b)
      }
      
      if(is.null(transf)){
              f.transf.inv<-attr(object$fitted.values, "transf.inv")
              } else {
                if(!is.character(transf)) stop(" 'transf' should be NULL or a character")
                f.transf.inv<- eval(parse(text=paste("function(y){", transf, "}")))
              } 
      if(!is.null(f.transf.inv)) fit<-apply(fit, 2, f.transf.inv) #fit <- f.transf.inv(fit) #l$y <- eval(parse(text=transf), list(y=l$y))               
      
      if(se.fit){
             V<-vcov.gcrq(object)
             se<-sapply(V, function(x)sqrt(rowSums((xreg %*% x * xreg))))
             if(!is.null(f.transf.inv)) se <- se*abs(apply(fit, 2, f.transf.inv, deriv=1))
             fit<-list(fit=fit, se=se)
             }
      return(fit)
      }
