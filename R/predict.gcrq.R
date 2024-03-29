predict.gcrq <-
function(object, newdata, se.fit=FALSE, transf=NULL, xreg, type=c("sandw","boot"), ...){
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
    } else {
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
  n<-length(object$y)
  b<-as.matrix(object$coefficients)  #corretto in 0.3-1
  if(missing(xreg)){
    if(missing(newdata)) {
      fit<- object$fitted.values
      xreg <- object$x[1:n,]
    } else {
      if(is.null(object$info.smooth)){ #if the model is linear
        nomiVariabiliMod <- as.vector(unlist(sapply(object$BB, function(.x) attr(.x,"coef.names"))))
        id.var<-match(nomiVariabiliMod, names(newdata))
        #Qui controlla i nomi... come fare a controllare se factor?
        if(any(is.na(id.var))) stop("`newdata' does not include all the covariates in the model")
        xreg <-stats::model.matrix(as.formula(paste("~0+",paste(colnames(newdata),collapse = "+"))), data= newdata)
        nomiCoef <- rownames(b)
        if("(Intercept)" %in% nomiCoef ) xreg<-cbind("(Intercept)"=1, xreg)
        xreg <- xreg[,nomiCoef]
        fit<-drop(xreg%*%b) 
      } else {
        
        #browser()
        
        nomiCoef<-rownames(b)
        idVeriSmooth<-sapply(object$BB, function(.x) !is.null(attr(.x,"smoothName1")))
        nome.smooth <- sapply(object$BB, function(.x) attr(.x, "smoothName")) ##nomi variabili smooth (quelle vere.. non con ps()...) usa "smoothName1" 
        nome.smooth<- nome.smooth[idVeriSmooth] #solo i veri smooth, togli i termini lineari...
        n.smooth <-length(nome.smooth)
        
        nomiCoefPen <- sapply(object$BB, function(.x) attr(.x,"coef.names"))
        nomiCoefPen <- as.vector(unlist(nomiCoefPen[idVeriSmooth]))
        nomiCoefUnpen <- setdiff(nomiCoef, nomiCoefPen )
        id.coef.smooth <-match(nomiCoefPen, nomiCoef)
        b.smooth<-b[id.coef.smooth, ]
        
        info.smooth<-object$info.smooth #estrai 

        nomiVarModelloRispo<-all.vars(object$call$formula) #include la risposta che puo' includere piu nome se e scritta come y[g==1]
        nRispo<-length(all.vars(update.formula(object$call$formula,.~1)))
        nomiVarModello <-nomiVarModelloRispo[-(1:length(nRispo))]
        #"nomiVarModello" puo includere anche altri nomi che non sono variabili, cioe argomenti di ps(). Ad es., se la call e ps(x,center=T)
        #ci sara il nome T.. Quindi si deve eliminare..
        nomiID<-sapply(nomiVarModello, function(.x) grep(.x, nomiCoef),simplify =FALSE)
        nomiVarModello <-names(nomiID)[sapply(nomiID,length) !=0]
        id.var<-match(nomiVarModello, names(newdata))
        
        
        #browser()
        
        
        #Qui controlla i nomi... come fare a controllare se factor?
        if(any(is.na(id.var))) stop("`newdata' does not include all the covariates in the model")
        newdata <-newdata[,id.var,drop=FALSE] 
        id.var.smooth <- match(nome.smooth,names(newdata))
        x.new <-newdata[, id.var.smooth,drop=FALSE] #individua le variabili smooth
        newdata <-newdata[,-id.var.smooth, drop=FALSE] #togli le variabili smooth
        #costruire la base su x.new.. prendere il min max..
        B.new.list<- B.fixed.list<- vector("list", length=n.smooth)
        for(i in 1:n.smooth){
            m<-min(attr(object$BB[[i]], "covariate.35"))       #as.numeric(attr(object$BB[[1]], "covariate.n"))
            M<-max(attr(object$BB[[i]], "covariate.35"))
            B.new<-bspline(c(m, x.new[,i], M), ndx=info.smooth$ndx[i], deg=info.smooth$deg[i], knots=info.smooth$knots[[i]])
            if(!is.null(attr(object$BB[[i]],"name.fixed.params"))){
              nomi.okF <-attr(object$BB[[i]],"name.fixed.params")
              d<-length(nomi.okF)+1 #ricostruiamo l'ordine diff
              Di<-diff(diag(ncol(B.new)), diff=d)
              B.new <- B.new %*%t(solve(tcrossprod(Di),Di)) #sono ncol(B[[j]])-vDiff[j]
              B.fixed.list[[i]]<-poly(x.new[,i], degree=d-1, raw=TRUE)
              colnames(B.fixed.list[[i]])<- nomi.okF 
            } else {
              #modifica la base in funzione di drop, center, vc
              if(isTRUE(attr(object$BB[[i]],"drop"))) B.new<-B.new[,-1,drop=FALSE]
              if(isTRUE(attr(object$BB[[i]],"center"))){
                colmeansB<-attr(object$BB[[i]],"colmeansB")
                B.new0 <- B.new
                B.new <- sweep(B.new, 2, colmeansB)
              }
            }
            B.new<-B.new[-c(1,nrow(B.new)),,drop=FALSE] #rimuovi le righe relative al min e max aggiunte sopra!
    
            ## per VC models devi semplicemente moltiplicare B.new<- xvc*cbind(1,B.new)
            if(isTRUE(attr(object$BB[[i]],"vc"))) { 
              nomeVC <- attr(object$BB[[i]], "vcName")
              xVC <-newdata[, nomeVC ] 
              #if(is.factor(xVC)){ #
              if(!is.null(attr(object$BB[[i]], "vcLevels"))){
                xVC<-factor(xVC, levels=attr(object$BB[[i]], "vcLevels"))
                Mm<-stats::model.matrix(~0+xVC)
                colnames(Mm) <- attr(object$BB[[i]], "vcLevels")
                B.new<- Mm[,attr(object$BB[[i]], "vcCategory")]*B.new #cbind(1, B.new)
              } else {
                B.new<- xVC*B.new #cbind(1, B.new)
              }
            }
            colnames(B.new) <- attr(object$BB[[i]], "coef.names")
            B.new.list[[i]]<- B.new
        }
        #browser()

        B.new <- do.call(cbind, B.new.list)
        B.fixed.new <- do.call(cbind, B.fixed.list)
        B.new <- cbind(B.new, B.fixed.new)
        if(ncol(newdata)>0) {
          if("(Intercept)" %in% nomiCoef ) {
            newdata<-stats::model.matrix(as.formula(paste("~1+",paste(colnames(newdata),collapse = "+"))), #perche' era paste("~0+",?? mha <=1.6-2
                    data= newdata, contrasts = object$contrasts, xlev= object$xlevels)
          } else {
            newdata<-stats::model.matrix(as.formula(paste("~0+",paste(colnames(newdata),collapse = "+"))), #perche' era paste("~0+",?? mha <=1.6-2
                                  data= newdata, contrasts = object$contrasts, xlev= object$xlevels)
          }
          xreg <-as.matrix(cbind(newdata,B.new))
        } else {
          xreg<-B.new
        }
        #browser()
        #se non ci sono altre esplicative (e quindi ncol(newdata)=0) devi comunque aggiungere l'intercetta
        if("(Intercept)" %in% nomiCoef && !"(Intercept)"%in%colnames(xreg)) xreg<-cbind("(Intercept)"=1, xreg)
        #if(!"(Intercept)" %in% nomiCoef ) xreg<-xreg[,-match("(Intercept)", colnames(xreg)),drop=FALSE] 
        
        xreg <- xreg[,nomiCoef]
        fit<-drop(xreg%*%b) #mettere b[nomiCoef,]? Non serve perche' nomiCoef<-rownames(b)
      } #end if there are smooths
    } #end if(missing(newdata)) else.. 
  } else { #se c'e' xreg
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
      V<-vcov.gcrq(object, type=type)
      se<-sapply(V, function(x)sqrt(rowSums((xreg %*% x * xreg))))
      if(!is.null(f.transf.inv)) se <- se*abs(apply(fit, 2, f.transf.inv, deriv=1))
      fit<-list(fit=fit, se.fit=drop(se))
  }
  return(fit)
}
