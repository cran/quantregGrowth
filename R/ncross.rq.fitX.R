ncross.rq.fitX <- function(y, X = NULL, taus, adjX.constr=TRUE, lambda.ridge = 0, eps = 1e-04, ...) {
#                           sgn.constr = 1, adjX.constr=TRUE, ...) {
  # Stima dei non-crossing rq con X lineari (la X dovrebbe avere una colonna di 1, se richiesta..) 
  #--------------------------------------------------------
      adj.middle = FALSE
      Rho <- function(u, tau) u * (tau - (u < 0))
      #-------------------------------------
      # require(quantreg)
      if (length(taus) <= 1) {
        o <- rq.fit(x = X, y = y, tau = taus, ...)
        o$rho <- sum(Rho(o$residuals[1:length(y)], taus))
        return(o)
      }
      n <- length(y)
      B <- X
      p <- ncol(B)
      taus <- sort(taus)
      
      id.start.tau <- which.min(abs(taus - 0.5))
      start.tau <- taus[id.start.tau]
      
      pos.taus <- taus[(taus - start.tau) > 0]
      neg.taus <- taus[(taus - start.tau) < 0]
      n.pos.taus <- length(pos.taus)
      n.neg.taus <- length(neg.taus)
      
      if (lambda.ridge > 0) {
        DD <- diag(p)
        B <- rbind(B, lambda.ridge * DD)
        y <- c(y, rep(0, nrow(DD)))
      }

      if (length(taus) <= 1) {
        o.start <- rq.fit(x = B, y = y, tau = start.tau)
        all.COEF <- o.start$coef
        # colnames(all.COEF)<-paste(taus)
        all.df <- sum(round(o.start$residuals[1:n], 2) == 0)
        all.rho <- sum(Rho(o.start$residuals[1:n], start.tau))
        r <- list(coefficients = all.COEF, B = B, df = all.df, rho = all.rho, 
                  fitted.values = o.start$fitted.values[1:n], residuals = o.start$residuals[1:n])
      } else { # se length(taus)>1
        #browser()
        DF.NEG <- DF.POS<- NULL
        Ident <- diag(p)
        
        colnamesB <- colnames(B)
        if("(Intercept)" %in% colnames(B)) {
          if(!ncol(B)>=2) stop("At leat one covariate")
          is.inter <- TRUE
          minX <- apply(B[,-1,drop=FALSE] ,2, min)
          names(minX)<-colnames(B)[-1]
          B<- cbind(B[,1], apply(B[,-1,drop=FALSE], 2, function(.x) .x- min(.x)))
        } else {
          if(!ncol(B)>=1) stop("At leat one covariate")
          is.inter <- FALSE
          minX <- apply(B ,2, min)
          names(minX)<-colnames(B)
          B<- apply(B, 2, function(.x) .x- min(.x))
        }
        colnames(B)<-colnamesB
        all.max <- apply(B, 2, max)
        
        

        
        
        # signB <- apply(B, 2, function(.x){prod(sign(range(.x)))} )
        # if(any(signB<0)) warning("covariate(s) taking pos AND neg values: noncrossing is not guaranteed.", call. = FALSE)
        # 
        # 
        # if (is.null(sgn.constr)) {
        #   # i segni per i vincoli..
        #   yres <- lm.fit(y = y, x = B)$residuals^2
        #   if("(Intercept)" %in% colnames(B)){
        #     sgn.constr <- c(1, sign(lm.fit(y = yres, x = B)$coefficients[-1]))
        #   } else {
        #     sgn.constr <- sign(lm.fit(y = yres, x = B)$coefficients)
        #   }
        # } else {
        #   sgn.constr<-rep(sgn.constr, l=ncol(B))
        # }
        # #traslare la variabile?
        # if(adjX.constr && ("(Intercept)" %in% colnames(B))){
        #   const<-c(NA, rep(0, ncol(B)-1))
        #   for(j in 2:ncol(B)){
        #     if((sgn.constr[j]<0) && (min(B[,j])>=0)) {const[j]<- max(B[,j]); B[,j]<- B[,j] - max(B[,j])}  
        #     if((sgn.constr[j]>0) && (max(B[,j])<=0)) {const[j]<- min(B[,j]); B[,j]<- B[,j] - min(B[,j])}
        #   }
        # } else {
        #   const<-rep(0, ncol(B))
        # }
        # names(const)<-colnames(B)

        o.start <- rq.fit(x = B, y = y, tau = start.tau)        
        
        COEF.POS <- COEF.NEG <- FIT.POS <- FIT.NEG <- RES.POS <- RES.NEG <- NULL
        df.pos.tau <- df.neg.tau <- rho.pos.tau <- rho.neg.tau <- NULL
        if (n.pos.taus > 0) {
          rho.pos.tau <- df.pos.tau <- vector(length = n.pos.taus)
          COEF.POS <- matrix(, ncol(B), n.pos.taus)
          colnames(COEF.POS) <- paste(pos.taus)
          b.start <- o.start$coef

          if(is.inter){
            RR <- diag(all.max)
            RR[,1]<-1
            rr <- c(b.start[1], b.start[1] + b.start[-1]*all.max[-1]) + eps
          } else {
            RR<-Ident
            rr <- b.start + eps
          }
          #diag(RR) <- diag(RR) #* sgn.constr
          #rr <- rr * sgn.constr
          
          # DF per NONcrossing
          DF.POS <- matrix(,nrow(RR),n.pos.taus)
          colnames(DF.POS) <- paste(pos.taus)

          FIT.POS <- RES.POS <- matrix(, n, n.pos.taus)
          for (i in 1:n.pos.taus) {
            o <- rq.fit(x = B, y = y, tau = pos.taus[i], method = "fnc", R = RR, r = rr)
            FIT.POS[, i] <- o$fitted.values[1:n]
            RES.POS[, i] <- o$residuals[1:n]
            # estrai la f. obiettivo
            df.pos.tau[i] <- sum(abs(o$residuals[1:n]) <= 1e-06)  #length(o$coef)
            rho.pos.tau[i] <- sum(Rho(o$residuals[1:n], pos.taus[i]))
            b.start <- o$coef
            COEF.POS[, i] <- b.start
            DF.POS[,i] <- 1*(drop(RR%*%b.start-rr)<=1e-8)
            if(is.inter) {
              rr <- c(b.start[1], b.start[1] + b.start[-1]*all.max[-1]) + eps
            } else {
              rr <- b.start + eps
            }
            #rr <- b.start + eps
            #rr <- rr * sgn.constr
          }  #end for
        }  #end if(n.pos.taus>0)
        
        if (n.neg.taus > 0) {
          rho.neg.tau <- df.neg.tau <- vector(length = n.pos.taus)
          COEF.NEG <- matrix(, ncol(B), n.neg.taus)
          colnames(COEF.NEG) <- paste(neg.taus)
          b.start <- o.start$coef
          neg.taus <- sort(neg.taus, TRUE)
          
          if(is.inter){
            RR <- -diag(all.max)
            RR[,1]<- -1
            rr <- -(c(b.start[1], b.start[1] + b.start[-1]*all.max[-1]) - eps)
          } else {
            RR<- -Ident
            rr <- -b.start + eps
          }
          
          # RR <- -Ident
          # rr <- -b.start + eps
          # diag(RR) <- diag(RR) * sgn.constr
          # rr <- rr * sgn.constr
          #DF per NONcrossing
          DF.NEG<- matrix(,nrow(RR),n.neg.taus)
          colnames(DF.NEG)<-paste(neg.taus)

          FIT.NEG <- RES.NEG <- matrix(, n, n.neg.taus)
          for (i in 1:n.neg.taus) {
            o <- rq.fit(x = B, y = y, tau = neg.taus[i], method = "fnc", 
                        R = RR, r = rr)
            FIT.NEG[, i] <- o$fitted.values[1:n]
            RES.NEG[, i] <- o$residuals[1:n]
            df.neg.tau[i] <- sum(abs(o$residuals[1:n]) <= 1e-06)  #length(o$coef)
            rho.neg.tau[i] <- sum(Rho(o$residuals[1:n], neg.taus[i]))
            b.start <- o$coef
            COEF.NEG[, i] <- b.start
            DF.NEG[,i] <- 1*(drop(RR%*%b.start-rr)<=1e-8)
            if(is.inter) {
              rr <- -(c(b.start[1], b.start[1] + b.start[-1]*all.max[-1]) - eps)
            } else {
              rr <- -b.start + eps
            }
            
            #rr <- -b.start + eps
            #rr <- rr * sgn.constr
          }  #end for
        }  #end if(n.neg.taus>0)
        #------------------------------
        monotone <- FALSE
        R = NULL
        if (adj.middle) {
          if (monotone) {
            RR <- rbind(Ident, -Ident, R)
            rr <- c(COEF.NEG[, 1], -COEF.POS[, 1], rep(0, p - 1))
          } else {
            RR <- rbind(Ident, -Ident)
            rr <- c(COEF.NEG[, 1], -COEF.POS[, 1])
          }
          o.start <- rq.fit(x = B, y = y, tau = start.tau, method = "fnc", 
                            R = RR, r = rr)
        }
        #-------------------------------
        all.COEF <- cbind(COEF.NEG[, n.neg.taus:1, drop = FALSE], o.start$coef, 
                          COEF.POS)
        colnames(all.COEF) <- paste(taus)
        all.FIT <- cbind(FIT.NEG[, n.neg.taus:1, drop = FALSE], o.start$fitted.values[1:n], 
                         FIT.POS)
        colnames(all.FIT) <- paste(taus)
        all.RES <- cbind(RES.NEG[, n.neg.taus:1, drop = FALSE], o.start$residuals[1:n], 
                         RES.POS)
        colnames(all.RES) <- paste(taus)
        all.df <- c(df.neg.tau[n.neg.taus:1], sum(abs(o.start$residuals[1:n]) <= 1e-06), df.pos.tau)
        all.rho <- c(rho.neg.tau[n.neg.taus:1], sum(Rho(o.start$residuals[1:n], start.tau)), rho.pos.tau)
        r <- list(coefficients = all.COEF, x = B, df = all.df, rho = all.rho, 
                  fitted.values = all.FIT, residuals = all.RES, DF.NEG=DF.NEG, DF.POS=DF.POS, minX=minX)
      }
      id.coef <- 1:ncol(X)
      attr(id.coef, "nomi") <- colnames(X)
      r$id.coef <- id.coef
      return(r)
}
