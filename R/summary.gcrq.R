summary.gcrq <-
function(object, type=c("sandw","boot"), digits = max(3, getOption("digits") - 3), signif.stars =getOption("show.signif.stars"), ...){
      edfS<-object$edf.j
      #browser()
      if(rownames(edfS)[1]=="Xlin") edfS<-edfS[-1,,drop=FALSE]
      #rownames(edfS)<-paste("ps(", names(object$BB), ")",sep="")
      colnames(edfS)<-rep("edf", ncol(edfS))
      #cat("\n***Noncrossing regression quantiles via P-splines***\n")
      cat("Call:\n")
      print(object$call)

      n<-nrow(as.matrix(object$fitted.values))
      p<-nrow(as.matrix(object$coefficients ))
      n.tau<-ncol(as.matrix(object$coefficients))
      sic<- sum(log(object$rho/n)) +log(n)*sum(object$edf.j)/(2*n)
      
      list.vcov<-vcov.gcrq(object, type=type)
      #if(!is.null(object$boot.coef)) list.vcov<-vcov.gcrq(object)
      for(j in 1:n.tau){
          est<-as.matrix(object$coefficients)[,j]
          se<-sqrt(diag(list.vcov[[j]])) 
          ris<-cbind(Est=est, StErr=se, "|z|"=round(abs(est)/se,2), 
                     "p-value"=pchisq((est/se)^2,df=1,lower.tail = FALSE) )
          rownames(ris)<-rownames(as.matrix(object$coefficients))
          nomi.p.spline<-unlist(unname(sapply(object$BB, function(.x) attr(.x,"coef.names"))))
          #nomi.p.spline<-as.vector(sapply(object$BB, function(.x) attr(.x,"coef.names")))
          nomi.param<-setdiff(rownames(ris), nomi.p.spline )
          ris<-ris[nomi.param,,drop=FALSE]
          cat("\n--------  Percentile:", object$taus[j], "  check function: ", round(object$rho[j], digits-1) , "-----\n")
          if(object$pLin>0){
             cat("\nparametric terms:\n")
             printCoefmat(ris, digits = digits, signif.stars = signif.stars)
         }
          cat("\nsmooth terms:\n")
          printCoefmat(edfS[,j,drop=FALSE], digits = digits)
          #cat("====================\n")
          #cat("\n")
      }          
      #cat("\nNo. of obs:", n, "  Check function =", round(sum(object$rho),digits-1), "  SIC =", round(sic,digits-1),"\n")
      cat("\n=====================\nNo. of obs:", n, "  Check function =", round(sum(object$rho),digits-1), 
              " SIC =", round(sic,digits), " ( on edf =",round(sum(object$edf.j),3),")" ,"\n")
      cat("No. of params:", p,"(for each curve);", p*n.tau,"(total)\n")       
   }
       