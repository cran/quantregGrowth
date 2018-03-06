vcov.gcrq <- function(object, term, ...){
        if(!("boot.coef" %in% names(object))) {
              object <- update(object, n.boot=100)
              warning("The fit misses boot replicates.. 100 runs have been performed now")
                         }
        #if(is.null(object$boot.coef)) stop(" 'vcov.gcrq' works only with boot")
        coef.boot<-object$boot.coef
        n.boot<-dim(coef.boot)[3]
        n.tau<- dim(coef.boot)[2]
        n.coef<-dim(coef.boot)[1]
#come calcolare la vcov complessiva???
#non saprei. Comunque fissato il percentile, puoi ottenere la vcov
#Una naive puoi ottenerla calcolando per ogni percentile la formula sandwich.. (vedi Koenker bandaid)
        VCOV<-NULL
        for(j in 1:n.tau) {
            m<-var(t(coef.boot[,j,]))
            colnames(m)<-rownames(m)<-rownames(as.matrix(object$coefficients))
            VCOV[[length(VCOV)+1]]<-m
            }
        if(!missing(term)){
            #if(!term %in% names(object$BB)) stop(" 'term' is not a smooth variable")
            if(!term %in% all.vars(formula(object))[-1]) stop(" 'term' is not a variable in the model")
            f.smooth<-function(x){
                id<-startsWith(colnames(x),paste(term, "ps.",sep="."))
                r<-x[id,id,drop=TRUE]
                r
                }
            f.lin<-function(x){
                id<-grep(term,colnames(x))
                r<-x[id,id,drop=TRUE]
                r
                }
            f.ok<- if(term %in% names(object$BB)) f.smooth else f.lin

            VCOV<- lapply(VCOV, f.ok)
            }
        names(VCOV)<-paste(object$tau) #colnames(object$coefficients)
        return(VCOV)
        }
        
        
        