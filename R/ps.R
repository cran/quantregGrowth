ps <-
function(..., lambda=-1, d=3, by=NULL, ndx=NULL, deg=3, knots=NULL, 
    monotone=0, concave=0, var.pen=NULL, pen.matrix=NULL, dropc=TRUE, center=TRUE,
    K=2, decompose=FALSE, constr.fit=TRUE, shared.pen=FALSE){
#arguments in ... such as 'a=3' are ignored
#dropc, if TRUE the first column of the basis is dropped..
#decompose: if TRUE the decomposition of the B-spline is employed. Notice the new basis becomes
#     [x, ..., x^{d-1}, Z] - the intercept is *not* included in the basis (but it is included in the model)
#     decompose is incompatible with monotonocity restrictions. Also if 'decompose=TRUE', dropc is set to FALSE.
#ndx: if NULL the empirical rule of Ruppert (2001) is used, min(n/4,40)
#knots: eventuale vettore di nodi; se fornito ndx e' ignorato
#d: the diff order (d=0 means a ridge penalty matrix) #prima era 'pdiff'
#se lambda<0 allora viene stimato, altrimenti deve essere un valore numerico
#   se lambda e' vettore, deve essere di length(lambda)=lenght(tau) causing a different amount of smoothing in 
#   each curve.
#var.pen: una stringa del tipo "1:k" per varying penalty
#monotone: 0: unconstrained; +1: non-descreasing; -1= non-increasing (NB sign(T)=1 and sign(F)=0)
#K fattore che regola la selezione dello spar (se lambda=-1)
#ridge: if TRUE, d is set to 0 and ndx is set to length(unique(x)).
#     ridge overwrites everything
#NB x NON deve essere un factor, anche se ridge=TRUE (sebbene si possa usare as.numeric(as.character(x)) )
#   e puo' essere una matrice??
#pen.matrix l'eventuale matrice A t.c. A'A e' la matrice di penalizzaz. Se fornita, e' questa che viene utilizzata! 
    #browser()
    #nomi<-sapply(as.list(substitute(list(...)))[-1], function(xx) as.character(xx)) #da problemi se il termine e' log(x)
    
    #nomi<-sapply(as.list(substitute(list(...)))[-1], function(xx) all.vars(xx)) #funziona anche con log(x).. FALSO al 14/10!!!
    nomi<-as.list(substitute(list(...)))[-1]
    
    vars<-list(...)

    id<-""==names(nomi) #e' TRUE solo se *non* e' relativa ad "altri" argomenti messi in ps(), e quindi indica realmente una variabile..
    if(length(id)<=0) id<-rep(TRUE, length(nomi))
    vars<-vars[id]
    nomi<-nomi[id]    
    names(vars)<-nomi
    

    if(is.factor(vars[[1]]) || is.matrix(vars[[1]]) || length(nomi)>1) stop("ps() accepts single numerical terms only")
    #if(length(nomi)>1) stop(" ps() accepts single terms only") #: numeric/factor vector or matrix")
    
    p<-1
    if(is.matrix(vars[[1]])) p<- ncol(vars[[1]])
    r<-matrix(unlist(vars), ncol=p, byrow = FALSE, dimnames = list(NULL, colnames(vars[[1]])))    
    
    #da mgcv::s
    #vars <- as.list(substitute(list(...)))[-1]
    #p<-length(vars)
    #by.var <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
    #term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
    #browser()
    r<- if(!is.null(by)) cbind(r,by) else cbind(r) #cbind(r,1)
    #r<- apply(r,2, function(.x)as.numeric(as.character(as.factor(.x)))) #non funziona se la prima colonna non ha numeri..
    #r<- apply(r,2, function(.x)as.numeric(as.factor(.x)))
    attr(r,"penMatrix")<-pen.matrix
    attr(r,"nodi")<-knots
    attr(r,"ndx")<-ndx
    attr(r,"deg")<-deg
    attr(r,"pdiff")<-d
    attr(r,"monot")<-monotone #isTRUE(monotone)
    attr(r,"conc")<-concave
    attr(r,"constr.fit")<-constr.fit
    attr(r,"lambda")<-lambda
    attr(r,"nomeX")<- (nomi[[1]])#deparse(substitute(x))
    attr(r,"var.pen")<-var.pen
    attr(r,"K")<-K
    attr(r,"ridge") <-if(is.matrix(vars[[1]]) || is.factor(vars[[1]])) TRUE else FALSE
    attr(r,"nomeBy")<-deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
    attr(r,"levelsBy")<-levels(by)
    attr(r,"dimSmooth")<-p
    attr(r,"decom")<-decompose
    attr(r,"shared.pen")<-shared.pen
    if(monotone!=0 && decompose) stop("'decom=TRUE' is incompatible with monotonicity restrictions") #se i vincoli sono sui fitted e' OK..
    if(concave!=0 && decompose) stop("'decom=TRUE' is incompatible with concavity restrictions")
    if(decompose) dropc<-FALSE
    
    if(!is.null(by)) colnames(r)<-c(nomi, "by") else nomi
    #colnames(r)<-c(nomi, "by")
    #if(!is.null(by)) colnames(r)<-c(deparse(substitute(x)), deparse(substitute(by)))
    attr(r,"dropc")<-dropc
    attr(r,"center")<-center
    if(dropc && decompose) stop("'decom=TRUE' is incompatible with 'dropc=TRUE' ")
	#browser()
	r
}

