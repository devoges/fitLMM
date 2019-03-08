##-- found online at
##-- https://rstudio-pubs-static.s3.amazonaws.com/305732_7a0cfc5b535c41c8b56ffdc2e322a51d.html

#' @importFrom Matrix Matrix bdiag
#' @importFrom stats rmultinom
#' @importFrom methods as
myReTrms<-function(ListZ,ListVar=NULL){
  reTrms<-list()
  reTrms$Zt    <- Matrix::Matrix(t(Reduce('cbind',ListZ)),sparse=TRUE)
  reTrms$Zt <- as(reTrms$Zt, "dgCMatrix")
  reTrms$theta <- rep(1,length(ListZ))  # Initial Value of the covariance parameters
  reTrms$Lind  <- rep(1:length(ListZ),unlist(lapply(ListZ,ncol))) #an integer vector of indices determining the mapping of the elements of the theta vector to the "x" slot of Lambdat
  reTrms$Gp    <- as.integer(unname(cumsum(c(0,unlist(lapply(ListZ,ncol))))))
  reTrms$lower   <- rep(0,length(ListZ)) # lower bounds on the covariance parameters
  if (is.null(ListVar)) {
    reTrms$Lambdat <- Matrix(diag(rep(1,sum(unlist(lapply(ListZ,ncol))))),sparse=TRUE)
  }  else {
    reTrms$Lambdat <- bdiag(lapply(ListVar,chol))
  }
  reTrms$Ztlist <- lapply(ListZ,function(Z) as(Matrix(t(Z),sparse=T), "dgCMatrix"))
  reTrms$cnms   <-  as.list(names(ListZ)) ; names(reTrms$cnms)<- names(ListZ)
  # Flist is Not very clean (to say the least... )
  reTrms$flist <- lapply(ListZ,function(Z) {flist.factor<- as.factor(colnames(Z)[apply(Z,1,function(x) which(rmultinom(n=1,size=1,prob =abs(x)+0.1)==1) )]);
  levels(flist.factor)<-colnames(Z); return(flist.factor)}) #NULL # list of grouping factors used in the random-effects terms (used for computing initial variance ??)
  return(reTrms)
}

#' @importFrom stats sigma
myranef<-function(model,condVar=TRUE){
  re.cond.mode<-tapply(model@u,mymod@pp$Lind,function(x) x)
  names(re.cond.mode)<- names(model@cnms)

  if (condVar) {
    Zt<-model@pp$Zt
    D  <- sigma(model)* t(model@pp$Lambdat) %*% model@pp$Lambdat
    Sigma<- t(Zt)%*%  D %*%Zt + sigma(model)*diag(rep(1,ncol(Zt)))
    var.cond <- D - Zt %*%solve(Sigma) %*% t(Zt)
    var.cond <- diag(var.cond)
    var.cond.mode <- tapply(var.cond,mymod@pp$Lind,function(x) x)
  }
  for (i in 1:length(re.cond.mode)) {
    re.cond.mode[[i]]<-data.frame(re.cond.mode[i])
    names(re.cond.mode[[i]])<-names(re.cond.mode[i])
    row.names(re.cond.mode[[i]]) <- levels(model@flist[[i]])

    if (condVar) attr(re.cond.mode[[i]],"postVar")=array(var.cond.mode[[i]],c(1,1,nrow(re.cond.mode[[i]])))
  }
  attr(re.cond.mode,"class") <- "ranef.mer"
  re.cond.mode
}


#' @importFrom stats model.frame setNames
#' @importFrom lme4 mkLmerDevfun optimizeLmer mkMerMod getME VarCorr
#' @importFrom dplyr case_when mutate
#' @importFrom magrittr %>%
#' @export
fit_lmer <- function(Y,
                     X,
                     Zlist,
                     REML = TRUE) {
  fr <- model.frame(Y ~ .,
                    data.frame(Response = Y,
                               X))
  reTrms <- myReTrms(Zlist, NULL)
  reTrms$Zt <- as(reTrms$Zt, "dgCMatrix")
  devfun <- mkLmerDevfun(fr, X, reTrms, REML = REML)
  opt <- optimizeLmer(devfun)
  model <- mkMerMod(environment(devfun), opt, reTrms, fr = fr)

  beta <- getME(model, name = "fixef")
  foo <- as.data.frame(VarCorr(model))[, c("var1", "vcov")] %>%
    mutate(var1 = case_when(is.na(var1) ~ "Error",
                            TRUE ~ var1))
  vcs <- setNames(foo$vcov, foo$var1)
  vcs <- vcs[order(names(vcs))]

  list("type" = "lm4",
       "model" = model,
       "beta" = beta,
       "vcs" = vcs)
}
