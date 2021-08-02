# y=as.matrix(y); factors=est.factors2;
# hyper=hyper1;
# groups=list(x=part$corr, y=0); method="eb";
# start=NULL; rescale=TRUE;
# constrained=list(x=TRUE, y=FALSE);
# control=list(maxit=1000, epsilon=1e-4,
#              trace=TRUE)
# y=as.matrix(y); factors=est.factors2;
# hyper=hyper1;
# groups=NULL; method="vb";
# start=NULL; rescale=TRUE;
# control=list(maxit=1000, epsilon=1e-4,
#              trace=TRUE)
# conditionally conjugate logistic Bayesian factor analysis
logbayesfactanal <- function(x, y, factors, hyper=NULL, groups=NULL, 
                             method="vb", start=NULL, rescale=FALSE,
                             constrained=list(x=FALSE, y=FALSE),
                             control=list(maxit=1000, nsamples=1000, nchains=1,
                                          ncores=1,
                                          epsilon=sqrt(.Machine$double.eps),
                                          trace=TRUE)) {
  
  # set default hyperparameters
  if(is.null(hyper)) {
    hyper <- list(kappa=rep(0.01, ncol(x)), nu=rep(0.01, ncol(x)), 
                  gamma=list(x=rep(1, ncol(x)), 
                             y=rep(1, ifelse(is.null(ncol(y)), 1, ncol(y)))))
  }
  
  # estimate posterior
  if(method=="vb") {
    res <- logvbfactanal(x=x, y=y, factors=factors, hyper=hyper, start=start, 
                         rescale=rescale, 
                         control=control[c("maxit", "epsilon", "trace")])
  } else if(method=="mcmc") {
    res <- logmcmcfactanal(x=x, y=y, factors=factors, hyper=hyper, start=start, 
                           control=control[c("nsamples", "nchains", "ncores",
                                             "trace")])
  } else if(method=="eb") {
    res <- logebfactanal(x=x, y=y, factors=factors, hyper=hyper, start=start, 
                         groups=groups, rescale=rescale, 
                         constrained=constrained,
                         control=control[c("maxit", "epsilon", "trace")])
  }
  
  return(res)
}