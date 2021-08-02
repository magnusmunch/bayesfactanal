# x=xy[1:(n + m[s]), ]; factors=est.factors; groups=c(groups, 0);
# hyper=hyper3; method="eb"; rescale=TRUE; constrained=TRUE
# start=NULL
# control=list(maxit=1000, nsamples=1000, nchains=1,
#              ncores=1,
#              epsilon=sqrt(.Machine$double.eps),
#              trace=FALSE)
# conditionally conjugate Bayesian factor analysis
bayesfactanal <- function(x, factors, hyper=NULL, groups=NULL, 
                          method="vb", start=NULL, rescale=FALSE,
                          constrained=NULL,
                          control=list(maxit=1000, nsamples=1000, nchains=1,
                                       ncores=1,
                                       epsilon=sqrt(.Machine$double.eps), 
                                       trace=FALSE)) {
  
  # set default hyperparameters
  if(is.null(hyper)) {
    hyper <- list(kappa=rep(0.01, ncol(x)), nu=rep(0.01, ncol(x)), 
                  gamma=rep(1, ncol(x)))
  }
  
  # estimate posterior
  if(method=="vb") {
    res <- vbfactanal(x=x, factors=factors, hyper=hyper, start=start, 
                      rescale=rescale, 
                      control=control[c("maxit", "epsilon", "trace")])
  } else if(method=="mcmc") {
    res <- mcmcfactanal(x=x, factors=factors, hyper=hyper, start=start, 
                        control=control[c("nsamples", "nchains", "ncores", 
                                          "trace")])
  } else if(method=="eb") {
    res <- ebfactanal(x=x, factors=factors, hyper=hyper, start=start, 
                      groups=groups, rescale=rescale, constrained=constrained,
                      control=control[c("maxit", "epsilon", "trace")])
  }
  
  return(res)
}