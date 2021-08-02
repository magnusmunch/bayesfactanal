# factanal for non-convergence issues
.factanal_conv <- function(factors, covmat, n.obs=NA, rotation="none", 
                           start=NULL, 
                           control=list(nstart=10, opt=list(), lower=0.005)) {
  
  # try fitting until convergence or maximum number of fits reached
  error <- TRUE
  ntries <- 0
  while(error & ntries <= control$nstart) {
    ntries <- ntries + 1
    fit <- tryCatch(factanal(factors=factors, covmat=covmat, n.obs=n.obs, 
                             rotation=rotation, start=start, 
                             control=list(nstart=1, opt=control$opt, 
                                          lower=control$lower)),
                    error=function(e) {return(NULL)})
    start <- runif(nrow(covmat))
    error <- is.null(fit)
  }  
  if(is.null(fit)) {
    fit <- list(factors=factors)
  }
  return(c(fit, list(tries=ntries, conv=!error)))
}

# calculates expectation and variance of data
.expfun <- function(x, param, miss) {
  
  # fixed things
  p <- ncol(x)
  
  # calculate expectation and variance
  variance <- matrix(0, nrow=p, ncol=p)
  expectation <- x
  for(g in 1:length(miss$groupid)) {
    Btilde <- t(t(param$B)[-miss$covid[[g]], , drop=FALSE]/
                  param$psi[-miss$covid[[g]]])
    mat <- Btilde %*% t(param$B)[-miss$covid[[g]], , drop=FALSE]
    diag(mat) <- diag(mat) + 1
    mat <- t(param$B)[miss$covid[[g]], , drop=FALSE] %*% solve(mat) %*% 
      Btilde
    variance[miss$covid[[g]], miss$covid[[g]]] <- 
      variance[miss$covid[[g]], miss$covid[[g]]] + 
      (diag(length(param$psi[miss$covid[[g]]]))*
         param$psi[miss$covid[[g]]] + 
         t(param$B)[miss$covid[[g]], , drop=FALSE] %*% 
         t(t(param$B)[miss$covid[[g]], , drop=FALSE]) - 
         mat %*% t(param$B)[-miss$covid[[g]], , drop=FALSE] %*% 
         t(t(param$B)[miss$covid[[g]], , drop=FALSE]))*
      length(miss$groupid[[g]])
    expectation[miss$groupid[[g]], miss$covid[[g]]] <- 
      as.numeric(x[miss$groupid[[g]], -miss$covid[[g]], drop=FALSE] %*% t(mat))
  }
  return(list(expectation=expectation, variance=variance))
}

# EM update
.emupdate <- function(x, param.old, hyper, miss, rotation, start, control) {
  
  # fixed things
  n <- nrow(x)
  d <- nrow(param.old$B)
  
  # calculate expectation and variance of missing values
  expvar <- .expfun(x=x, param=param.old, miss=miss) 
  
  # calculate expected covariance
  S <- (t(expvar$expectation) %*% expvar$expectation + expvar$variance)/n
  S <- t(sqrt(1 - hyper$gamma)*t(sqrt(1 - hyper$gamma)*S)) +
    hyper$gamma*diag(hyper$Gamma)
  
  # maximuse expected log likelihood and return object
  fit <- .factanal_conv(factors=factors, covmat=S, n.obs=n, rotation=rotation,
                        start=start,
                        control=list(nstart=control$nstart, opt=control$opt,
                                     lower=control$lower))
  param <- list(B=t(fit$loadings), psi=unname(fit$uniquenesses),
                Upsilon=expvar$expectation,
                ell=unname(fit$criteria[1]), fit=fit)
  class(param$B) <- "matrix"; dimnames(param$B) <- NULL
  return(param)
}

# x=xy; factors=est.factors2; 
# gamma=c(rep(est.gamma2, ncol(xy) - 1), 0); 
# Gamma=rep(1, ncol(xy)); rotation="none"; 
# start=NULL;
# control=list(maxit=1000, nstart=10, lower=0.005, 
#              epsilon=sqrt(.Machine$double.eps))
# factor analysis with missing data
emfactanal <- function(x, factors, gamma=0, Gamma=NULL, rotation="none", 
                       start=NULL, 
                       control=list(maxit=100, nstart=10, lower=0.005,
                                    epsilon=sqrt(.Machine$double.eps))) {
  
  # locate missing values
  sampleid <- which(apply(is.na(x), 1, any))
  dups <- duplicated(is.na(x[sampleid, ]))
  miss <- list(groupid=lapply(1:length(which(!dups)), function(s) {
                 which(apply(is.na(x), 1, identical, is.na(x[sampleid, ])[
                   which(!dups)[s], ]))}),
               covid=apply(is.na(x[sampleid, , drop=FALSE])[
                 !dups, , drop=FALSE], 1, which))
  
  # if no penalty parameter supplied, we use identity target and CV
  if(is.null(Gamma)) {
    Gamma <- rep(1, ncol(x))
  }
  if(length(gamma)==1) {
    gamma <- rep(gamma, length.out=ncol(x))
    gamma[Gamma==0] <- 0
  }
  hyper <- list(gamma=gamma, Gamma=Gamma)
  
  # initial values
  param.old <- list(B=NULL, psi=NULL)
  fit <- .factanal_conv(factors=factors, 
                        covmat=t(sqrt(1 - hyper$gamma)*
                                   t(sqrt(1 - hyper$gamma)*cov(na.omit(x)))) +
                          hyper$gamma*diag(hyper$Gamma), n.obs=nrow(na.omit(x)), 
                        rotation=rotation, start=start,
                        control=list(nstart=control$nstart, opt=list(),
                                     lower=control$lower))
  param.old <- c(list(B=t(fit$loadings), psi=unname(fit$uniquenesses),
                      Upsilon=na.omit(x), ell=as.numeric(fit$criteria[1])), fit)
  class(param.old$B) <- "matrix"; dimnames(param.old$B) <- NULL
  
  # only iterate if there are missing values
  ell <- as.numeric(fit$criteria[1])
  conv <- TRUE
  iter <- 0
  if(length(miss$covid)!=0) {
    
    # iterations
    check <- FALSE
    while(!check) {
      iter <- iter + 1
      param <- .emupdate(x=x, param.old=param.old, hyper=hyper, miss=miss, 
                         rotation=rotation, start=start, 
                         control=list(nstart=control$nstart, opt=list(),
                                      lower=control$lower))
      ell <- c(ell, param$ell)
      
      # check convergence
      conv <- abs((param$ell - param.old$ell)/
                    ifelse(is.finite(param.old$ell), param.old$ell, 1)) < 
        control$epsilon
      check <- conv | iter >= control$maxit
      
      # reassign values
      param.old <- param
      
    }
  }
  return(c(param.old$fit, 
           list(param=list(B=param.old$B, psi=param.old$psi, 
                           Upsilon=param.old$Upsilon),
                gamma=hyper$gamma, Gamma=hyper$Gamma, ell=ell, iter=iter, 
                conv=conv)))
}



