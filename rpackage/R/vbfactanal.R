# variational Bayes update
.vbupdate <- function(x, param.old, hyper, const) {
  
  # auxiliary variables
  p <- length(param.old$Omega)
  n <- nrow(param.old$Phi)
  d <- nrow(param.old$Xi)
  m <- colSums(is.na(x))
  
  # parameter updates
  Xi <- solve(Reduce("+", lapply(1:p, function(j) {
    (n/2 + d/2 + hyper$kappa[j])*
      (param.old$Omega[[j]] + as.vector(param.old$Mu[, j]) %*% 
         t(as.vector(param.old$Mu[, j])))/param.old$zeta[j]})) + diag(d))
  Phi <- param.old$Upsilon %*% 
    (t(param.old$Mu)*(n/2
                      # + d/2
                      + hyper$kappa)/param.old$zeta) %*% Xi
  aux <- n*Xi + t(Phi) %*% Phi
  Omega <- lapply(1:p, function(j) {
    param.old$zeta[j]*solve(aux + diag(d)/hyper$gamma[j])/
      (n/2
       # + d/2
       + hyper$kappa[j])})
  Mu <- sapply(1:p, function(j) {
    ((n/2
      # + d/2
      + hyper$kappa[j])*Omega[[j]]/param.old$zeta[j]) %*% 
      t(Phi) %*% param.old$Upsilon[, j]})
  zeta <- sapply(1:p, function(j) {
    sum(param.old$Upsilon[, j]^2)/2 + m[j]*param.old$chi[j]/2 - 
      sum(colSums(Mu[, j]*t(Phi))*param.old$Upsilon[, j]) +
      sum(diag(aux %*% Omega[[j]]))/2 + 
      sum(colSums(Mu[, j]*aux)*Mu[, j])/2 + 
      sum(diag(Omega[[j]]))/(2*hyper$gamma[j]) + 
      sum(Mu[, j]^2)/(2*hyper$gamma[j]) + hyper$nu[j]})
  Upsilon <- x
  if(any(is.na(Upsilon))) {
    Upsilon[is.na(Upsilon)] <- apply(
      which(is.na(Upsilon), arr.ind=TRUE), 1, function(ind) {
        sum(Phi[ind[1], ]*Mu[, ind[2]])})
  }
  chi <- zeta/(n/2 + d/2 + hyper$kappa)
  
  # return values
  param <- list(Phi=Phi, Xi=Xi, Mu=Mu, Omega=Omega, zeta=zeta, Upsilon=Upsilon,
                chi=chi)
  
  # add evidence lower bound and return
  param <- c(param, elbo=.elbo(param, hyper, const, m))
  return(param)
  
}

.vbstart <- function() {
  
}

# variational Bayes factor analysis
vbfactanal <- function(x, factors, hyper, start=NULL, rescale, control) {
  
  # save call
  call <- match.call()
  
  # starting values
  param.old <- start
  if(is.null(start)) {
    param.old <- list(Phi=matrix(rnorm(factors*nrow(x)), nrow=nrow(x)), 
                      Xi=diag(factors), 
                      Mu=matrix(rnorm(factors*ncol(x)), ncol=ncol(x)), 
                      Omega=rep(list(diag(factors)), ncol(x)), 
                      zeta=hyper$nu/hyper$kappa, 
                      Upsilon=replace(x, which(is.na(x)), 0), 
                      chi=rep(1, ncol(x)), elbo=-Inf)
  } 
  
  # calculate constant part of ELBO
  const <- .const(hyper, n=nrow(x), d=factors, colSums(is.na(x)))
  
  # iterations
  iter <- 0
  check <- FALSE
  elbo <- numeric(0)
  while(!check) {
    iter <- iter + 1
    if(control$trace) {
      cat("\r", "iteration ", iter, "      ", sep="")
    }
    
    # update the parameters
    param <- .vbupdate(x=x, param.old=param.old, hyper=hyper, const=const)  
    elbo <- c(elbo, param$elbo)
    
    # check convergence and iteration number
    conv <- abs((param$elbo - param.old$elbo)/
                  ifelse(is.finite(param.old$elbo), param.old$elbo, 1)) < 
      control$epsilon
    check <- conv | iter >= control$maxit
    
    param.old <- param
  }
  
  # rescale if correlation matrix is modelled and return
  if(rescale) {
    fact <- colSums(param$Mu^2) + 
      sapply(param$Omega, function(om) {sum(diag(om))}) + 
      param$zeta/(nrow(x)/2
                  # + factors/2
                  + hyper$kappa - 1)
    param$Mu <- t(t(param$Mu)/sqrt(fact))
    param$Omega <- sapply(1:length(param$Omega), function(j) {
      param$Omega[[j]]/fact[j]}, simplify=FALSE)
    param$zeta <- param$zeta/fact
  }
  return(list(call=call, param=param, hyper=hyper, elbo=elbo, iter=iter, 
              conv=conv, const=const))
}

vbfitupdate <- function(fit, x, y, start=NULL, simple=TRUE, rescale, control) {
  
}
