# MCMC Bayes update
.mcmcupdate <- function(x, sample, hyper) {
  
  # auxiliary variables
  p <- ncol(sample$B)
  n <- nrow(sample$Lambda)
  d <- nrow(sample$B)
  
  # sample new parameters
  Xi <- sample$B %*% (t(sample$B)/sample$psi)
  diag(Xi) <- diag(Xi) + 1
  Xi <- solve(Xi)
  Lambda <- t(sapply(1:n, function(i) {
    rmvnorm(1, colSums(Xi*colSums(t(sample$B)*(sample$x[i, ]/
                                                 sample$psi))), Xi)}))
  B <- sapply(1:p, function(j) {
    Omega <- sample$psi[j]*
      solve(t(Lambda) %*% Lambda + diag(d)/hyper$gamma[j])
    rmvnorm(1, colSums(Omega*colSums(Lambda*sample$x[, j]))/sample$psi[j], 
            Omega)})
  psi <- 1/rgamma(p, n/2 + hyper$kappa, 
                  colSums((sample$x - Lambda %*% B)^2)/2 + hyper$nu)
  if(any(is.na(x))) {
    x[is.na(x)] <- apply(
      which(is.na(x), arr.ind=TRUE), 1, function(ind) {
        rnorm(1, sum(Lambda[ind[1], ]*B[, ind[2]]), psi[ind[2]])})
  }
  
  # return values
  sample <- list(Lambda=Lambda, B=B, psi=psi, x=x)
  return(sample)
  
}

# MCMC Bayes factor analysis
mcmcfactanal <- function(x, factors, hyper, start, control) {
  
  # setup parallel computations
  cluster <- makeForkCluster(nnodes=control$ncores)
  if(control$ncores > 1) {
    registerDoParallel(cluster)
  } else {
    registerDoSEQ()
  }
  
  # combine chains function
  comb <- function(a, b) {
    out <- lapply(names(a), function(nm) {
      unname(abind(a[[nm]], b[[nm]], along=length(dim(a[[nm]])) + 1))})
    names(out) <- names(a)
    return(out)
  }
  
  # run the separate chains
  res <- foreach(k=c(1:control$nchains), .combine=comb) %dopar% {
    
    # starting values
    sample <- start[[k]]
    if(is.null(start)) {
      sample <- list(B=matrix(rnorm(ncol(x)*factors), nrow=factors, 
                              ncol=ncol(x)), 
                     Lambda=matrix(rnorm(nrow(x)*factors), nrow=nrow(x), 
                                   ncol=factors), 
                     psi=1/rgamma(ncol(x), hyper$kappa, hyper$nu),
                     x=replace(x, which(is.na(x)), 0))
    } 
    
    # object to store samples
    param <- list(B=array(NA, c(dim(sample$B), control$nsamples)),
                  Lambda=array(NA, c(dim(sample$Lambda), control$nsamples)),
                  psi=array(NA, dim=c(length(sample$psi), control$nsamples)),
                  x=array(NA, dim=c(dim(x), control$nsamples)))
    
    # iterations
    for(k in 1:control$nsamples) {
      sample <- .mcmcupdate(x, sample, hyper)
      param$Lambda[, , k] <- sample$Lambda
      param$B[, , k] <- sample$B
      param$psi[, k] <- sample$psi
      param$x[, , k] <- sample$x
    }
    
    # return value
    param
  }
  if(control$ncores > 1) {stopCluster(cluster)}
  
  # return object
  return(list(param=res, nsamples=control$nsamples))
}




