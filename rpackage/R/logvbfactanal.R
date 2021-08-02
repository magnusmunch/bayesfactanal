# variational Bayes update
.logvbupdate <- function(y, param.old, hyper, const) {
  
  # auxiliary variables
  px <- length(param.old$Omega$x)
  py <- length(param.old$Omega$y)
  n <- nrow(param.old$Phi)
  mx <- colSums(is.na(x))
  my <- colSums(is.na(y))
  d <- ncol(param.old$Phi)
  
  # parameter updates
  delta <- sqrt(
    sapply(1:py, function(jy) {
    mu <- param.old$Mu$y[, jy]
    omega <- param.old$Omega$y[[jy]]
    sapply(param.old$Xi, function(xi) {
      sum(colSums(mu[-1]*xi)*mu[-1]) +
      sum(xi*omega[-1, -1])}) + 
    colSums(mu[-1]*t(param.old$Phi))^2 + 
      rowSums((param.old$Phi %*% omega[-1, -1])*param.old$Phi) +
    2*colSums(t(param.old$Phi)*omega[-1, 1]) +
    2*colSums(t(param.old$Phi)*mu[-1])*mu[1] + mu[1]^2 + omega[1, 1]}))
  aux1 <- tanh(delta/2)/(2*delta)
  aux2 <- (n/2 
           # + d/2
           + hyper$kappa)/param.old$zeta
  Xi <- lapply(1:n, function(i) {
    mat1 <- Reduce("+", lapply(1:px, function(jx) {
      aux2[jx]*(param.old$Omega$x[[jx]] + as.vector(param.old$Mu$x[, jx]) %*% 
                  t(as.vector(param.old$Mu$x[, jx])))}))
    mat2 <- Reduce("+", lapply(1:py, function(jy) {
      aux1[i, jy]*(param.old$Omega$y[[jy]][-1, -1] + 
                     as.vector(param.old$Mu$y[-1, jy]) %*% 
                     t(as.vector(param.old$Mu$y[-1, jy])))}))
    solve(mat1 + mat2 + diag(d))})
  Phi <- t(sapply(1:n, function(i) {
    mat <- colSums(t(param.old$Mu$x)*aux2*param.old$Upsilon$x[i, ]) +
      colSums(t(param.old$Mu$y[-1, ])*
                (param.old$Upsilon$y[i, ] - 0.5 - 
                   aux1[i, ]*param.old$Mu$y[1, ])) -
      colSums(t(sapply(param.old$Omega$y, function(om) {om[1, -1]}))*aux1[i, ])
    mat %*% Xi[[i]]}))
  aux3 <- Reduce("+", Xi) + t(Phi) %*% Phi
  Omega <- list()
  Omega$x <- lapply(1:px, function(jx) {
    solve(aux3 + diag(d)/hyper$gamma$x[jx])/aux2[jx]})
  Mu <- list()
  Mu$x <- sapply(1:px, function(jx) {
    (aux2[jx]*Omega$x[[jx]]) %*% t(Phi) %*% param.old$Upsilon$x[, jx]})
  Omega$y <- lapply(1:py, function(jy) {
    om <- matrix(nrow=d + 1, ncol=d + 1)
    om[1, 1] <- sum(aux1[, jy])
    om[1, -1] <- om[-1, 1] <- colSums(aux1[, jy]*Phi)
    om[2:(d + 1), 2:(d + 1)] <- Reduce("+", lapply(1:n, function(i) {
    aux1[i, jy]*(Xi[[i]] + t(Phi[i, , drop=FALSE]) %*% 
                   Phi[i, , drop=FALSE])})) + diag(d)/hyper$gamma$y[jy]
    solve(om)})
  Mu$y <- as.matrix(sapply(1:py, function(jy) {
    colSums(c(sum(param.old$Upsilon$y[, jy] - 0.5), 
              colSums(Phi*(param.old$Upsilon$y[, jy] - 0.5)))*Omega$y[[jy]])}))
  zeta <- sapply(1:px, function(jx) {
    sum(param.old$Upsilon$x[, jx]^2)/2 + mx[jx]*param.old$chi[jx] - 
      sum(colSums(Mu$x[, jx]*t(Phi))*param.old$Upsilon$x[, jx]) +
      sum(aux3*Omega$x[[jx]])/2 + 
      sum(colSums(Mu$x[, jx]*aux3)*Mu$x[, jx])/2 + 
      sum(diag(Omega$x[[jx]]))/(2*hyper$gamma$x[jx]) + 
      sum(Mu$x[, jx]^2)/(2*hyper$gamma$x[jx]) + hyper$nu[jx]})
  Upsilon <- list()
  Upsilon$x <- x
  if(any(is.na(Upsilon$x))) {
    Upsilon$x[is.na(Upsilon$x)] <- apply(
      which(is.na(Upsilon$x), arr.ind=TRUE), 1, function(ind) {
        sum(Phi[ind[1], ]*Mu$x[, ind[2]])})
  }
  chi <- zeta/(n/2
               # + d/2
               + hyper$kappa)
  Upsilon$y <- y
  if(any(is.na(Upsilon$y))) {
    Upsilon$y[is.na(Upsilon$y)] <- apply(
      which(is.na(Upsilon$y), arr.ind=TRUE), 1, function(ind) {
        1/(1 + exp(-sum(t(Phi[ind[1], ])*Mu$y[-1, ind[2]]) - 
                     Mu$y[1, ind[2]]))})
  }
  
  # return values
  param <- list(delta=delta, Phi=Phi, Xi=Xi, Mu=Mu, Omega=Omega, zeta=zeta, 
                Upsilon=Upsilon, chi=chi)
  
  # add evidence lower bound and return
  linparam <- list(Phi=Phi, Xi=Xi, Mu=Mu$x, Omega=Omega$x, zeta=zeta, 
                   Upsilon=Upsilon$x, chi=chi)
  linhyper <- list(kappa=hyper$kappa, nu=hyper$nu, gamma=hyper$gamma$x)
  linelbo <- .elbo(linparam, linhyper, 0, mx)
  param <- c(param, elbo=.logelbo(param, hyper, const, linelbo, y))
  return(param)
  
}

# factors=est.factors; 
# start=NULL; rescale=TRUE; 
# control=list(maxit=100, epsilon=1e-4)
# variational Bayes factor analysis
logvbfactanal <- function(x, y, factors, hyper, start=NULL, rescale, control) {
  
  # convert data to matrix
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  # starting values
  param.old <- start
  if(is.null(start)) {
    param.old <- list(delta=matrix(0, nrow=nrow(y), ncol=ncol(y)),
                      Phi=matrix(rnorm(factors*nrow(x)), nrow=nrow(x)), 
                      Xi=rep(list(diag(factors)), nrow(x)), 
                      Mu=list(x=matrix(rnorm(factors*ncol(x)), ncol=ncol(x)),
                              y=rbind(log(colMeans(y, na.rm=TRUE)) - 
                                        log(1 - colMeans(y, na.rm=TRUE)), 
                                      matrix(rnorm(factors*ncol(y)), 
                                             nrow=factors, ncol=ncol(y)))), 
                      Omega=c(x=list(rep(list(diag(factors)), ncol(x))), 
                              y=list(rep(list(diag(factors + 1)), ncol(y)))), 
                      zeta=hyper$nu/hyper$kappa, 
                      Upsilon=list(x=replace(x, which(is.na(x)), 0),
                                   y=replace(y, which(is.na(y)), 0.5)), 
                      chi=rep(1, ncol(x)),
                      elbo=-Inf)
  } 
  
  # calculate constant part of ELBO
  linconst <- .const(hyper=list(kappa=hyper$kappa, nu=hyper$nu, 
                                gamma=hyper$gamma$x), 
                     n=nrow(x), d=factors, m=colSums(is.na(x)))
  const <- .logconst(hyper=hyper, d=factors, linconst=linconst)
  
  # iterations
  iter <- 0
  check <- FALSE
  elbo <- numeric(0)
  diff <- setNames(vector(mode="list", length=length(param.old) - 1),
                   names(param.old)[names(param.old)!="elbo"])
  # maxabsdiff <- setNames(vector(mode="list", length=length(param.old)),
  #                       names(param.old))
  while(!check) {
    iter <- iter + 1
    if(control$trace) {
      cat("\r", "iteration ", iter, "      ", sep="")
    }
    
    # update the parameters
    param <- .logvbupdate(y=y, param.old=param.old, hyper=hyper, const=const)
    elbo <- c(elbo, param$elbo)
    diff <- sapply(names(diff), function(name) {
      new <- unlist(param[[name]])
      old <- unlist(param.old[[name]])
      if(name=="Upsilon") {
        new <- new[is.na(unlist(list(x, y)))]
        old <- old[is.na(unlist(list(x, y)))]
      } 
      c(diff[[name]], max(abs(new - old)/ifelse(old!=0, abs(old), 1)))},
      simplify=FALSE)
    # maxabsdiff <- sapply(names(maxabsdiff), function(name) {
    #   if(name=="Upsilon") {
    #     midx <- is.na(x)
    #     midy <- is.na(y)
    #     diff <- abs(c(param$Upsilon$x[midx] - param.old$Upsilon$x[midx],
    #                   param$Upsilon$y[midy] - param.old$Upsilon$y[midy]))/
    #       abs(ifelse(c(param.old$Upsilon$x[midx], param.old$Upsilon$y[midy])!=0, 
    #                  c(param.old$Upsilon$x[midx], param.old$Upsilon$y[midy]), 
    #                  1))
    #     
    #   } else {
    #     diff <- abs(unlist(param[[name]]) - unlist(param.old[[name]]))/
    #       abs(ifelse(unlist(param.old[[name]])!=0, unlist(param.old[[name]]), 
    #                  1))
    #   }
    #   c(maxabsdiff[[name]], max(diff))}, simplify=FALSE)
    
    
    # check convergence and iteration number
    conv1 <- all(sapply(diff, tail, n=1) < control$epsilon)
    conv2 <- abs((param$elbo - param.old$elbo)/
                  ifelse(is.finite(param.old$elbo), param.old$elbo, 1)) <
      control$epsilon
    conv <- conv1 | conv2
    # conv <- abs((param$elbo - param.old$elbo)/
    #               ifelse(is.finite(param.old$elbo), param.old$elbo, 1)) < 
    #   control$epsilon
    check <- conv | iter >= control$maxit
    
    param.old <- param
  }
  
  # rescale if correlation matrix is modelled and return
  if(rescale) {
    fact <- colSums(param$Mu$x^2) + 
      sapply(param$Omega$x, function(om) {
        sum(diag(om))}) + param$zeta/(nrow(x)/2
                                      # + factors/2
                                      + hyper$kappa - 1)
    param$Mu$x <- t(t(param$Mu$x)/sqrt(fact))
    param$Omega$x <- sapply(1:length(param$Omega$x), function(jx) {
      param$Omega$x[[jx]]/fact[jx]}, simplify=FALSE)
    param$zeta <- param$zeta/fact
  }
  return(list(param=param, hyper=hyper, elbo=elbo, iter=iter, conv=conv, 
              const=const, 
              # maxabsdiff=maxabsdiff,
              diff=diff))
}

# Bayesian updating of logistic VB model
logvbfitupdate <- function(fit, x, y, start=NULL, simple=TRUE, rescale, 
                           control) {
  
}
