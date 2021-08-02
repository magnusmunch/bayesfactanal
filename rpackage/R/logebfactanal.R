# param.old=param; hyper.old=hyper.old; groups=groups;
# mx=colSums(is.na(x)); y=y; constrained=constrained
# empirical Bayes update
.logebupdate <- function(param.old, hyper.old, groups, mx, y, constrained) {
  
  # auxiliary variables
  n <- nrow(param.old$Phi)
  d <- ncol(param.old$Phi)
  
  # parameter updates
  hyper <- hyper.old
  
  # parameter updates for linear variables
  if(any(groups$x!=0)) {
    sizes <- rle(sort(groups$x[groups$x!=0]))$length
    G <- length(sizes)
    coefs <- sapply(1:G, function(g) {
      id <- which(groups$x==g)
      sum((sapply(param.old$Omega$x[id], function(om) {
        sum(diag(om))}) + colSums(param.old$Mu$x[, id]^2))*
          (n/2 + hyper.old$kappa[id])/param.old$zeta[id])})
    if(constrained$x) {
      hyper[["gammagroups"]][["x"]] <- 
        .constrebupdate(hyper$gammagroups$x, coefs=coefs, 
                        gamma=hyper$gammatotal$x, sizes=sizes, d=d)
      gamma <- hyper[["gammagroups"]][["x"]]*hyper[["gammatotal"]][["x"]]
    } else {
      gamma <- coefs/(d*sizes)
    }
    hyper[["gamma"]][["x"]][groups$x!=0] <- gamma[groups$x[groups$x!=0]]
  }
  
  # parameter updates for logistic variables
  if(any(groups$y!=0)) {
    sizes <- rle(sort(groups$y[groups$y!=0]))$length
    G <- length(sizes)
    coefs <- sapply(1:G, function(g) {
      id <- which(groups$y==g)
      sum(sapply(param.old$Omega$y[id], function(om) {
        sum(diag(om))}) + colSums(param.old$Mu$x[, id]^2))})
    if(constrained$y) {
      hyper[["gammagroups"]][["y"]] <- 
        .constrebupdate(hyper$gammagroups$y, coefs=coefs, 
                        gamma=hyper$y$gammatotal, sizes=sizes, d=d)
      gamma <- hyper[["gammagroups"]][["y"]]*hyper[["gammatotal"]][["y"]]
    } else {
      gamma <- coefs/(d*sizes)
    }
    hyper[["gamma"]][["y"]][groups$y!=0] <- gammay[groups$y[groups$y!=0]]  
  }
  
  # update elbo and return hyper parameters
  if(any(groups$x!=0) | any(groups$y!=0)) {
    linconst <- .const(hyper=list(kappa=hyper$kappa, nu=hyper$nu, 
                                  gamma=hyper$gamma$x), n=n, d=d, m=mx)
    const <- .logconst(hyper=hyper, d=d, linconst=linconst)
    linparam <- list(Phi=param.old$Phi, Xi=param.old$Xi, Mu=param.old$Mu$x, 
                     Omega=param.old$Omega$x, zeta=param.old$zeta, 
                     Upsilon=param.old$Upsilon$x, chi=param.old$chi)
    linhyper <- list(kappa=hyper$kappa, nu=hyper$nu, gamma=hyper$gamma$x)
    linelbo <- .elbo(linparam, linhyper, 0, mx)
    hyper[["const"]] <- const
    hyper[["elbo"]] <- .logelbo(param.old, hyper, const, linelbo, y)
  }
  return(hyper)
  
}

# x=x; y=y; factors=factors; hyper=hyper; start=start;
# groups=groups; rescale=rescale;
# constrained=constrained;
# control=control[c("maxit", "epsilon", "trace")]
# empirical Bayes factor analysis
logebfactanal <- function(x, y, factors, hyper, start=NULL, 
                          groups=list(x=NULL, y=NULL), rescale, 
                          constrained=list(x=FALSE, y=FALSE), control) {
  
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
  hyper.old <- c(hyper, elbo=-Inf, const=0)
  
  # use one group of features for x and none for y if not given
  if(is.null(groups$x)) {
    groups$x <- rep(1, ncol(x))
  }
  if(is.null(groups$y)) {
    groups$y <- rep(0, ncol(y))
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
  #                        names(param.old))
  while(!check) {
    iter <- iter + 1
    if(control$trace) {
      cat("\r", "iteration ", iter, "      ", sep="")
    }
    
    # update the parameters
    param <- logvbfactanal(x=x, y=y, factors=factors, hyper=hyper, 
                           start=param.old, rescale=FALSE,
                           control=list(maxit=1, epsilon=control$epsilon, 
                                        trace=FALSE))$param
    elbo <- c(elbo, param$elbo)
    
    # track parameter changes
    diff <- sapply(names(diff), function(name) {
      new <- unlist(param[[name]])
      old <- unlist(param.old[[name]])
      if(name=="Upsilon") {
        new <- new[is.na(unlist(list(x, y)))]
        old <- old[is.na(unlist(list(x, y)))]
      } 
      c(diff[[name]], max(abs(new - old)/ifelse(old!=0, abs(old), 1)))},
      simplify=FALSE)
    
    # update the EB parameters
    hyper <- .logebupdate(param.old=param, hyper.old=hyper.old, groups=groups, 
                          mx=colSums(is.na(x)), y=y, constrained=constrained)
    elbo <- c(elbo, hyper$elbo)
    
    # check convergence and iteration number
    conv.diff <- all(sapply(diff, tail, n=1) < control$epsilon)
    conv.elbo <- abs((hyper$elbo - hyper.old$elbo)/
                       ifelse(is.finite(hyper.old$elbo), hyper.old$elbo, 1)) <
      control$epsilon
    conv <- conv.diff | conv.elbo
    check <- conv | iter >= control$maxit
    
    param.old <- param
    hyper.old <- hyper
  }
  
  # rescale if correlation matrix is modelled and return
  if(rescale) {
    fact <- colSums(param$Mu$x^2) + 
      sapply(param$Omega$x, function(om) {
        sum(diag(om))}) + param$zeta/(nrow(x)/2
                                      # + d/2
                                      + hyper$kappa - 1)
    param$Mu$x <- t(t(param$Mu$x)/sqrt(fact))
    param$Omega$x <- sapply(1:length(param$Omega$x), function(jx) {
      param$Omega$x[[jx]]/fact[jx]}, simplify=FALSE)
    param$zeta <- param$zeta/fact
  }
  return(list(param=param, hyper=hyper, elbo=elbo, iter=iter, conv=conv, 
              const=hyper[["const"]], diff=diff))
}

