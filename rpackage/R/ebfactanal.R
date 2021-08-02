# objective function for constrained penalty multiplier estimation
.fopt <- function(par, gamma, sizes, d, coefs) {
  
  s <- par[1]
  loggammagroups <- par[-1]
  
  return(sum(((s + 0.5*d)*sizes - 0.5*coefs*exp(-loggammagroups)/gamma)^2) + 
           sum(loggammagroups*sizes)^2)
}

# constrained EB update in logistic model
.constrebupdate <- function(gammagroups.old, coefs, gamma, sizes, d) {
  
  # optimse constrained problem
  opt <- optim(par=c(0, log(gammagroups.old)), fn=.fopt, gamma=gamma, 
               sizes=sizes, d=d, coefs=coefs, method="BFGS")
  
  # assigning new penalty multipliers
  return(exp(opt$par[-1]))
  
}

# param.old=param; hyper.old=hyper.old; groups=groups;
# m=colSums(is.na(x)); constrained=constrained
# empirical Bayes update
.ebupdate <- function(param.old, hyper.old, groups, m, constrained) {
  
  # auxiliary variables
  sizes <- rle(sort(groups[groups!=0]))$length
  G <- length(sizes)
  n <- nrow(param.old$Phi)
  d <- nrow(param.old$Xi)
  
  # parameter updates
  coefs <- sapply(1:G, function(g) {
    id <- which(groups==g)
    sum((sapply(param.old$Omega[id], function(om) {
      sum(diag(om))}) + colSums(param.old$Mu[, id]^2))*
        (n/2 + hyper.old$kappa[id])/param.old$zeta[id])})
  hyper <- hyper.old
  if(constrained) {
    hyper[["gammagroups"]] <- 
      .constrebupdate(hyper$gammagroups, coefs=coefs, 
                      gamma=hyper$gammatotal, sizes=sizes, d=d)
    gamma <- hyper[["gammagroups"]]*hyper[["gammatotal"]]
  } else {
    gamma <- coefs/(d*sizes)
  }
  
  # update and return hyper parameters
  hyper[["gamma"]][groups!=0] <- gamma[groups[groups!=0]]
  hyper[["const"]] <- .const(hyper, n=n, d=d, m=m)
  hyper[["elbo"]] <- .elbo(param.old, hyper, hyper[["const"]], m)
  return(hyper)
  
}

# x=x; factors=factors; hyper=hyper; start=start;
# groups=groups; rescale=rescale; constrained=constrained;
# control=control[c("maxit", "epsilon", "trace")]
# empirical Bayes factor analysis
ebfactanal <- function(x, factors, hyper, start, groups, rescale, constrained,
                       control) {
  
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
  hyper.old <- c(hyper, elbo=-Inf, const=0)
  
  # use one group of features if not given
  if(is.null(groups)) {
    groups <- rep(1, ncol(x))
  }
  
  # check whether constrained EB required
  if(is.null(constrained)) {
    constrained <- FALSE
  }
  
  # iterations
  iter <- 0
  check <- FALSE
  elbo <- numeric(0)
  diff <- setNames(vector(mode="list", length=length(param.old) - 1),
                   names(param.old)[names(param.old)!="elbo"])
  while(!check) {
    iter <- iter + 1
    if(control$trace) {
      cat("\r", "iteration ", iter, "      ", sep="")
    }
    
    # update the VB parameters
    param <- vbfactanal(x=x, factors=factors, hyper=hyper, 
                        start=param.old, rescale=FALSE,
                        control=list(maxit=1, epsilon=control$epsilon,
                                     trace=FALSE))$param
    elbo <- c(elbo, param$elbo)
    
    # track parameter changes
    diff <- sapply(names(diff), function(name) {
      new <- unlist(param[[name]])
      old <- unlist(param.old[[name]])
      c(diff[[name]], max(abs(new - old)/ifelse(old!=0, abs(old), 1)))},
      simplify=FALSE)
    
    # update the EB parameters
    hyper <- .ebupdate(param.old=param, hyper.old=hyper.old, groups=groups, 
                       m=colSums(is.na(x)), constrained=constrained)
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
  return(list(param=param, hyper=hyper, elbo=elbo, iter=iter, conv=conv, 
              const=hyper[["const"]], diff=diff))
  
}