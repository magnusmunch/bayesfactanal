# simulate from p(eta | beta)
.reta <- function(n, a, b0) {
  x <- numeric(n)
  accept <- 0
  samples <- 0
  while(accept < n) {
    samples <- samples + 1
    y <- rpg(1)
    u <- runif(1, 0, 1)
    prob <- exp(-y*a^2/(8 + 8*y*a))/sqrt(1 + y*a)
    if(b0!=0) {
      prob <- prob*exp(-b0^2*y/(2 + 2*y*a))*cosh(b0/(2 + 2*y*a))/
        cosh(b0/2)
    }
    if(u < prob) {
      accept <- accept + 1
      x[accept] <- y
    } 
  }
  return(list(x=x, rate=accept/samples))
}

# compute posterior expectation of regression estimates
coef.logvbfactreg <- function(fit, type="mc", nsamples=1000) {
  p <- ncol(fit$param$Mu$x)
  d <- nrow(fit$param$Mu$x)
  n <- nrow(fit$param$Phi)
  
  if(type=="mc") {
    coef <- matrix(nrow=p + 1, ncol=nsamples)
    tcholx <- lapply(fit$param$Omega$x, chol)
    tcholy <- chol(fit$param$Omega$y[[1]])
    for(k in 1:nsamples) {
      psiinv <- rgamma(p, shape=fit$hyper$kappa + n/2, 
                       rate=fit$param$zeta)
      B <- matrix(rnorm(d*p), nrow=d, ncol=p)
      B <- sapply(1:p, function(j) {colSums(tcholx[[j]]*B[, j])})
      B <- B + fit$param$Mu$x
      beta <- colSums(tcholy*rnorm(d + 1)) + as.numeric(fit$param$Mu$y)
      eta <- .reta(n=1, a=sum(beta[-1]^2), b0=beta[1])$x
      mat <- B %*% (t(B)*psiinv) + eta*matrix(beta[-1]) %*% t(beta[-1])
      diag(mat) <- diag(mat) + 1
      mat <- solve(mat)
      coef[1, k] <- beta[1]*(1 - sum(beta[-1]*colSums(mat*beta[-1]))*eta)
      coef[-1, k] <- colSums(B*colSums(mat*beta[-1]))*psiinv
    }
    vcoef <- cov(t(coef))
    coef <- rowMeans(coef)
  } else if(type=="taylor") {
    
  }
  return(list(coef=coef, vcoef=vcoef))
}

# predict new data in logistic model
predict.logvbfactreg <- function(newx, fit=NULL, coef=NULL, type="mc", 
                                 taylor=TRUE, nsamples=1000) {
  if(taylor) {
    if(is.null(coef)) {
      coef <- coef.logvbfactreg(fit=fit, type=type, nsamples=nsamples)
    }
    newx <- cbind(1, newx)
    c <- newx %*% coef$coef
    s <- rowSums((newx %*% coef$vcoef)*newx)
    pred <- .expit(c) + s*.expit(c)*(1 - .expit(c))^2/2 - 
      s*.expit(c)^2*(1 - .expit(c))/2
    vpred <- s*.expit(c)^2*(1 - .expit(c))^2
  }
  
  return(list(pred=pred, vpred=vpred))
}

# predict data points in logistic VB model
predict.logvbfactanal <- function(fit, newx=NULL, newy=NULL, simple=TRUE) {
  if(is.null(newx) & is.null(newy)) {
    predx <- fit$param$Phi %*% fit$param$Mu$x
    vpredx <- matrix(rep(fit$param$chi, nrow(fit$param$Phi)), 
                     ncol=length(fit$param$chi), nrow=nrow(fit$param$Phi),
                     byrow=TRUE)
    predy <- 1/(1 + exp(-cbind(1, fit$param$Phi) %*% fit$param$Mu$y))
    vpredy <- predy*(1 - predy)
  } else if(simple) {
    # newXi <- solve(Reduce("+", lapply(1:ncol(newx), function(j) {
    #   (nrow(fit$param$Upsilon)/2 + fit$hyper$kappa[j])*
    #     (fit$param$Omega[[j]] + as.vector(fit$param$Mu[, j]) %*% 
    #        t(as.vector(fit$param$Mu[, j])))/fit$param$zeta[j]})) + 
    #     diag(nrow(fit$param$Mu)))
    # newPhi <- newx %*% 
    #   (t(fit$param$Mu)*(nrow(fit$param$Upsilon)/2 + fit$hyper$kappa)/
    #      fit$param$zeta) %*% newXi
    # pred <- newPhi %*% fit$param$Mu
    # vpred <- matrix(rep(fit$param$chi, nrow(newx)), 
    #                 ncol=length(fit$param$chi), nrow=nrow(newx),
    #                 byrow=TRUE)
  } else {
    fit <- logvbfitupdate(fit, newx)
  }
  return(list(x=list(pred=predx, vpred=vpredx), 
              y=list(pred=predy, vpred=vpredy)))
}

# compute posterior expectation of regression estimates
coef.vbfactreg <- function(fit, outind=NULL, type="mc", nsamples=1000) {
  p <- ncol(fit$param$Mu) - 1
  d <- nrow(fit$param$Mu)
  n <- nrow(fit$param$Phi)
  if(is.null(outind)) {
    outind <- p + 1
  }
  if(type=="mc") {
    coef <- matrix(nrow=p, ncol=nsamples)
    tchols <- lapply(fit$param$Omega, chol)
    for(k in 1:nsamples) {
      psiinv <- rgamma(p, shape=fit$hyper$kappa[-outind] + n/2, 
                       rate=fit$param$zeta[-outind])
      B <- matrix(rnorm(d*(p + 1)), nrow=d, ncol=p + 1)
      B <- sapply(1:(p + 1), function(j) {colSums(tchols[[j]]*B[, j])})
      B <- B + fit$param$Mu
      beta <- B[, outind]
      B <- B[, -outind]
      mat <- B %*% (t(B)*psiinv)
      diag(mat) <- diag(mat) + 1
      coef[, k] <- colSums(B*colSums(solve(mat)*beta))*psiinv
    }
    vcoef <- cov(t(coef))
    coef <- rowMeans(coef)
  } else if(type=="taylor") {
    
  }
  return(list(coef=coef, vcoef=vcoef))
}

# compute posterior expectation of regression estimates
predict.vbfactreg <- function(newx, fit=NULL, coef=NULL, outind=NULL,
                              type="mc", nsamples=1000) {
  if(is.null(coef)) {
    coef <- coef.vbfactreg(fit=fit, outind=outind, type=type, nsamples=nsamples)
  }
  pred <- newx %*% coef$coef
  vpred <- NULL
  if(!is.null(coef$vcoef)) {
    vpred <- rowSums(newx*(newx %*% coef$vcoef))
  }
  return(list(pred=pred, vpred=vpred))
}

# compute regression EM estimates
coef.emfactreg <- function(fit, outind=NULL) {
  
  d <- nrow(fit$param$B)
  
  if(is.null(outind)) {
    outind <- ncol(fit$param$B)
  }
  B <- fit$param$B
  psi <- fit$param$psi[-outind]
  beta <- B[, outind]
  B <- B[, -outind]
  mat <- B %*% (t(B)/psi)
  diag(mat) <- diag(mat) + 1
  
  coef <- colSums(B*colSums(solve(mat)*beta))/psi
  vcoef <- NULL
  
  return(list(coef=coef, vcoef=vcoef))
}

# predict from EM regression
predict.emfactreg <- function(newx, fit=NULL, coef=NULL, outind=NULL) {
  if(is.null(coef)) {
    coef <- coef.emfactreg(fit=fit, outind=outind)
  }
  pred <- newx %*% coef$coef
  vpred <- NULL
  if(!is.null(coef$vcoef)) {
    vpred <- rowSums(newx*(newx %*% coef$vcoef))
  }
  return(list(pred=pred, vpred=vpred))
}

# predict data points in VB model
predict.vbfactanal <- function(fit, newx=NULL, simple=TRUE) {
  if(is.null(newx)) {
    pred <- fit$param$Phi %*% fit$param$Mu
    vpred <- matrix(rep(fit$param$chi, nrow(fit$param$Phi)), 
                    ncol=length(fit$param$chi), nrow=nrow(fit$param$Phi),
                    byrow=TRUE)
  } else if(simple) {
    newXi <- solve(Reduce("+", lapply(1:ncol(newx), function(j) {
      (nrow(fit$param$Upsilon)/2 + fit$hyper$kappa[j])*
        (fit$param$Omega[[j]] + as.vector(fit$param$Mu[, j]) %*% 
           t(as.vector(fit$param$Mu[, j])))/fit$param$zeta[j]})) + 
        diag(nrow(fit$param$Mu)))
    newPhi <- newx %*% 
      (t(fit$param$Mu)*(nrow(fit$param$Upsilon)/2 + fit$hyper$kappa)/
         fit$param$zeta) %*% newXi
    pred <- newPhi %*% fit$param$Mu
    vpred <- matrix(rep(fit$param$chi, nrow(newx)), 
                    ncol=length(fit$param$chi), nrow=nrow(newx),
                    byrow=TRUE)
  } else {
    fit <- vbfitupdate(fit, newx)
  }
  return(list(pred=pred, vpred=vpred))
}

# predict data points in EM model
predict.emfactanal <- function(fit, newx=NULL, type=1) {
  
  if(type==1) {
    sampleid <- which(apply(is.na(newx), 1, any))
    dups <- duplicated(is.na(newx[sampleid, , drop=FALSE]))
    miss <- list(groupid=lapply(1:length(which(!dups)), function(s) {
      which(apply(is.na(newx), 1, identical, 
                  is.na(newx[sampleid, , drop=FALSE])[
                    which(!dups)[s], ]))}),
      covid=apply(is.na(newx[sampleid, , drop=FALSE])[
        !dups, , drop=FALSE], 1, which))
    expvar <- .expfun(newx, fit$param, miss)
    pred <- expvar$expectation
    vpred <- matrix(rep(diag(expvar$variance), nrow(newx)), 
                    ncol=ncol(expvar$variance), nrow=nrow(newx),
                    byrow=TRUE)
  } else if(type==2) {
    mat1 <- t(fit$param$B)/fit$param$psi
    variance <- expectation <- fit$param$B %*% mat1
    diag(variance) <- diag(variance) + 1
    variance <- t(fit$param$B) %*% solve(variance) %*% fit$param$B
    diag(variance) <- diag(variance) + fit$param$psi
    expectation <- 2*expectation
    diag(expectation) <- diag(expectation) + 1
    expectation <- mat1 %*% solve(expectation) %*% t(mat1) %*% variance
    
    if(is.null(newx)) {
      newx <- fit$param$Upsilon
    }
    pred <- newx %*% expectation
    vpred <- matrix(rep(diag(variance), nrow(newx)), ncol=ncol(variance),
                    nrow=nrow(newx), byrow=TRUE)
  }
  
  return(list(pred=pred, vpred=vpred))
}

# "regularized" correlation matrix
regcor <- function(x, Gamma=NULL, nfolds=5, lower=0, upper=1, seed=NULL,
                   return.cor=TRUE) {
  
  # set the seed
  set.seed(seed)
  
  # foldids
  n <- nrow(x)
  foldid <- sample(rep(1:nfolds, times=round(c(rep(
    n %/% nfolds + as.numeric((n %% nfolds)!=0), times=n %% nfolds),
    rep(n %/% nfolds, times=nfolds - n %% nfolds)))))
  nk <- rle(sort(foldid))$lengths
  p <- ncol(x)
  
  # set penalty matrix
  if(is.null(Gamma)) {
    Gamma <- rep(1, ncol(x))
  }
  
  # calculate some needed values
  dval <- dvec <- ldetGamma <- NULL
  if(all(Gamma!=0)) {
    ldetGamma <- sum(log(Gamma))
    dval <- vector("list", length=nfolds)
    dvec <- vector("list", length=nfolds)
    for(k in 1:nfolds) {
      svd.x <- svd(t(t(scale(x[foldid!=k, , drop=FALSE]))/sqrt(Gamma)), 
                   nv=max(n - nk[k], p))  
      dvec[[k]] <- colSums((t(t(scale(x[foldid==k, , drop=FALSE]))/
                                sqrt(Gamma)) %*% svd.x$v)^2)
      dval[[k]] <- svd.x$d
    }
  }
  
  # Determine optimal penalty value and return
  opt <- optim(par=0.5, fn=.cvllcor, dval=dval, dvec=dvec, n=n, nk=nk, p=p, 
               ldetGamma=ldetGamma, x=x, Gamma=Gamma, foldid=foldid, 
               method="Brent", lower=lower, upper=upper)
  gamma <- opt$par
  return(list(gamma=gamma, cor=(1 - gamma)*cor(x) + gamma*diag(Gamma), opt=opt))
  
}

# calculate factor regression coefficients
factreg <- function(fit, type, burnin=0) {
  if(type=="vb") {
    res <- .vbbeta(param=fit$param, hyper=fit$hyper)
  } else if(type=="mcmc") {
    res <- .mcmcbeta(fit$param, burnin=burnin)
  } else if(type=="em") {
    res <- .embeta(param=fit$param)
  }
  return(res)
}
