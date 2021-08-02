### libraries
library(bayesfactanal)
library(FMradio)
library(mvtnorm)
library(glmnet)
library(RColorBrewer)

## functions
# simulation functions
sim1 <- function(nreps, ntest, n, m, p, d, groups, communality, trace=TRUE,
                 seed=NULL) {

  # set seed
  set.seed(seed)

  # setup results objects
  methods <- c("true", "ridge", "lasso",
               paste0("FMradio", 1:length(m)),
               paste0("emfactanal1.", 1:length(m)),
               # paste0("emfactanal2.", 1:length(m)),
               paste0("vbfactanal1.", 1:length(m)),
               # paste0("vbfactanal2.", 1:length(m)),
               paste0("ebfactanal1.", 1:length(m)),
               # paste0("ebfactanal2.", 1:length(m)),
               "null")
  fact <- numeric(nreps)
  conv <- pmse <- emse <-
    matrix(NA, nrow=nreps, ncol=length(methods),
           dimnames=list(NULL, methods))
  coln <- apply(expand.grid(1:length(unique(groups)), 1:length(m),
                            c("ebfactanal1."
                              # , "ebfactanal2."
                              )), 1, function(s) {
                              paste0(s[3], s[2], ".", s[1])})
  gamma <- matrix(NA, nrow=nreps, ncol=1*length(unique(groups))*length(m),
                  dimnames=list(NULL, coln))

  # calculate loadings matrix
  psi <- 1
  sigma2 <- 1
  B <- model.matrix(~ -1 + factor(rep(1:d, each=p/d)))
  B <- t(B*sqrt(psi*communality[-length(communality)][groups]/
                  (rowSums(B)*(1 - communality[-length(communality)][groups]))))
  B <- cbind(B, sqrt(sigma2*communality[length(communality)]/
                       (d*(1 - communality[length(communality)]))))
  beta <- colSums(B[, p + 1]*(solve(B[, -(p + 1)] %*% t(B[, -(p + 1)]) +
                                      diag(d)/psi) %*% B[, -(p + 1)]))*
    sqrt((colSums(B[, -(p + 1)]^2) + psi)/(sum(B[, p + 1]^2) + sigma2))

  # simulations
  for(r in 1:nreps) {
    if(trace) {cat("rep ", r, "\r")}

    # simulate data
    xy <- scale(rmvnorm(n + tail(m, 1), rep(0, (p + 1)),
                        t(B) %*% B + diag(c(rep(psi, p), sigma2))))

    # fit ridge and lasso models
    fit.ridge <- cv.glmnet(scale(xy[1:n, -(p + 1)]),
                           scale(xy[1:n, p + 1])[, 1],
                           alpha=0, intercept=FALSE)
    fit.lasso <- cv.glmnet(scale(xy[1:n, -(p + 1)]),
                           scale(xy[1:n, p + 1])[, 1],
                           alpha=1, intercept=FALSE)

    # estimate models
    for(s in 1:length(m)) {

      # estimate penalty parameters
      fit.cormat <- bayesfactanal::regcor(xy[1:(n + m[s]), -(p + 1)])
      cormat <- fit.cormat$cor
      est.gamma <- fit.cormat$gamma

      # estimate number of factors
      est.factors <- sum(eigen(cor(xy[1:(n + m[s]), -(p + 1)]),
                               only.values=TRUE)$values > 1)
      est.factors <- max(est.factors, 2)

      # set hyperparameters
      hyper1 <- hyper2 <-
        # hyper3 <- hyper4 <-
        list(kappa=rep(9, p + 1), nu=rep(4, p + 1))
      hyper1 <- c(hyper1,
                  list(gamma=(hyper1$kappa - 1)/(hyper1$nu*2*est.factors)))
      # hyper2 <- c(hyper2,
      #             list(gamma=c((hyper2$kappa[-(p + 1)] - 1)/
      #                            (hyper2$nu[-(p + 1)]*2*est.factors), y=Inf)))
      # hyper3 <- c(hyper3,
      #             list(gamma=(hyper3$kappa - 1)/(hyper3$nu*2*est.factors)),
      #             list(gammagroups=c(1, 1)),
      #             list(gammatotal=c((mean(hyper3$kappa[-(p + 1)]) - 1)/
      #                                 (mean(hyper3$nu[-(p + 1)])*2*
      #                                    est.factors))))
      hyper2 <- c(hyper1,
                  list(gammagroups=c(1, 1)),
                  list(gammatotal=c((mean(hyper1$kappa[-(p + 1)]) - 1)/
                                      (mean(hyper1$nu[-(p + 1)])*2*
                                         est.factors))))
      # hyper4 <- c(hyper4,
      #             list(gamma=c((hyper4$kappa[-(p + 1)] - 1)/
      #                            (hyper3$nu[-(p + 1)]*2*est.factors), y=Inf)),
      #             list(gammagroups=c(1, 1)),
      #             list(gammatotal=c((mean(hyper4$kappa[-(p + 1)]) - 1)/
      #                                 (mean(hyper4$nu[-(p + 1)])*2*
      #                                    est.factors))))

      # fit FMradio models
      fit.FMradio <- tryCatch({
        fit <- mlFA(cormat, est.factors);
        fit <- c(fit, cv.glmnet(as.matrix(
          facScore(scale(xy[1:n, -(p + 1)]), LM=fit$Loadings,
                   UM=fit$Uniqueness)), scale(xy[1:n, p + 1])[, 1],
          alpha=0, intercept=FALSE));
        c(fit, list(param=list(
          B=cbind(t(fit$Loadings),
                  coef(structure(fit, class="cv.glmnet"),
                       s="lambda.min")[-1]),
          psi=c(diag(fit$Uniqueness), 0))))},
        error=function(e) {list()})
      assign(paste0("fit.FMradio", s), fit.FMradio)

      # fit EM factanal models
      fit.emfactanal1 <- tryCatch({
        emfactanal(xy[1:(n + m[s]), ], factors=est.factors,
                   gamma=est.gamma, Gamma=rep(1, p + 1))},
        error=function(e) {list()})
      assign(paste0("fit.emfactanal1.", s), fit.emfactanal1)
      # fit.emfactanal2 <- emfactanal(xy[1:(n + m[s]), ], factors=est.factors,
      #                               gamma=est.gamma, Gamma=c(rep(1, p), 0))
      # assign(paste0("fit.emfactanal2.", s), fit.emfactanal2)

      # fit Bayesian factanal models
      fit.vbfactanal1 <- bayesfactanal(x=xy[1:(n + m[s]), ],
                                       factors=est.factors, hyper=hyper1,
                                       method="vb", rescale=TRUE)
      assign(paste0("fit.vbfactanal1.", s), fit.vbfactanal1)
      # fit.vbfactanal2 <- bayesfactanal(x=xy[1:(n + m[s]), ],
      #                                  factors=est.factors, hyper=hyper2,
      #                                  method="vb", rescale=TRUE)
      # assign(paste0("fit.vbfactanal2.", s), fit.vbfactanal2)

      # fit empirical Bayesian factanal models
      fit.ebfactanal1 <- bayesfactanal(
        x=xy[1:(n + m[s]), ], factors=est.factors, groups=c(groups, 0),
        hyper=hyper2, method="eb", rescale=TRUE, constrained=TRUE)
      assign(paste0("fit.ebfactanal1.", s), fit.ebfactanal1)
      # fit.ebfactanal2 <- bayesfactanal(
      #   x=xy[1:(n + m[s]), ], factors=est.factors,  groups=c(groups, 0),
      #   hyper=hyper4, method="eb", rescale=TRUE, constrained=TRUE)
      # assign(paste0("fit.ebfactanal2.", s), fit.ebfactanal2)
    }

    # simulate test data
    xytest <- scale(rmvnorm(ntest, rep(0, (p + 1)),
                            t(B) %*% B + diag(c(rep(psi, p), sigma2))))
    xtest <- xytest[, -(p + 1)]
    ytest <- xytest[, p + 1]

    # calculate regression parameters
    best <- cbind(true=beta,
                  ridge=coef(fit.ridge, s="lambda.min")[-1],
                  lasso=coef(fit.lasso, s="lambda.min")[-1],
                  matrix(sapply(1:length(m), function(s) {
                    fit <- eval(parse(text=paste0("fit.FMradio", s)))
                    if(length(fit)==0) {
                      rep(NA, p)
                    } else {
                      coef.emfactreg(fit, outind=p + 1)$coef
                    }}), ncol=length(m),
                    dimnames=list(NULL, paste0("FMradio", 1:length(m)))),
                  matrix(sapply(1:length(m), function(s) {
                    fit <- eval(parse(text=paste0("fit.emfactanal1.", s)))
                    if(length(fit)==0) {
                      rep(NA, p)
                    } else {
                      coef.emfactreg(fit, outind=p + 1)$coef
                    }}), ncol=length(m),
                    dimnames=list(NULL, paste0("emfactanal1.", 1:length(m)))),
                  # matrix(sapply(1:length(m), function(s) {
                  #   fit <- eval(parse(text=paste0("fit.emfactanal2.", s)))
                  #   coef.emfactreg(fit, outind=p + 1)$coef}), ncol=length(m),
                  #   dimnames=list(NULL, paste0("emfactanal2.", 1:length(m)))),
                  matrix(sapply(1:length(m), function(s) {
                    fit <- eval(parse(text=paste0("fit.vbfactanal1.", s)))
                    coef.vbfactreg(fit, outind=p + 1)$coef}), ncol=length(m),
                    dimnames=list(NULL, paste0("vbfactanal1.", 1:length(m)))),
                  # matrix(sapply(1:length(m), function(s) {
                  #   fit <- eval(parse(text=paste0("fit.vbfactanal2.", s)))
                  #   coef.vbfactreg(fit, outind=p + 1)$coef}), ncol=length(m),
                  #   dimnames=list(NULL, paste0("vbfactanal2.", 1:length(m)))),
                  matrix(sapply(1:length(m), function(s) {
                    fit <- eval(parse(text=paste0("fit.ebfactanal1.", s)))
                    coef.vbfactreg(fit, outind=p + 1)$coef}), ncol=length(m),
                    dimnames=list(NULL, paste0("ebfactanal1.", 1:length(m))))
                  # ,
                  # matrix(sapply(1:length(m), function(s) {
                  #   fit <- eval(parse(text=paste0("fit.ebfactanal2.", s)))
                  #   coef.vbfactreg(fit, outind=p + 1)$coef}), ncol=length(m),
                  #   dimnames=list(NULL, paste0("ebfactanal2.", 1:length(m))))
                  )

    # calculate performance measures
    pred <- cbind(xtest %*% best, null=rep(mean(ytest), ntest))
    emse[r, ] <- c(colMeans((best - beta)^2), null=NA)
    pmse[r, ] <- colMeans((pred - ytest)^2)
    fact[r] <- est.factors
    conv[r, ] <- c(
      rep(NA, 3 + length(m)),
      sapply(1:length(m), function(s) {
        fit <- eval(parse(text=paste0("fit.emfactanal1.", s)))
        if(length(fit)==0) {
          NA
        } else {
          fit$conv
        }}),
      # sapply(1:length(m), function(s) {
      #   eval(parse(text=paste0("fit.emfactanal2.", s, "$conv")))}),
      sapply(1:length(m), function(s) {
        eval(parse(text=paste0("fit.vbfactanal1.", s, "$conv")))}),
      # sapply(1:length(m), function(s) {
      #   eval(parse(text=paste0("fit.vbfactanal2.", s, "$conv")))}),
      sapply(1:length(m), function(s) {
        eval(parse(text=paste0("fit.ebfactanal1.", s, "$conv")))}),
      # sapply(1:length(m), function(s) {
      #   eval(parse(text=paste0("fit.ebfactanal2.", s, "$conv")))}),
      NA)
    gamma[r, ] <- c(sapply(1:length(m), function(s) {
      eval(parse(text=paste0("fit.ebfactanal1.", s, "$hyper$gammagroups")))})
      # ,
      # sapply(1:length(m), function(s) {
      #   eval(parse(text=paste0("fit.ebfactanal2.", s, "$hyper$gammagroups")))})
      )

  }
  return(list(emse=emse, pmse=pmse, fact=fact, conv=conv, gamma=gamma,
              best=best))
}

# simulation function
sim2 <- function(nreps, ntest, n, m, p, d, groups, trace=TRUE, 
                 seed=NULL, simpar) {
  
  # set seed
  set.seed(seed)
  
  # setup results objects
  methods <- c("true", "ridge", "lasso",
               paste0("FMradio", 1:length(m)),
               paste0("emfactanal", 1:length(m)),
               paste0("vbfactanal", 1:length(m)),
               paste0("ebfactanal", 1:length(m)),
               "null")
  fact <- numeric(nreps)
  conv <- pcor <- pmse <- emse <-
    matrix(NA, nrow=nreps, ncol=length(methods), 
           dimnames=list(NULL, methods))
  coln <- apply(expand.grid(1:length(unique(groups)), 1:length(m)), 
                1, function(s) {paste0("ebfactanal", s[2], ".", s[1])})
  mat.gamma <- matrix(NA, nrow=nreps, ncol=length(unique(groups))*length(m),
                      dimnames=list(NULL, coln))
  
  # simulations
  for(r in 1:nreps) {
    if(trace) {cat("rep ", r, "\r")}
    
    # simulate parameters
    par <- eval(parse(text=simpar))
    Bx <- par$Bx
    By <- par$By
    psi <- par$psi
    sigma2 <- par$sigma2
    beta <- colSums(By*(solve(Bx %*% t(Bx) + diag(d)/psi) %*% Bx))*
      sqrt((colSums(Bx^2) + psi)/(sum(By^2) + sigma2))
    
    # simulate data
    Lambda <- matrix(rnorm((n + max(m))*d), nrow=n + max(m), ncol=d)
    x <- scale(t(apply(Lambda, 1, function(lam) {
      rmvnorm(1, colSums(lam*Bx), diag(psi))})))
    y <- c(scale(rnorm(n, as.numeric(Lambda %*% By), sigma2))[, 1], 
           rep(NA, max(m)))
    xy <- cbind(x, y)
    
    # fit ridge and lasso models
    fit.ridge <- cv.glmnet(scale(x[1:n, ]), y[1:n], alpha=0, intercept=FALSE)
    fit.lasso <- cv.glmnet(scale(x[1:n, ]), y[1:n], alpha=1, intercept=FALSE)
    
    # estimate models
    for(s in 1:length(m)) {
      
      # estimate penalty parameters
      fit.cormat <- bayesfactanal::regcor(xy[1:(n + m[s]), -(p + 1)])
      cormat <- fit.cormat$cor  
      est.gamma <- fit.cormat$gamma
      
      # estimate number of factors
      est.factors <- sum(eigen(cor(xy[1:(n + m[s]), -(p + 1)]),
                               only.values=TRUE)$values > 1)
      est.factors <- max(est.factors, 2)
      
      # set hyperparameters
      hyper1 <- list(kappa=rep(9, p + 1), nu=rep(4, p + 1))
      # hyper1 <- c(hyper1, 
      #             list(gamma=c((hyper1$kappa[-(p + 1)] - 1)/(hyper1$nu[-(p + 1)]*
      #                                                          est.factors),
      #                          Inf)))
      hyper1 <- c(hyper1, 
                  list(gamma=(hyper1$kappa - 1)/(hyper1$nu*est.factors)))
      hyper2 <- c(hyper1, 
                  list(gammagroups=c(1, 1)),
                  list(gammatotal=c((mean(hyper1$kappa[-(p + 1)]) - 1)/
                                      (mean(hyper1$nu[-(p + 1)])*est.factors),
                                    (hyper1$kappa[p + 1] - 1)/
                                      (hyper1$nu[p + 1]*est.factors))))
      
      # fit FMradio models
      fit.FMradio <- tryCatch({
        fit <- mlFA(cormat, est.factors);
        fit <- c(fit, cv.glmnet(as.matrix(
          facScore(scale(xy[1:n, -(p + 1)]), LM=fit$Loadings, 
                   UM=fit$Uniqueness)), scale(xy[1:n, p + 1])[, 1],
          alpha=0, intercept=FALSE));
        c(fit, list(param=list(
          B=cbind(t(fit$Loadings), 
                  coef(structure(fit, class="cv.glmnet"), 
                       s="lambda.min")[-1]),
          psi=c(diag(fit$Uniqueness), 0))))},
        error=function(e) {list()})   
      assign(paste0("fit.FMradio", s), fit.FMradio)
      
      # fit EM factanal models
      fit.emfactanal <- tryCatch({
        emfactanal(xy[1:(n + m[s]), ], factors=est.factors, 
                   gamma=est.gamma, Gamma=rep(1, p + 1))},
        error=function(e) {list()})
      assign(paste0("fit.emfactanal", s), fit.emfactanal)
      
      # fit Bayesian factanal models
      fit.vbfactanal <- bayesfactanal(x=xy[1:(n + m[s]), ], 
                                      factors=est.factors, hyper=hyper1, 
                                      method="vb", rescale=TRUE,
                                      control=list(
                                        maxit=2000,
                                        epsilon=sqrt(.Machine$double.eps), 
                                        trace=FALSE))
      assign(paste0("fit.vbfactanal", s), fit.vbfactanal)
      
      # fit empirical Bayesian factanal models
      fit.ebfactanal <- bayesfactanal(
        x=xy[1:(n + m[s]), ], factors=est.factors, groups=c(groups, 0), 
        hyper=hyper2, method="eb", rescale=TRUE, constrained=TRUE,
        control=list(maxit=2000, epsilon=sqrt(.Machine$double.eps), 
                     trace=FALSE))
      assign(paste0("fit.ebfactanal", s), fit.ebfactanal)
    }
    
    # simulate test data
    Lambdatest <- matrix(rnorm(ntest*d), nrow=ntest, ncol=d)
    xtest <- scale(t(apply(Lambdatest, 1, function(lam) {
      rmvnorm(1, colSums(lam*Bx), diag(psi))})))
    ytest <- scale(rnorm(ntest, as.numeric(Lambdatest %*% By), sigma2))[, 1]
    
    # calculate regression parameters
    best <- cbind(true=beta,
                  ridge=coef(fit.ridge, s="lambda.min")[-1], 
                  lasso=coef(fit.lasso, s="lambda.min")[-1], 
                  matrix(sapply(1:length(m), function(s) {
                    fit <- eval(parse(text=paste0("fit.FMradio", s)))
                    if(length(fit)==0) {
                      rep(NA, p)
                    } else {
                      coef.emfactreg(fit, outind=p + 1)$coef
                    }}), ncol=length(m), 
                    dimnames=list(NULL, paste0("FMradio", 1:length(m)))),
                  matrix(sapply(1:length(m), function(s) {
                    fit <- eval(parse(text=paste0("fit.emfactanal", s)))
                    if(length(fit)==0) {
                      rep(NA, p)
                    } else {
                      coef.emfactreg(fit, outind=p + 1)$coef
                    }}), ncol=length(m), 
                    dimnames=list(NULL, paste0("emfactanal", 1:length(m)))),
                  matrix(sapply(1:length(m), function(s) {
                    fit <- eval(parse(text=paste0("fit.vbfactanal", s)))
                    coef.vbfactreg(fit, outind=p + 1)$coef}), ncol=length(m), 
                    dimnames=list(NULL, paste0("vbfactanal", 1:length(m)))),
                  matrix(sapply(1:length(m), function(s) {
                    fit <- eval(parse(text=paste0("fit.ebfactanal", s)))
                    coef.vbfactreg(fit, outind=p + 1)$coef}), ncol=length(m), 
                    dimnames=list(NULL, paste0("ebfactanal", 1:length(m)))))
    
    # calculate performance measures
    pred <- xtest %*% best
    emse[r, ] <- c(colMeans((best - beta)^2), null=NA)
    pmse[r, ] <- c(colMeans((pred - ytest)^2), null=mean(ytest^2))
    pcor[r, ] <- c(apply(pred, 2, function(pr) {cor(pr, ytest)}), null=0)
    fact[r] <- est.factors
    conv[r, ] <- c(
      rep(NA, 3 + length(m)), 
      sapply(1:length(m), function(s) {
        fit <- eval(parse(text=paste0("fit.emfactanal", s)))
        if(length(fit)==0) {
          NA
        } else {
          fit$conv
        }}),
      sapply(1:length(m), function(s) {
        eval(parse(text=paste0("fit.vbfactanal", s, "$conv")))}),
      sapply(1:length(m), function(s) {
        eval(parse(text=paste0("fit.ebfactanal", s, "$conv")))}),
      NA)
    mat.gamma[r, ] <- c(sapply(1:length(m), function(s) {
      eval(parse(text=paste0("fit.ebfactanal", s, "$hyper$gammagroups")))}))
    
  }
  return(list(emse=emse, pmse=pmse, pcor=pcor, fact=fact, conv=conv, 
              gamma=mat.gamma, best=best))
}

# res=res18; m=c(0, 50, 100, 200, 500); measure="median";
# lwd=2; legend=TRUE
# EMSE and PMSE figures
simfig1 <- function(res, methods=c("true", "ridge", "lasso", "FMradio",
                                   "emfactanal", "vbfactanal", 
                                   "ebfactanal", "null"),
                    m=NULL, plotid=NULL, measure="median",
                    legend=TRUE, col=NULL, lty=NULL, labels=NULL, ...) {
  
  # set some visual
  if(is.null(col)) {
    col <- brewer.pal(n=length(methods), name="Dark2")
  }
  if(is.null(lty)) {
    lty <- rep(1, length(methods))
  }
  if(is.null(labels)) {
    labels <- methods
  }
  
  # calculate number of fits
  lengths <- sapply(methods, function(mt) {
    sum(substr(colnames(res$emse), 1, nchar(mt))==mt)})
  if(is.null(m)) {
    m <- 1:max(lengths)
  }
  if(is.null(plotid)) {
    plotid <- 1:length(m)
  }
  
  # create plot tables
  temse <- lapply(methods, function(mt) {
    vals <- apply(res$emse[, substr(colnames(res$emse), 1, nchar(mt))==mt, 
                           drop=FALSE], 2, measure, na.rm=TRUE)
    cbind(m, rep(vals, length.out=length(m)))[plotid, , drop=FALSE]})
  tpmse <- lapply(methods, function(mt) {
    vals <- apply(res$pmse[, substr(colnames(res$pmse), 1, nchar(mt))==mt, 
                           drop=FALSE], 2, measure, na.rm=TRUE)
    cbind(m, rep(vals, length.out=length(m)))[plotid, , drop=FALSE]})
  tpcor <- lapply(methods, function(mt) {
    vals <- apply(res$pcor[, substr(colnames(res$pcor), 1, nchar(mt))==mt, 
                           drop=FALSE], 2, measure, na.rm=TRUE)
    cbind(m, rep(vals, length.out=length(m)))[plotid, , drop=FALSE]})
  
  # create plots
  opar <- par(no.readonly=TRUE)
  par(mar=opar$mar*c(1, 1.3, 1, 1))
  layout(matrix(c(1:3), nrow=1, ncol=3, byrow=TRUE))
  plot(temse[[1]], type="l", main="a)", xlab="Number of unlabeled", 
       ylab="EMSE", ylim=range(sapply(temse, function(x) {x[, 2]}), na.rm=TRUE), 
       col=col[1], lty=lty[1], ...)
  for(mt in 2:length(methods)) {
    lines(temse[[mt]], col=col[mt], lty=lty[mt], ...)
  }
  plot(tpmse[[1]], type="l", main="b)", xlab="Number of unlabeled", 
       ylab="PMSE", ylim=range(sapply(tpmse, function(x) {x[, 2]}), na.rm=TRUE), 
       col=col[1], lty=lty[1], ...)
  for(mt in 2:length(methods)) {
    lines(tpmse[[mt]], col=col[mt], lty=lty[mt], ...)
  }
  plot(tpcor[[1]], type="l", main="b)", xlab="Number of unlabeled", 
       ylab=expression("corr(y,"~hat(y)), 
       ylim=range(sapply(tpcor, function(x) {x[, 2]}), na.rm=TRUE), 
       col=col[1], lty=lty[1], ...)
  for(mt in 2:length(methods)) {
    lines(tpcor[[mt]], col=col[mt], lty=lty[mt], ...)
  }
  if(legend) {
    legend("bottomright", legend=labels, col=col, lty=lty, ...)
  }
  par(opar)
}

# log multiplier figure
simfig2 <- function(res, methods=c("ebfactanal"),
                    m=NULL, measure="median",
                    legend=TRUE, col=NULL, lty=NULL, labels=NULL, ...) {
  
  # find number of groups
  G <- length(unique(sapply(strsplit(colnames(res$gamma), ".", fixed=TRUE), "[", 
                            2)))
  
  # calculate number of fits
  lengths <- sapply(methods, function(mt) {
    sum(substr(colnames(res$gamma), 1, nchar(mt))==mt)})/G
  
  # set some visual
  if(is.null(col)) {
    col <- brewer.pal(n=max(3, G), name="Dark2")
  }
  if(is.null(lty)) {
    lty <- rep(1, length(methods)) + 1
  }
  if(is.null(labels)) {
    labels <- c(methods, paste0("group ", 1:G))
  }
  
  # set number of unlabeled
  if(is.null(m)) {
    m <- 1:max(lengths)
  }
  
  # create plot tables
  tlgamma <- lapply(methods, function(mt) {
    lapply(1:G, function(g) {
      vals <- apply(
        log(res$gamma[, substr(colnames(res$gamma), 1, nchar(mt))==mt & 
                        sapply(strsplit(colnames(res$gamma), ".", fixed=TRUE), 
                               "[", 2)==as.character(g), drop=FALSE]), 
        2, measure)
      cbind(m, vals)})})
  
  # create plots
  opar <- par(no.readonly=TRUE)
  par(mar=opar$mar*c(1, 1.3, 1, 1))
  for(mt in 1:length(methods)) {
    for(g in 1:G) {
      if(mt==1 & g==1) {
        plot(tlgamma[[mt]][[g]], type="l", xlab="Number of unlabeled", 
             ylab=expression(log~hat(gamma)[g]),
             ylim=range(sapply(tlgamma, function(x) {
               sapply(x, function(y) {y[, 2]})}), na.rm=TRUE), col=col[g], 
             lty=lty[mt], ...)
      } else {
        lines(tlgamma[[mt]][[g]], col=col[g], lty=lty[mt], ...)
      }
    }
  }
  if(legend) {
    legend("topright", legend=labels, col=c(rep(1, length(methods)), col), 
           lty=c(lty, rep(1, G)), ...)
  }
  par(opar)
}



res1 <- sim1(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500),
             p=20, d=10, groups=rep(c(1, 2), each=10),
             communality=c(0.3, 0.7, 0.9), trace=TRUE, seed=NULL)
res2 <- sim1(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
             p=100, d=10, groups=rep(c(1, 2), each=50),
             communality=c(0.3, 0.7, 0.9), trace=TRUE, seed=NULL)
res3 <- sim1(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
             p=100, d=10, groups=rep(c(1, 2), each=50),
             communality=c(0.3, 0.7, 0.7), trace=TRUE, seed=NULL)
res4 <- sim1(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
             p=100, d=10, groups=rep(c(1, 2), each=50),
             communality=c(0.2, 0.5, 0.5), trace=TRUE, seed=NULL)
res5 <- sim1(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500),
             p=100, d=10, groups=rep(c(1, 2), each=50),
             communality=c(0.2, 0.8, 0.7), trace=TRUE, seed=NULL)
res6 <- sim1(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
             p=100, d=10, groups=rep(c(1, 2), each=50),
             communality=c(0.2, 0.8, 0.6), trace=TRUE, seed=NULL)
res7 <- sim1(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
             p=100, d=20, groups=rep(c(1, 2), each=50),
             communality=c(0.2, 0.8, 0.6), trace=TRUE, seed=NULL)
res8 <- sim1(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500),
             p=120, d=40, groups=rep(c(1, 2), each=60),
             communality=c(0.2, 0.8, 0.7), trace=TRUE, seed=NULL)
save(res1, res2, res3, res4, res5, res6, res7, res8,
     file="results/simulations_bayesfactreg_res2.Rdata")

simpar <- paste0("list(psi=rep(1, p), sigma2=1,",
                 "Bx=cbind(matrix(rnorm(d*p/2, 0, 0.1), nrow=d,",
                 "ncol=p/2), matrix(rnorm(d*p/2, 0, 10), nrow=d,",
                 "ncol=p/2)), By=rnorm(d))")
res9 <- sim2(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), p=100, d=40, 
             groups=rep(c(1, 2), each=50), trace=TRUE, 
             seed=NULL, simpar=simpar)
res10 <- sim2(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), p=100, 
              d=30, groups=rep(c(1, 2), each=50), trace=TRUE, 
              seed=NULL, simpar=simpar)
simpar <- paste0("list(psi=rep(1, p), sigma2=1,",
                 "Bx=cbind(rbind(matrix(rnorm(d*p/4, 0, 0.1), nrow=d/2,", 
                 "ncol=p/2), matrix(0, nrow=d/2, ncol=p/2)),", 
                 "rbind(matrix(0, nrow=d/2, ncol=p/2),", 
                 "matrix(rnorm(d*p/4, 0, 10), nrow=d/2, ncol=p/2))),", 
                 "By=rnorm(d))")
res11 <- sim2(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), p=100, 
              d=30, groups=rep(c(1, 2), each=50), trace=TRUE, 
              seed=NULL, simpar=simpar)
simpar <- paste0("list(psi=rep(1, p), sigma2=1,",
                 "Bx=cbind(rbind(matrix(rnorm(d*p/4, 0, 0.1), nrow=d/2,", 
                 "ncol=p/2), matrix(0, nrow=d/2, ncol=p/2)),", 
                 "rbind(matrix(0, nrow=d/2, ncol=p/2),", 
                 "matrix(rnorm(d*p/4, 0, 10), nrow=d/2, ncol=p/2))),", 
                 "By=c(rnorm(d/2), rep(0, d/2)))")
res12 <- sim2(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), p=100, 
              d=30, groups=rep(c(1, 2), each=50), trace=TRUE, 
              seed=NULL, simpar=simpar)
simpar <- paste0("list(psi=rep(1, p), sigma2=1,", 
                 "Bx=t(model.matrix(~ -1 + factor(rep(1:d, each=p/d)))*",
                 "sqrt(rep(c(0.2, 0.8), each=p/2)/",
                 "(rowSums(model.matrix(~ -1 + factor(rep(1:d, each=p/d))))*",
                 "(1 - rep(c(0.2, 0.8), each=p/2))))),",
                 "By=rep(sqrt(0.7/(d*(1 - 0.7))), d))")
res13 <- sim2(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
              p=100, d=10, groups=rep(c(1, 2), each=50), 
              trace=TRUE, seed=NULL, simpar=simpar)
simpar <- paste0("list(psi=rep(1, p), sigma2=1,", 
                 "Bx=t(model.matrix(~ -1 + factor(rep(1:d, each=p/d)))*",
                 "sqrt(rep(c(0.3, 0.7), each=p/2)/",
                 "rowSums(model.matrix(~ -1 + factor(rep(1:d, each=p/d))))*",
                 "(1 - rep(c(0.3, 0.7), each=p/2/)))),",
                 "By=rep(sqrt(0.9/(d*(1 - 0.9))), d))")
res14 <- sim2(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
              p=20, d=10, groups=rep(c(1, 2), each=10), 
              trace=TRUE, seed=NULL, simpar=simpar)
simpar <- paste0("list(psi=rep(1, p), sigma2=1,", 
                 "Bx=t(model.matrix(~ -1 + factor(rep(1:d, each=p/d)))*",
                 "sqrt(rep(c(0.2, 0.8), each=p/2)/",
                 "(rowSums(model.matrix(~ -1 + factor(rep(1:d, each=p/d))))*",
                 "(1 - rep(c(0.2, 0.8), each=p/2))))),",
                 "By=rep(sqrt(0.6/(d*(1 - 0.6))), d))")
res15 <- sim2(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
              p=100, d=20, groups=rep(c(1, 2), each=50),
              trace=TRUE, seed=NULL, simpar=simpar)
simpar <- paste0("list(psi=rep(1, p), sigma2=1,", 
                 "Bx=t(model.matrix(~ -1 + factor(rep(1:d, each=p/d)))*",
                 "sqrt(rep(c(0.2, 0.8), each=p/2)/",
                 "(rowSums(model.matrix(~ -1 + factor(rep(1:d, each=p/d))))*",
                 "(1 - rep(c(0.2, 0.8), each=p/2))))),",
                 "By=rep(sqrt(0.9/(d*(1 - 0.9))), d))")
res16 <- sim2(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
              p=100, d=10, groups=rep(c(1, 2), each=50),
              trace=TRUE, seed=NULL, simpar=simpar)
simpar <- paste0("list(psi=rep(1, p), sigma2=1,", 
                 "Bx=t(model.matrix(~ -1 + factor(rep(1:d, each=p/d)))*",
                 "sqrt(rep(c(0.3, 0.7), each=p/2)/",
                 "(rowSums(model.matrix(~ -1 + factor(rep(1:d, each=p/d))))*",
                 "(1 - rep(c(0.3, 0.7), each=p/2))))),",
                 "By=rep(sqrt(0.5/(d*(1 - 0.5))), d))")
res17 <- sim2(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
              p=100, d=10, groups=rep(c(1, 2), each=50),
              trace=TRUE, seed=NULL, simpar=simpar)
simpar <- paste0("Bx <- Reduce('cbind', lapply(1:d, function(h)",
                 "{m <- matrix(0, nrow=d, ncol=p/d);",
                 "m[c(1:d, 1)[h:(h + 1)], ] <- 1; m}));",
                 "list(psi=rep(1, p), sigma2=1,", 
                 "Bx=t(t(Bx)*sqrt(rep(c(0.3, 0.7), each=p/2)/",
                 "(colSums(Bx)*(1 - rep(c(0.2, 0.8), each=p/2))))),",
                 "By=rep(sqrt(0.8/(d*(1 - 0.8))), d))")
res18 <- sim2(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
              p=100, d=10, groups=rep(c(1, 2), each=50),
              trace=TRUE, seed=NULL, simpar=simpar)
simpar <- paste0("Bx <- Reduce('cbind', lapply(1:d, function(h)",
                 "{m <- matrix(0, nrow=d, ncol=p/d);",
                 "m[c(1:d, 1)[h:(h + 1)], ] <- 1; m}));",
                 "list(psi=rep(1, p), sigma2=1,", 
                 "Bx=cbind(Bx[, 1:(p/2)]*rnorm(d*p/2, 0, 0.1),", 
                 "Bx[, (p/2 + 1):p]*rnorm(d*p/2, 0, 1)), ",
                 "By=rep(sqrt(0.7/(d*(1 - 0.7))), d))")
res19 <- sim2(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
              p=100, d=10, groups=rep(c(1, 2), each=50),
              trace=TRUE, seed=NULL, simpar=simpar)
simpar <- paste0("Bx <- Reduce('cbind', lapply(1:d, function(h)",
                 "{m <- matrix(0, nrow=d, ncol=p/d);",
                 "m[c(1:d, 1)[h:(h + 1)], ] <- 1; m}));",
                 "list(psi=rep(1, p), sigma2=1,", 
                 "Bx=cbind(Bx[, 1:(p/2)]*rnorm(d*p/2, 0, 0.1),", 
                 "Bx[, (p/2 + 1):p]*rnorm(d*p/2, 0, 1)), ",
                 "By=rnorm(d))")
res20 <- sim2(nreps=50, ntest=1000, n=50, m=c(0, 50, 100, 200, 500), 
              p=100, d=10, groups=rep(c(1, 2), each=50),
              trace=TRUE, seed=NULL, simpar=simpar)

save(res9, res10, res11, res12, res13, res14, res15, res16, res17, res18,
     res19, res20, file="results/simulations_bayesfactreg_res3.Rdata")
# load(file="results/simulations_bayesfactreg_res3.Rdata")

simfig1(res9, m=c(0, 50, 100, 200, 500), measure="median",
        lwd=2, legend=TRUE)
simfig1(res19, m=c(0, 50, 100, 200, 500), measure="median",
        lwd=2, legend=TRUE)
