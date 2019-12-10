### testing FusedLasso
library(FusedLasso)
library(glmnet)
library(mvtnorm)

family="binomial"
nlambda=10
addIntercept=FALSE
seed=NULL
lambda1=NULL
lambda2=NULL
exp.grid=TRUE
nfolds=10
foldid=NULL
trace=FALSE
lambda.min.ratio=ifelse(nrow(x) < ncol(x), 0.01, 1e-04)
cv.fusedlasso <- function(x, y, graph, family, lambda1=NULL, lambda2=NULL, 
                          nlambda=100, exp.grid=TRUE, addIntercept=TRUE, 
                          nfolds=10, foldid=NULL, seed=NULL, trace=FALSE, 
                          lambda.min.ratio=ifelse(nrow(x) < ncol(x), 0.01, 
                                                  1e-04)) {
  
  if(is.null(lambda1) | is.null(lambda2)) {
    invisible(capture.output(
      fit <- fusedlassoMaxLambdas(x, y, family=family, graph=graph, 
                                  addIntercept=addIntercept)))
    lim.lambda1 <- c(fit$maxLambda1, lambda.min.ratio*fit$maxLambda1)
    lim.lambda2 <- c(fit$maxLambda2, lambda.min.ratio*fit$maxLambda2)
    if(is.null(lambda1)) {
      if(!is.null(lambda2) & length(lambda2)!=1) {
        lambda1 <- exp(seq(log(lim.lambda1[1]), log(lim.lambda1[2]), 
                           length.out=length(lambda2)))
      } else {
        lambda1 <- exp(seq(log(lim.lambda1[1]), log(lim.lambda1[2]), 
                           length.out=nlambda))
      }
    }
    if(is.null(lambda2)) {
      if(!is.null(lambda1) & length(lambda1)!=1) {
        lambda2 <- exp(seq(log(lim.lambda2[1]), log(lim.lambda2[2]), 
                           length.out=length(lambda1)))
      } else {
        lambda2 <- exp(seq(log(lim.lambda2[1]), log(lim.lambda2[2]), 
                           length.out=nlambda))
      }
    }
  }
  if(exp.grid) {
    grid <- expand.grid(lambda1=lambda1, lambda2=lambda2)
    seq.lambda1 <- grid$lambda1
    seq.lambda2 <- grid$lambda2
  }
  
  n <- nrow(x)
  p <- ncol(x)
  if(is.null(foldid)) {
    if(is.null(seed)) {
      seed <- Sys.time()
    }
    set.seed(seed)
    if(is.null(nfolds)) {
      nfolds <- n
    }
    remainder <- n %% nfolds
    if(remainder==0) {
      foldid <- rep(1:nfolds, each=n/nfolds)
    } else {
      foldid <- c(rep(1:remainder, each=n %/% nfolds + 1),
                  rep((remainder + 1):nfolds, each=n %/% nfolds))
    }
    foldid <- sample(foldid)
  }
  
  cvl <- matrix(NA, nrow=nfolds, ncol=length(lambda1)*length(lambda1))
  for(fold in 1:nfolds) {
    if(trace) {
      cat("\r", "fold ", fold, sep="")
    }
    xtrain <- as.matrix(x[foldid!=fold, ], ncol=p)
    xtest <- as.matrix(x[foldid==fold, ], ncol=p)
    ytrain <- y[foldid!=fold]
    ytest <- y[foldid==fold]
    for(l2 in 1:length(lambda2)) {
      invisible(capture.output(
        fit <- fusedlasso(xtrain, ytrain, lambda1=lambda1, lambda2=lambda2[l2], 
                          family=family, graph=graph, addIntercept=addIntercept,
                          verbose=FALSE)))
      
      best <- as.matrix(rbind(fit$Intercept[1:ncol(fit$beta)], fit$beta))
      lpred <- cbind(1, xtest) %*% best
      if(family=="gaussian") {
        cvl[fold, ((l2 - 1)*length(lambda1) + 1):((l2 - 1)*length(lambda1) + 
                                                    length(fit$lambda1))] <-  
          0.5*colSums((ytest - lpred)^2)
      } else if(family=="binomial") {
        cvl[fold, ((l2 - 1)*length(lambda1) + 1):((l2 - 1)*length(lambda1) + 
                                                    length(fit$lambda1))] <- 
          -colSums(ytest*lpred - log(1 + exp(lpred)))
      }
    }
  }
  cvm <- colMeans(cvl, na.rm=TRUE)
  cvsd <- apply(cvl, 2, sd, na.rm=TRUE)
  lambda.min <- c(seq.lambda1[which.min(cvm)], seq.lambda2[which.min(cvm)])
  fit.list <- lapply(1:length(lambda2), function(l2) {
    invisible(capture.output(
      fit <- fusedlasso(x, y, lambda1=lambda1, lambda2=lambda2[l2], 
                        family=family, graph=graph, addIntercept=addIntercept, 
                        verbose=FALSE))); return(fit)})
  
  fit <- list(beta=Reduce("cbind", lapply(fit.list, "[[", "beta")), 
              success=unlist(sapply(fit.list, "[[", "success")), 
              lambda1=unlist(sapply(fit.list, "[[", "lambda1")), 
              lambda2=unlist(sapply(fit.list, "[[", "lambda2")),
              outerIterNum=unlist(sapply(fit.list, "[[", "outerIterNum")), 
              innerIterNum=unlist(sapply(fit.list, "[[", "innerIterNum")),
              family=family, 
              Intercept=unlist(sapply(1:length(fit.list), function(l) {
                fit.list[[l]]$Intercept[1:length(fit.list[[l]]$lambda1)]})), 
              maxLambda1=lim.lambda1[1],
              maxLambda2=lim.lambda2[1])
  class(fit) <- "fusedlasso"
  
  nzero <- colSums(fit$beta!=0) + as.numeric(fit$Intercept!=0)
  
  out <- list(lambda1=fit$lambda1, lambda2=fit$lambda2, cvm=cvm, cvsd=cvsd, 
              nzero=nzero, fit=fit, lambda.min=lambda.min)
  return(out)
  
}


# set.seed(2019)
n <- 200
p <- 100
G <- 2
ne <- choose(p, 2)/5
A <- matrix(0, p, p)
A[upper.tri(A)] <- sample(rep(c(1, 0), times=c(ne, p*(p-1)/2 - ne)))
A[lower.tri(A)] = t(A)[lower.tri(A)]
graph <- list(apply(A, 1, function(a) {which(a!=0)}), 
              apply(A, 1, function(a) {rep(1, sum(a))}))
rho <- 0.4
Sigma <- A*rho
diag(Sigma) <- rowSums(Sigma) + 1
Sigma <- cov2cor(Sigma)
beta <- rmvnorm(1, rep(0, p), Sigma)
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
lpred <- as.numeric(x %*% t(beta))
y <- rbinom(n, 1, exp(lpred)/(1 + exp(lpred)))

fit.fusedlasso1 <- cv.fusedlasso(x, y, graph, family="binomial", 
                                 nlambda=10, addIntercept=FALSE, seed=NULL)
fit.lasso1 <- cv.glmnet(x, y, family="binomial", alpha=1, intercept=FALSE,
                        standardize=FALSE)
fit.ridge1 <- cv.glmnet(x, y, family="binomial", alpha=0, intercept=FALSE,
                        standardize=FALSE)

best.fusedlasso1 <- fit.fusedlasso1$fit$beta[, which(
  fit.fusedlasso1$lambda1==fit.fusedlasso1$lambda.min[1] &
    fit.fusedlasso1$lambda2==fit.fusedlasso1$lambda.min[2])]
best.lasso1 <- as.numeric(coef(fit.lasso1, s="lambda.min"))[-1]
best.ridge1 <- as.numeric(coef(fit.ridge1, s="lambda.min"))[-1]

nzero.fusedlasso1 <- fit.fusedlasso1$nzero[which(
  fit.fusedlasso1$lambda1==fit.fusedlasso1$lambda.min[1] &
    fit.fusedlasso1$lambda2==fit.fusedlasso1$lambda.min[2])]
nzero.lasso1 <- unname(fit.lasso1$nzero[which(
  fit.lasso1$lambda==fit.lasso1$lambda.min)])
nzero.ridge1 <- unname(fit.ridge1$nzero[which(
  fit.ridge1$lambda==fit.ridge1$lambda.min)])

nzero <- c(fusedlasso1=nzero.fusedlasso1, lasso1=nzero.lasso1, 
           ridge1=nzero.ridge1)
emse <- c(fusedlasso1=mean((beta - best.fusedlasso1)^2), 
          lasso1=mean((beta - best.lasso1)^2),
          ridge1=mean((beta - best.ridge1)^2))
cbind(nzero, emse)

plot(beta, best.fusedlasso1)
plot(beta, best.lasso1)
plot(beta, best.ridge1)

fit.fusedlasso1$lambda.min
fit.lasso1$lambda.min*n


#### variational bayes posterior linear regression with inverse gamma prior
n <- 100
p <- 10
beta <- c(1:p)
sigma2 <- 1
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
y <- as.numeric(x %*% beta) + rnorm(n, 0, sqrt(sigma2))

vbpost <- function(x, y, a0, b0, lambda2, eps=0.0001, maxiter=100) {
  p <- ncol(x)
  fit.start <- lm(y ~ -1 + x)
  A <- t(x) %*% x + lambda2*diag(p)
  Ainv <- solve(A)
  const <- t(mu) %*% A %*% mu + t(y) %*% y
  
  mu <- Ainv %*% t(x) %*% y
  a <- p/2 + n/2 + a0
  
  sigmaold <- sigma <- var(residuals(fit.start))*Ainv
  bold <- b <- as.numeric(b0 + 0.5*(sum(diag(A %*% sigma)) - 
                                      2*t(y) %*% x %*% mu + const))
  
  iter <- 0
  conv <- FALSE
  while(!conv & (iter <= maxiter)) {
    iter <- iter + 1
    sigma <- bold*Ainv/a
    b <- as.numeric(b0 + 0.5*(sum(diag(A %*% sigma)) - 2*t(y) %*% x %*% mu + 
                                const))
    
    conv <- max(abs(c(b - bold, sigma - sigmaold))) < eps
    
    bold <- b
    sigmaold <- sigma
    
  }
  
  return(list(mu=mu, sigma=sigma, a=a, b=b, niter=iter))
  
}

mcmcpost <- function(x, y, a0, b0, lambda2, K) {
  p <- ncol(x)
  n <- nrow(x)
  A <- t(x) %*% x + lambda2*diag(p)
  Ainv <- solve(A)
  const <- t(y) %*% y
  
  mu <- Ainv %*% t(x) %*% y
  a <- n/2 + p/2 + a0
  
  sigma2 <- var(residuals(lm(y ~ -1 + x)))
  
  sigma2mat <- numeric(K)
  betamat <- matrix(ncol=p, nrow=K)
  k <- 1
  while(k <= K) {
    beta <- t(rmvnorm(1, mu, sigma2*Ainv))
    sigma2 <- rgamma(1, shape=a, rate=b0 + 0.5*(const - 2*t(y) %*% x %*% beta + 
                                                  t(beta) %*% A %*% beta))
    betamat[k, ] <- beta
    sigma2mat[k] <- sigma2
    k <- k + 1
  }
  
  return(list(beta=betamat, sigma2=sigma2mat))
}


vbpost1 <- vbpost(x, y, a0=1, b0=1, lambda2=1)
mcmcpost1 <- mcmcpost(x, y, a0=1, b0=1, lambda2=1, K=10000)



hist(mcmcpost1$beta[, 1], freq=FALSE, breaks=20)
curve(dnorm(x, mean=vbpost1$mu[1], sd=sqrt(vbpost1$sigma[1, 1])), add=TRUE)
curve(dnorm(x, mean=vbpost1$mu[1], sd=sqrt(vbpost1$sigma[1, 1])), from=0, 2)

mean(mcmcpost1$beta[, 1])
vbpost1$mu[1]
sqrt(vbpost1$sigma[1, 1])
sd(mcmcpost1$beta[, 1])


#### network propagation for prediction
library(mvtnorm)
set.seed(123)
n <- 100
p <- 50
rho <- 0.7
# Sigma <- matrix(rbinom(p^2, 1, 0.2), ncol=p, nrow=p)*rho
# Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
Sigma <- matrix(rho, ncol=p, nrow=p)
diag(Sigma) <- 1
x <- rmvnorm(n, rep(0, p), Sigma)
beta <- rnorm(p)
y <- as.numeric(x %*% beta + rnorm(n))

seq.lambda <- seq(0.01, 2, length.out=100)
seq.cvll1 <- seq.cvll2 <- numeric(length(seq.lambda))
for(l in 1:length(seq.lambda)) {
  for(i in 1:n) {
    xtrain <- x[-i, ]
    xtest <- matrix(x[i, ], ncol=p)
    ytrain <- y[-i]
    ytest <- y[i]

    pred1 <- xtest %*% solve(cov(xtrain) + seq.lambda[l]*diag(p)) %*% 
      apply(x, 2, function(j) {cov(j, y)})
    pred2 <- xtest %*% solve(t(xtrain) %*% xtrain + seq.lambda[l]*diag(p)) %*% 
      t(xtrain) %*% ytrain
    
    seq.cvll1[l] <- seq.cvll1[l] + mean((ytest - pred1)^2)
    seq.cvll2[l] <- seq.cvll2[l] + mean((ytest - pred2)^2)
    
  }
}

set.seed(213)
ntest <- 1000
xtest <- rmvnorm(ntest, rep(0, p), Sigma)
ytest <- as.numeric(xtest %*% beta + rnorm(n))
pred1 <- xtest %*% solve(cov(x) + seq.lambda[which.min(seq.cvll1)]*diag(p)) %*% 
  apply(x, 2, function(j) {cov(j, y)})
pred2 <- xtest %*% solve(t(x) %*% x + seq.lambda[which.min(seq.cvll2)]*
                           diag(p)) %*% t(x) %*% y
pred3 <- xtest %*% solve(cov(x) + seq.lambda[seq.lambda >= 0.5][
  which.min(seq.cvll1[seq.lambda >=0.5 ])]*diag(p)) %*% 
  apply(x, 2, function(j) {cov(j, y)})

mse0 <- mean((ytest - xtest %*% coef(lm(y ~ -1 + x)))^2)
mse1 <- mean((ytest - pred1)^2)
mse2 <- mean((ytest - pred2)^2)
mse3 <- mean((ytest - pred3)^2)
c(mse0, mse1, mse2, mse3)

plot(ytest, pred1)
plot(ytest, pred2)
plot(ytest, pred3)





n <- 100
p <- 40
rho <- 0.7
Sigma <- matrix(0.2, ncol=p, nrow=p); diag(Sigma) <- 1
x <- rmvnorm(n, rep(0, p), Sigma)
beta <- (1:p - mean(1:p))
y <- x %*% beta + rnorm(n)

xs <- apply(x, 2, function(j) {(j - mean(j))/sd(j)})
ys <- y - mean(y)

sxy <- apply(x, 2, function(j) {cov(j, y)})*(n-1)/n
sy <- sd(y)*(n-1)/n
Sx <- cov(x)*(n-1)/n

cbind(coef(lm(ys ~ -1 + xs)), coef(lm(y ~ -1 + x))*sqrt(diag(Sx)),
      solve(cov(xs)) %*% apply(xs, 2, function(j) {cov(j, ys)}))

cbind((t(xs) %*% ys)/(n - 1), apply(xs, 2, function(j) {cov(j, ys)}))

set.seed(1)
L <- matrix(rbinom())
Tm <- diag(rchisq(10, 2))
Tm %*% solve(Tm - diag(10))

### Bayesian formulation
p <- 10
ne <- 30
A <- matrix(0, p, p)
A[upper.tri(A)] <- sample(rep(c(1, 0), times=c(ne, p*(p-1)/2 - ne)))
A[lower.tri(A)] = t(A)[lower.tri(A)]
D <- diag(rowSums(A))
L <- D - A
Lnorm <- diag(p) - sqrt(solve(D)) %*% A %*% sqrt(solve(D))
S <- diag(rchisq(p, 1) + 1)
eigL <- eigen(L)
QL <- eigL$vectors
ML <- diag(eigL$values)
eigLnorm <- eigen(Lnorm)
MLnorm <- eigLnorm$vectors
QLnorm <- diag(eigLnorm$values)

t(MLnorm)
MLnorm
solve(MLnorm)

det(S + QL %*% ML %*% t(QL))
det(QL %*% (S + ML) %*% t(QL))
prod(diag(S) + diag(ML))

test <- rnorm(100)
sum(abs(test[1] - test[-1]))
abs(sum(test[1] - test[-1]))

x <- rnorm(10)
y <- rnorm(10)
cbind(x, y, (abs(x - y) - abs(x + y))/2)
cbind(x, y, (- abs(x - y) + abs(x + y))/2)
cbind(x, y, apply(cbind(x, y), 1, function(r) {min(abs(r[1]), abs(r[2]))}))
cbind(x, y, apply(cbind(x, y), 1, function(r) {max(-r[1], -r[2])}))
cbind(abs(x) + abs(y) + abs(abs(x) - abs(y)), abs(x - y) - abs(x + y))




### variational bayes network lasso (with intercept)
vbnetlasso <- function(x, y, A, lambda1, lambda2, epsilon=0.001) {
  
  n <- length(y)
  p <- ncol(x)
  m <- rep(1, n)
  alpha <- lambda1/(2*lambda2 + lambda1)
  lambda <- (lambda2 + 0.5*lambda1)/n
  
  fit.start <- glmnet(x=x, y=y, family="binomial", alpha=0, lambda=lambda,
                      standardize=FALSE, intercept=intercept)
  pred.start <- as.numeric(predict(fit.start, newx=x, type="response"))
  Xw <- diag(sqrt(pred.start*(1 - pred.start))) %*% x
  inv.mat <- solve(t(Xw) %*% Xw + 2*lambda*diag(p)*lambda)
  sigmaold <- inv.mat %*% t(Xw) %*% Xw %*% inv.mat
  muold <- as.numeric(sigmaold %*% t(x) %*% as.matrix(y - 0.5*m))
  ciold <- sqrt(as.numeric(rowSums((x %*% sigmaold)*x)) + 
                  as.numeric(x %*% muold)^2)
  chiold <- as.numeric(diag(sigmaold) + muold^2)
  Psiold <- matrix(0, p, p)
  for(j in 2:p) {
    for(k in 1:(j - 1)) {
      if(A[j, k]!=0) {
        Psiold[j, k] <- -1/sqrt(sigmaold[j, j] + sigmaold[k, k] + muold[j]^2 + 
                                  muold[k]^2 - 2*sigmaold[j, k] - 
                                  2*muold[j]*muold[k])
      }
    }
  }
  Psiold[upper.tri(Psiold)] <- t(Psiold)[upper.tri(Psiold)]
  diag(Psiold) <- -rowSums(Psiold)
    
  Psi <- matrix(0, p, p)
  niter <- 0
  conv <- FALSE
  while(!conv) {
    niter <- niter + 1
    sigma <- solve(t(x) %*% diag(0.5*m/ciold*tanh(0.5*ciold)) %*% x +
                     lambda1*diag(1/sqrt(chiold)) + lambda2*Psiold/sqrt(2))
    mu <- sigma %*% t(x) %*% (y - 0.5*m)
    ci <- sqrt(as.numeric(rowSums((x %*% sigma)*x)) + as.numeric(x %*% mu)^2)
    chi <- as.numeric(diag(sigma) + mu^2)
    for(j in 2:p) {
      for(k in 1:(j - 1)) {
        if(A[j, k]!=0) {
          Psi[j, k] <- -1/sqrt(sigma[j, j] + sigma[k, k] + mu[j]^2 + mu[k]^2 - 
                                 2*sigma[j, k] - 2*mu[j]*mu[k])
        }
      }
    }
    Psi[upper.tri(Psi)] <- t(Psi)[upper.tri(Psi)]
    diag(Psi) <- 0
    diag(Psi) <- -rowSums(Psi)
    
    conv <- max(c(abs((mu - muold)/ifelse(muold==0, muold + 0.00001, muold)), 
                  abs((diag(sigma) - diag(sigmaold))/
                        ifelse(diag(sigmaold)==0, diag(sigmaold) + 0.00001, 
                               diag(sigmaold))))) < epsilon
   
    muold <- mu
    sigmaold <- sigma
    ciold <- ci
    chiold <- chi
    Psiold <- Psi
     
  }
  
  psi <- Psi
  diag(psi) <- 0
  psi[A!=0] <- 1/psi[A!=0]^2
  return(list(mu=mu, sigma=sigma, ci=ci, chi=chi, psi=psi))
  
}

library(glmnet)
n <- 40
p <- 100
x <- matrix(rnorm(p*n), nrow=n, ncol=p)
beta <- rnorm(p)
prob <- exp(x %*% beta)/(1 + exp(x %*% beta))
y <- rbinom(n, 1, prob)
lambda1 <- 0.01
lambda2 <- 1
intercept <- FALSE
A <- matrix(1, p, p)
diag(A) <- 0
res <- vbnetlasso(x, y, A, lambda1, lambda2, epsilon=0.001)
hist(res$chi)
plot(beta, res$mu)


heatmap(res$sigma, Rowv=NA, Colv=NA, revC=TRUE, labRow=FALSE, labCol=FALSE)
hist(res$psi[lower.tri(res$psi)])

f <- function(x, c1, c2, c3) {
  exp(-c1*abs(x) - c2*x^2 + c3*x)
    
}
c1 <- rchisq(1, 1)
c2 <- rchisq(1, 1)
c3 <- rchisq(1, 1)
curve(f(x, 100, c2, c3), from=-10, to=10, n=1000)





