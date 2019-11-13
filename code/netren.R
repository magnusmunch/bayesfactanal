

library(Rcpp)
sourceCpp("/Users/magnusmunch/Documents/OneDrive/PhD/netren/code/est_param.cpp")
source("/Users/magnusmunch/Documents/OneDrive/PhD/netren/code/fopt_groups.R")


netren <- function(x, y, A, m=rep(1, nrow(x)), alpha=c(0.5, 0.5), lambda=NULL, 
                   gamma=NULL) {
  
  
  # auxiliary variables
  igroup <- rep(1:G, times=sizes)
  
  # fixed variables
  kappa <- y - 0.5*m
  L <- -A
  diag(L) <- diag(L) + rowSums(A)
  
  
  
  # starting values
  ciold <- rchisq(n, 1)
  etaold <- rchisq(p, 1)
  Zold <- matrix(0, p, p)
  Zold[upper.tri(A) & A==1] <- -rchisq(sum(A)/2, 1)
  Zold[lower.tri(Zold)] <- t(Zold)[lower.tri(Zold)]
  diag(Zold) <- -rowSums(Zold)
  
  # iterative algorithm
  conv <- FALSE
  while(!conv & (iter1 < control$maxit)) {
    iter1 <- iter1 + 1
    
    
    dg <- sapply(1:G, function(g) {
      sum((diag(sigma)[igroup==g] + mu[igroup==g]^2)*
            (1 + lambda1/sqrt(lambda2*eta[igroup==g])))})
    
    hg <- sapply(1:G, function(g) {
      cA <- A[igroup==g, igroup==g]
      cmu <- as.numeric(mu)[igroup==g]
      csigma <- sigma[igroup==g, igroup==g] 
      czeta <- zeta[igroup==g, igroup==g, drop=FALSE]
      czeta[which(czeta!=0)] <- (1 + gamma1/sqrt(gamma2*czeta[czeta!=0]))
      sum((cA*cmu + t(cA*cmu) - 2*t(cA*cmu)*cmu + cA*diag(csigma) + 
             t(cA*diag(csigma)) - 2*csigma*(cA!=0))*czeta)})
    
    esizes <- sapply(1:G, function(g) {
      sum(A[igroup==g, igroup==g])/2})
    if(u==FALSE) {
      opt <- optim(par=c(s, log(lambdag)), fn=fopt_groups1, sizes=sizes, 
                   sum1=lambda2*dg, control=list(maxit=1000), method="BFGS")
      lambdagnew <- exp(opt$par[-1])
      gammagnew <- exp(sum(esizes*(log(hg) - log(esizes)))/sum(esizes))*
        esizes/hg
      
    } else {
      if(u==TRUE) {
        u <- 0.5
      }
      opt <- optim(par=c(s, log(lambdag)), fn=fopt_groups2, 
                   sizes=list(sizes, esizes), sum1=lambda2*dg + u*gamma2*hg, 
                   u=u, control=list(maxit=5000), method="BFGS")
      lambdagnew <- exp(opt$par[-c(1, 2)])
      gammagnew <- u*lambdagnew
    }
  }
  
}



### testing
library(FusedLasso)
library(glmnet)
library(pROC)
library(igraph)
library(mvtnorm)
library(Matrix)
library(R.utils)

# simulation 1
n <- 100
ntest <- 1000
p <- 12
rho <- 0.5

beta <- c(rep(0.3, 3), rep(0, p - 3))
A1 <- matrix(0, p, p)
A1[c(4:6), 1] <- A1[c(7:9), 2] <- A1[c(10:12), 3] <- 1
A1[upper.tri(A1)] <- t(A1)[upper.tri(A1)]
rownames(A1) <- colnames(A1) <- paste("X", c(1:12), sep="")

graph.list1 <- list(apply(A1, 1, function(a) {which(a!=0)}), 
                    apply(A1, 1, function(a) {rep(1, sum(a))}))

alpha <- 0.9
seq.lambda <- exp(seq(log(0.30), log(0.001), length.out=50))

set.seed(2018)
nreps <- 100
psel1.1 <- psel1.2 <- auc1.1 <- auc1.2 <- vector(mode="list", length=nreps)

for(rep in 1:nreps) {
  x <- matrix(0, n, p)
  x[, c(4:12)] <- rnorm(9*n)
  x[, 1] <- rowSums(rho*x[, c(4:6)]) + rnorm(n, 0, sqrt(1 - 3*rho^2))
  x[, 2] <- rowSums(rho*x[, c(7:9)]) + rnorm(n, 0, sqrt(1 - 3*rho^2))
  x[, 3] <- rowSums(rho*x[, c(10:12)]) + rnorm(n, 0, sqrt(1 - 3*rho^2))
  y <- rbinom(n, 1, 1/(1 + exp(-x %*% beta)))
  
  xtest <- matrix(0, ntest, p)
  xtest[, c(4:12)] <- rnorm(9*ntest)
  xtest[, 1] <- rowSums(rho*xtest[, c(4:6)]) + 
    rnorm(ntest, 0, sqrt(1 - 3*rho^2))
  xtest[, 2] <- rowSums(rho*xtest[, c(7:9)]) + 
    rnorm(ntest, 0, sqrt(1 - 3*rho^2))
  xtest[, 3] <- rowSums(rho*xtest[, c(10:12)]) + 
    rnorm(ntest, 0, sqrt(1 - 3*rho^2))
  ytest <- rbinom(ntest, 1, 1/(1 + exp(-xtest %*% beta)))
  
  test1.1 <- glmnet(x, y, alpha=1, family="binomial", standardize=FALSE)
  pred1.1 <- predict(test1.1, xtest)
  psel1.1[[rep]] <- test1.1$df
  auc1.1[[rep]] <- apply(pred1.1, 2, function(pred) {
    pROC::roc(ytest, pred)$auc})
  
  psel1.2[[rep]] <- auc1.2[[rep]] <- numeric(length(seq.lambda))
  for(l in 1:length(seq.lambda)) {
    test1.2 <- fusedlasso(x, y, lambda1=seq.lambda[l]*alpha, 
                          lambda2=seq.lambda[l]*(1 - alpha), family="binomial", 
                          graph=graph.list1)
    pred1.2 <- as.numeric(1/(1 + exp(-cbind(1, xtest) %*% 
                                       rbind(test1.2$Intercept, test1.2$beta))))
    psel1.2[[rep]][l] <- sum(test1.2$beta!=0)
    auc1.2[[rep]][l] <- pROC::roc(ytest, pred1.2)$auc
  }
}

aucd1.1 <- cbind(Reduce(c, psel1.1), Reduce(c, auc1.1))
aucm1.1 <- sapply(sort(unique(aucd1.1[, 1])), function(psel) {
  mean(aucd1.1[aucd1.1[, 1]==psel, 2], na.rm=TRUE)})
aucv1.1 <- sapply(sort(unique(aucd1.1[, 1])), function(psel) {
  var(aucd1.1[aucd1.1[, 1]==psel, 2], na.rm=TRUE)})
pselm1.1 <- sort(unique(aucd1.1[, 1]))

aucd1.2 <- cbind(Reduce(c, psel1.2), Reduce(c, auc1.2))
aucm1.2 <- sapply(sort(unique(aucd1.2[, 1])), function(psel) {
  mean(aucd1.2[aucd1.2[, 1]==psel, 2], na.rm=TRUE)})
aucv1.2 <- sapply(sort(unique(aucd1.2[, 1])), function(psel) {
  var(aucd1.2[aucd1.2[, 1]==psel, 2], na.rm=TRUE)})
pselm1.2 <- sort(unique(aucd1.2[, 1]))

par(mfrow=c(1, 2))
xlim <- range(pselm1.1, pselm1.2)
ylim <- range(aucm1.1 - aucv1.1, aucm1.1 + aucv1.1, aucm1.2 - aucv1.2, 
              aucm1.2 + aucv1.2)
plot(pselm1.1, aucm1.1, type="l", xlim=xlim, ylim=ylim, 
     xlab="Number of selected variables", ylab="AUC")
arrows(pselm1.1[aucv1.1!=0], (aucm1.1 - aucv1.1)[aucv1.1!=0], 
       pselm1.1[aucv1.1!=0], (aucm1.1 + aucv1.1)[aucv1.1!=0], length=0.05, 
       angle=90, code=3)
lines(pselm1.2, aucm1.2, type="l", col=2)
arrows(pselm1.2[aucv1.2!=0], (aucm1.2 - aucv1.2)[aucv1.2!=0], 
       pselm1.2[aucv1.2!=0], (aucm1.2 + aucv1.2)[aucv1.2!=0], length=0.05, 
       angle=90, code=3, col=2)
legend("bottomright", legend=c("regular", "fused"), lty=c(1, 1), col=c(1, 2))

tempA <- A1 
tempA[upper.tri(tempA)] <- 0
graph1 <- (graph_from_adjacency_matrix(tempA, "directed") + "Y") %>% 
  add_edges(c("X1", "Y", "X2", "Y", "X3", "Y"))
plot(graph1)
text(1, -1, labels=bquote(atop(rho== .(rho)~",", alpha==.(alpha))))
par(mfrow=c(1, 1))



# simulation 2
n <- 100
ntest <- 1000
p <- 9
rho <- 0.7
beta <- rep(c(1, 0, 0), 3)
Sigma <- matrix(rho, 3, 3)
diag(Sigma) <- 1

A2 <- as.matrix(bdiag(matrix(1, 3, 3), matrix(1, 3, 3), matrix(1, 3, 3)))
diag(A2) <- 0
rownames(A2) <- colnames(A2) <- paste("X", c(1:9), sep="")

graph.list2 <- list(lapply(seq_len(p), function(i) apply(A2, 1, function(a) {
  which(a!=0)})[, i]), lapply(seq_len(p), function(i) apply(A2, 1, function(a) {
    rep(1, sum(a))})[, i]))

alpha <- 0.5

set.seed(2018)
nreps <- 100
psel2.1 <- psel2.2 <- auc2.1 <- auc2.2 <- vector(mode="list", length=nreps)
for(rep in 1:nreps) {
  x <- matrix(0, n, p)
  x[, c(1:3)] <- rmvnorm(n, rep(0, 3), Sigma)
  x[, c(4:6)] <- rmvnorm(n, rep(0, 3), Sigma)
  x[, c(7:9)] <- rmvnorm(n, rep(0, 3), Sigma)
  y <- rbinom(n, 1, 1/(1 + exp(-x %*% beta)))

  xtest <- matrix(0, ntest, p)
  xtest[, c(1:3)] <- rmvnorm(n, rep(0, 3), Sigma)
  xtest[, c(4:6)] <- rmvnorm(n, rep(0, 3), Sigma)
  xtest[, c(7:9)] <- rmvnorm(n, rep(0, 3), Sigma)
  ytest <- rbinom(ntest, 1, 1/(1 + exp(-xtest %*% beta)))

  test2.1 <- glmnet(x, y, alpha=1, family="binomial", standardize=FALSE)
  pred2.1 <- predict(test2.1, xtest)
  psel2.1[[rep]] <- test2.1$df
  auc2.1[[rep]] <- apply(pred2.1, 2, function(pred) {
    pROC::roc(ytest, pred)$auc})
  
  seq.lambda <- test2.1$lambda*2
  
  psel2.2[[rep]] <- auc2.2[[rep]] <- numeric(length(seq.lambda))
  for(l in 1:length(seq.lambda)) {
    test2.2 <- fusedlasso(x, y, lambda1=seq.lambda[l]*alpha, 
                          lambda2=seq.lambda[l]*(1 - alpha), family="binomial", 
                          graph=graph.list2)
    pred2.2 <- as.numeric(1/(1 + exp(-cbind(1, xtest) %*% 
                                       rbind(test2.2$Intercept, test2.2$beta))))
    psel2.2[[rep]][l] <- sum(test2.2$beta!=0)
    auc2.2[[rep]][l] <- pROC::roc(ytest, pred2.2)$auc
  }
}

aucd2.1 <- cbind(Reduce(c, psel2.1), Reduce(c, auc2.1))
aucm2.1 <- sapply(sort(unique(aucd2.1[, 1])), function(psel) {
  mean(aucd2.1[aucd2.1[, 1]==psel, 2], na.rm=TRUE)})
aucv2.1 <- sapply(sort(unique(aucd2.1[, 1])), function(psel) {
  var(aucd2.1[aucd2.1[, 1]==psel, 2], na.rm=TRUE)})
pselm2.1 <- sort(unique(aucd2.1[, 1]))

aucd2.2 <- cbind(Reduce(c, psel2.2), Reduce(c, auc2.2))
aucm2.2 <- sapply(sort(unique(aucd2.2[, 1])), function(psel) {
  mean(aucd2.2[aucd2.2[, 1]==psel, 2], na.rm=TRUE)})
aucv2.2 <- sapply(sort(unique(aucd2.2[, 1])), function(psel) {
  var(aucd2.2[aucd2.2[, 1]==psel, 2], na.rm=TRUE)})
pselm2.2 <- sort(unique(aucd2.2[, 1]))

par(mfrow=c(1, 2))
xlim <- range(pselm2.1, pselm2.2)
ylim <- range(aucm2.1 - aucv2.1, aucm2.1 + aucv2.1, aucm2.2 - aucv2.2, 
              aucm2.2 + aucv2.2)
plot(pselm2.1, aucm2.1, type="l", xlim=xlim, ylim=ylim, 
     xlab="Number of selected variables", ylab="AUC")
arrows(pselm2.1[aucv2.1!=0], (aucm2.1 - aucv2.1)[aucv2.1!=0], 
       pselm2.1[aucv2.1!=0], (aucm2.1 + aucv2.1)[aucv2.1!=0], length=0.05, 
       angle=90, code=3)
lines(pselm2.2, aucm2.2, type="l", col=2)
arrows(pselm2.2[aucv2.2>1e-5], (aucm2.2 - aucv2.2)[aucv2.2>1e-5], 
       pselm2.2[aucv2.2>1e-5], (aucm2.2 + aucv2.2)[aucv2.2>1e-5], 
       length=0.05, 
       angle=90, code=3, col=2)
legend("bottomright", legend=c("regular", "fused"), lty=c(1, 1), col=c(1, 2))

tempA <- A2 
tempA[upper.tri(tempA)] <- 0
graph2 <- (graph_from_adjacency_matrix(A2, "undirected") + "Y") %>% 
  add_edges(c("X1", "Y", "X4", "Y", "X7", "Y"))
plot(graph2)
text(1, -1, labels=bquote(atop(rho== .(rho)~",", alpha==.(alpha))))
par(mfrow=c(1, 1))



# simulation 3
n <- 100
ntest <- 1000
p <- 9
rho <- 0.5
beta <- rep(0.5, 9)
Sigma <- matrix(rho, 3, 3)
diag(Sigma) <- 1

A3 <- as.matrix(bdiag(matrix(1, 3, 3), matrix(1, 3, 3), matrix(1, 3, 3)))
diag(A3) <- 0
rownames(A3) <- colnames(A3) <- paste("X", c(1:9), sep="")

graph.list3 <- list(lapply(seq_len(p), function(i) apply(A3, 1, function(a) {
  which(a!=0)})[, i]), lapply(seq_len(p), function(i) apply(A3, 1, function(a) {
    rep(1, sum(a))})[, i]))

alpha <- 0.3

set.seed(2018)
nreps <- 100
psel3.1 <- psel3.2 <- auc3.1 <- auc3.2 <- vector(mode="list", length=nreps)
for(rep in 1:nreps) {
  x <- matrix(0, n, p)
  x[, c(1:3)] <- rmvnorm(n, rep(0, 3), Sigma)
  x[, c(4:6)] <- rmvnorm(n, rep(0, 3), Sigma)
  x[, c(7:9)] <- rmvnorm(n, rep(0, 3), Sigma)
  y <- rbinom(n, 1, 1/(1 + exp(-x %*% beta)))
  
  xtest <- matrix(0, ntest, p)
  xtest[, c(1:3)] <- rmvnorm(n, rep(0, 3), Sigma)
  xtest[, c(4:6)] <- rmvnorm(n, rep(0, 3), Sigma)
  xtest[, c(7:9)] <- rmvnorm(n, rep(0, 3), Sigma)
  ytest <- rbinom(ntest, 1, 1/(1 + exp(-xtest %*% beta)))
  
  test3.1 <- glmnet(x, y, alpha=1, family="binomial", standardize=FALSE)
  pred3.1 <- predict(test3.1, xtest)
  psel3.1[[rep]] <- test3.1$df
  auc3.1[[rep]] <- apply(pred3.1, 2, function(pred) {
    pROC::roc(ytest, pred)$auc})
  
  seq.lambda <- test3.1$lambda*2
  
  psel3.2[[rep]] <- auc3.2[[rep]] <- numeric(length(seq.lambda))
  for(l in 1:length(seq.lambda)) {
    test3.2 <- fusedlasso(x, y, lambda1=seq.lambda[l]*alpha, 
                          lambda2=seq.lambda[l]*(1 - alpha), family="binomial", 
                          graph=graph.list3)
    pred3.2 <- as.numeric(1/(1 + exp(-cbind(1, xtest) %*% 
                                       rbind(test3.2$Intercept, test3.2$beta))))
    psel3.2[[rep]][l] <- sum(test3.2$beta!=0)
    auc3.2[[rep]][l] <- pROC::roc(ytest, pred3.2)$auc
  }
}

aucd3.1 <- cbind(Reduce(c, psel3.1), Reduce(c, auc3.1))
aucm3.1 <- sapply(sort(unique(aucd3.1[, 1])), function(psel) {
  mean(aucd3.1[aucd3.1[, 1]==psel, 2], na.rm=TRUE)})
aucv3.1 <- sapply(sort(unique(aucd3.1[, 1])), function(psel) {
  var(aucd3.1[aucd3.1[, 1]==psel, 2], na.rm=TRUE)})
pselm3.1 <- sort(unique(aucd3.1[, 1]))

aucd3.2 <- cbind(Reduce(c, psel3.2), Reduce(c, auc3.2))
aucm3.2 <- sapply(sort(unique(aucd3.2[, 1])), function(psel) {
  mean(aucd3.2[aucd3.2[, 1]==psel, 2], na.rm=TRUE)})
aucv3.2 <- sapply(sort(unique(aucd3.2[, 1])), function(psel) {
  var(aucd3.2[aucd3.2[, 1]==psel, 2], na.rm=TRUE)})
pselm3.2 <- sort(unique(aucd3.2[, 1]))

par(mfrow=c(1, 2))
xlim <- range(pselm3.1, pselm3.2)
ylim <- range(aucm3.1 - aucv3.1, aucm3.1 + aucv3.1, aucm3.2 - aucv3.2, 
              aucm3.2 + aucv3.2)
plot(pselm3.1, aucm3.1, type="l", xlim=xlim, ylim=ylim, 
     xlab="Number of selected variables", ylab="AUC")
arrows(pselm3.1[aucv3.1!=0], (aucm3.1 - aucv3.1)[aucv3.1!=0], 
       pselm3.1[aucv3.1!=0], (aucm3.1 + aucv3.1)[aucv3.1!=0], length=0.05, 
       angle=90, code=3)
lines(pselm3.2, aucm3.2, type="l", col=2)
arrows(pselm3.2[aucv3.2>1e-5], (aucm3.2 - aucv3.2)[aucv3.2>1e-5], 
       pselm3.2[aucv3.2>1e-5], (aucm3.2 + aucv3.2)[aucv3.2>1e-5], 
       length=0.05, 
       angle=90, code=3, col=2)
legend("bottomright", legend=c("regular", "fused"), lty=c(1, 1), col=c(1, 2))

tempA <- A3
tempA[upper.tri(tempA)] <- 0
graph3 <- (graph_from_adjacency_matrix(A3, "undirected") + "Y") %>% 
  add_edges(insert(rep("Y", 9), c(1:9), paste("X", c(1:9), sep="")))
plot(graph3)
text(1, -1, labels=bquote(atop(rho== .(rho)~",", alpha==.(alpha))))
par(mfrow=c(1, 1))



# simulation 4
n <- 100
ntest <- 1000
p <- 18
rho <- 0.7
beta <- rep(c(1, rep(0, 5)), 3)
Sigma <- matrix(rho, 6, 6)
diag(Sigma) <- 1

A4 <- as.matrix(bdiag(matrix(1, 6, 6), matrix(1, 6, 6), matrix(1, 6, 6)))
diag(A4) <- 0
rownames(A4) <- colnames(A4) <- paste("X", c(1:18), sep="")

graph.list4 <- list(lapply(seq_len(p), function(i) apply(A4, 1, function(a) {
  which(a!=0)})[, i]), lapply(seq_len(p), function(i) apply(A4, 1, function(a) {
    rep(1, sum(a))})[, i]))

alpha <- 0.7

set.seed(2018)
nreps <- 100
psel4.1 <- psel4.2 <- auc4.1 <- auc4.2 <- vector(mode="list", length=nreps)
for(rep in 1:nreps) {
  x <- matrix(0, n, p)
  x[, c(1:6)] <- rmvnorm(n, rep(0, 6), Sigma)
  x[, c(7:12)] <- rmvnorm(n, rep(0, 6), Sigma)
  x[, c(13:18)] <- rmvnorm(n, rep(0, 6), Sigma)
  y <- rbinom(n, 1, 1/(1 + exp(-x %*% beta)))
  
  xtest <- matrix(0, ntest, p)
  xtest[, c(1:6)] <- rmvnorm(n, rep(0, 6), Sigma)
  xtest[, c(7:12)] <- rmvnorm(n, rep(0, 6), Sigma)
  xtest[, c(13:18)] <- rmvnorm(n, rep(0, 6), Sigma)
  ytest <- rbinom(ntest, 1, 1/(1 + exp(-xtest %*% beta)))
  
  test4.1 <- glmnet(x, y, alpha=1, family="binomial", standardize=FALSE)
  pred4.1 <- predict(test4.1, xtest)
  psel4.1[[rep]] <- test4.1$df
  auc4.1[[rep]] <- apply(pred4.1, 2, function(pred) {
    pROC::roc(ytest, pred)$auc})
  
  seq.lambda <- test4.1$lambda*2
  
  psel4.2[[rep]] <- auc4.2[[rep]] <- numeric(length(seq.lambda))
  for(l in 1:length(seq.lambda)) {
    test4.2 <- fusedlasso(x, y, lambda1=seq.lambda[l]*alpha, 
                          lambda2=seq.lambda[l]*(1 - alpha), family="binomial", 
                          graph=graph.list4)
    pred4.2 <- as.numeric(1/(1 + exp(-cbind(1, xtest) %*% 
                                       rbind(test4.2$Intercept, test4.2$beta))))
    psel4.2[[rep]][l] <- sum(test4.2$beta!=0)
    auc4.2[[rep]][l] <- pROC::roc(ytest, pred4.2)$auc
  }
}

aucd4.1 <- cbind(Reduce(c, psel4.1), Reduce(c, auc4.1))
aucm4.1 <- sapply(sort(unique(aucd4.1[, 1])), function(psel) {
  mean(aucd4.1[aucd4.1[, 1]==psel, 2], na.rm=TRUE)})
aucv4.1 <- sapply(sort(unique(aucd4.1[, 1])), function(psel) {
  var(aucd4.1[aucd4.1[, 1]==psel, 2], na.rm=TRUE)})
pselm4.1 <- sort(unique(aucd4.1[, 1]))

aucd4.2 <- cbind(Reduce(c, psel4.2), Reduce(c, auc4.2))
aucm4.2 <- sapply(sort(unique(aucd4.2[, 1])), function(psel) {
  mean(aucd4.2[aucd4.2[, 1]==psel, 2], na.rm=TRUE)})
aucv4.2 <- sapply(sort(unique(aucd4.2[, 1])), function(psel) {
  var(aucd4.2[aucd4.2[, 1]==psel, 2], na.rm=TRUE)})
pselm4.2 <- sort(unique(aucd4.2[, 1]))

par(mfrow=c(1, 2))
xlim <- range(pselm4.1, pselm4.2)
ylim <- range(aucm4.1 - aucv4.1, aucm4.1 + aucv4.1, aucm4.2 - aucv4.2, 
              aucm4.2 + aucv4.2, na.rm=TRUE)
plot(pselm4.1, aucm4.1, type="l", xlim=xlim, ylim=ylim, 
     xlab="Number of selected variables", ylab="AUC")
arrows(pselm4.1[aucv4.1!=0], (aucm4.1 - aucv4.1)[aucv4.1!=0], 
       pselm4.1[aucv4.1!=0], (aucm4.1 + aucv4.1)[aucv4.1!=0], length=0.05, 
       angle=90, code=3)
lines(pselm4.2, aucm4.2, type="l", col=2)
arrows(pselm4.2[aucv4.2>1e-5], (aucm4.2 - aucv4.2)[aucv4.2>1e-5], 
       pselm4.2[aucv4.2>1e-5], (aucm4.2 + aucv4.2)[aucv4.2>1e-5], 
       length=0.05, 
       angle=90, code=3, col=2)
legend("bottomright", legend=c("regular", "fused"), lty=c(1, 1), col=c(1, 2))

tempA <- A4 
tempA[upper.tri(tempA)] <- 0
graph4 <- (graph_from_adjacency_matrix(A4, "undirected") + "Y") %>% 
  add_edges(c("X1", "Y", "X7", "Y", "X13", "Y"))
plot(graph4)
text(1, -1, labels=bquote(atop(rho== .(rho)~",", alpha==.(alpha))))
par(mfrow=c(1, 1))




# simulation 5
n <- 100
ntest <- 1000
p <- 30
rho <- 0.3

beta <- rep(c(0.5, rep(0, 9)), 3)
A5 <- matrix(0, p, p)
A5[c(2:10), 1] <- A5[c(12:20), 11] <- A5[c(22:30), 21] <- 1
A5[upper.tri(A5)] <- t(A5)[upper.tri(A5)]
rownames(A5) <- colnames(A5) <- paste("X", c(1:p), sep="")

graph.list5 <- list(apply(A5, 1, function(a) {which(a!=0)}), 
                    apply(A5, 1, function(a) {rep(1, sum(a))}))

alpha <- 0.5
seq.lambda <- exp(seq(log(0.30), log(0.001), length.out=50))

set.seed(2018)
nreps <- 100
psel5.1 <- psel5.2 <- auc5.1 <- auc5.2 <- vector(mode="list", length=nreps)

for(rep in 1:nreps) {
  x <- matrix(0, n, p)
  x[, -c(1, 11, 21)] <- rnorm((p - 3)*n)
  x[, 1] <- rowSums(rho*x[, c(2:10)]) + rnorm(n, 0, sqrt(1 - 9*rho^2))
  x[, 11] <- rowSums(rho*x[, c(12:20)]) + rnorm(n, 0, sqrt(1 - 9*rho^2))
  x[, 21] <- rowSums(rho*x[, c(22:30)]) + rnorm(n, 0, sqrt(1 - 9*rho^2))
  y <- rbinom(n, 1, 1/(1 + exp(-x %*% beta)))
  
  xtest <- matrix(0, ntest, p)
  xtest[, -c(1, 11, 21)] <- rnorm((p - 3)*ntest)
  xtest[, 1] <- rowSums(rho*xtest[, c(2:10)]) + 
    rnorm(ntest, 0, sqrt(1 - 3*rho^2))
  xtest[, 11] <- rowSums(rho*xtest[, c(12:20)]) + 
    rnorm(ntest, 0, sqrt(1 - 3*rho^2))
  xtest[, 21] <- rowSums(rho*xtest[, c(22:30)]) + 
    rnorm(ntest, 0, sqrt(1 - 3*rho^2))
  ytest <- rbinom(ntest, 1, 1/(1 + exp(-xtest %*% beta)))
  
  test5.1 <- glmnet(x, y, alpha=1, family="binomial", standardize=FALSE)
  pred5.1 <- predict(test5.1, xtest)
  psel5.1[[rep]] <- test5.1$df
  auc5.1[[rep]] <- apply(pred5.1, 2, function(pred) {
    pROC::roc(ytest, pred)$auc})
  
  psel5.2[[rep]] <- auc5.2[[rep]] <- numeric(length(seq.lambda))
  for(l in 1:length(seq.lambda)) {
    test5.2 <- fusedlasso(x, y, lambda1=seq.lambda[l]*alpha, 
                          lambda2=seq.lambda[l]*(1 - alpha), family="binomial", 
                          graph=graph.list5)
    pred5.2 <- as.numeric(1/(1 + exp(-cbind(1, xtest) %*% 
                                       rbind(test5.2$Intercept, test5.2$beta))))
    psel5.2[[rep]][l] <- sum(test5.2$beta!=0)
    auc5.2[[rep]][l] <- pROC::roc(ytest, pred5.2)$auc
  }
}

aucd5.1 <- cbind(Reduce(c, psel5.1), Reduce(c, auc5.1))
aucm5.1 <- sapply(sort(unique(aucd5.1[, 1])), function(psel) {
  mean(aucd5.1[aucd5.1[, 1]==psel, 2], na.rm=TRUE)})
aucv5.1 <- sapply(sort(unique(aucd5.1[, 1])), function(psel) {
  var(aucd5.1[aucd5.1[, 1]==psel, 2], na.rm=TRUE)})
pselm5.1 <- sort(unique(aucd5.1[, 1]))

aucd5.2 <- cbind(Reduce(c, psel5.2), Reduce(c, auc5.2))
aucm5.2 <- sapply(sort(unique(aucd5.2[, 1])), function(psel) {
  mean(aucd5.2[aucd5.2[, 1]==psel, 2], na.rm=TRUE)})
aucv5.2 <- sapply(sort(unique(aucd5.2[, 1])), function(psel) {
  var(aucd5.2[aucd5.2[, 1]==psel, 2], na.rm=TRUE)})
pselm5.2 <- sort(unique(aucd5.2[, 1]))

par(mfrow=c(1, 2))
xlim <- range(pselm5.1, pselm5.2)
ylim <- range(aucm5.1 - aucv5.1, aucm5.1 + aucv5.1, aucm5.2 - aucv5.2, 
              aucm5.2 + aucv5.2)
plot(pselm5.1, aucm5.1, type="l", xlim=xlim, ylim=ylim, 
     xlab="Number of selected variables", ylab="AUC")
arrows(pselm5.1[aucv5.1!=0], (aucm5.1 - aucv5.1)[aucv5.1!=0], 
       pselm5.1[aucv5.1!=0], (aucm5.1 + aucv5.1)[aucv5.1!=0], length=0.05, 
       angle=90, code=3)
lines(pselm5.2, aucm5.2, type="l", col=2)
arrows(pselm5.2[aucv5.2!=0], (aucm5.2 - aucv5.2)[aucv5.2!=0], 
       pselm5.2[aucv5.2!=0], (aucm5.2 + aucv5.2)[aucv5.2!=0], length=0.05, 
       angle=90, code=3, col=2)
legend("bottomright", legend=c("regular", "fused"), lty=c(1, 1), col=c(1, 2))

tempA <- A5 
tempA[upper.tri(tempA)] <- 0
graph5 <- (graph_from_adjacency_matrix(tempA, "directed") + "Y") %>% 
  add_edges(c("X1", "Y", "X11", "Y", "X21", "Y"))
plot(graph5)
text(1, -1, labels=bquote(atop(rho== .(rho)~",", alpha==.(alpha))))
par(mfrow=c(1, 1))

# # simulated data 1
# ntrain <- 100
# ngenes <- 1000
# ntrans <- 100
# p <- ngenes + ntrans
# G <- ntrans
# msizes <- ngenes/ntrans + 1
# A <- matrix(0, p, p)
# for(g in 1:G) {
#   A[((g - 1)*msizes + 2):(g*msizes), ((g - 1)*msizes + 1)] <- 1
# }
# beta <- rep(0, p)
# nactmod <- 4
# effsizes <- c(2, -2, 4, -4)
# for(g in 1:nactmod) {
#   beta[((g - 1)*msizes + 1)] <- effsizes[g]
#   beta[((g - 1)*msizes + 2):(g*msizes)] <- effsizes[g]/sqrt(ngenes/ntrans)
# }
# rho <- 0.5
# set.seed(123)
# x <- matrix(0, n, p)
# for(g in 1:G) {
#   x[, ((g - 1)*msizes + 1)] <- rnorm(n, 0, 1)
#   x[, ((g - 1)*msizes + 2):(g*msizes)] <- 
#     0.7*matrix(rep(x[, ((g - 1)*msizes + 1)], 10), nrow=n, byrow=FALSE) + 
#     rnorm(n*ngenes/ntrans, 0, sqrt(1 - rho^2))
# }
# y <- as.numeric(rbinom(n, 1, exp(-x %*% beta)/(1 + exp(-x %*% beta))))
# m <- rep(1, n)
# 
# # estimation
# lambda1 <- 0.01
# lambda2 <- 0.01
# gamma1 <- 0.01
# gamma2 <- 0.01
# lambdag <- exp(seq(-1, 1, length.out=G))
# gammag <- exp(seq(-1, 1, length.out=G))
# 
# graph.list <- list(apply(A, 1, function(a) {which(a!=0)}), 
#                    apply(A, 1, function(a) {rep(1, sum(a))}))
# ### fused lasso
# library(FusedLasso)
# library(glmnet)
# test1 <- fusedlasso(x, y, lambda1=0.01, lambda2=0.000001, family="binomial", 
#                     graph=graph.list)
# test2 <- glmnet(x, y, lambda1=0.01, lambda2=0.000001, family="binomial", 
#                     graph=graph.list)
# 
# 
# which(test1$beta!=0)
# test1$beta[which(test1$beta!=0)]

# seq.lambda1 <- exp(seq(-5, -1, length.out=50))
# seq.lambda2 <- exp(seq(-10, -5, length.out=50))
# seq.logll <- matrix(0, length(seq.lambda1), length(seq.lambda2))
# for(l1 in 1:length(seq.lambda1)) {
#   for(l2 in 1:length(seq.lambda2)) {
#     for(i in 1:nrow(x)) {
#       fit <- fusedlasso(x[-i, ], y[-i], lambda1=seq.lambda1[l1], 
#                         lambda2=seq.lambda2[l2], family="binomial", 
#                         graph=graph.list)
#       pred <- 1/(1 + exp(-sum(c(1, x[i, ])*
#                                 as.numeric(rbind(fit$Intercept, fit$beta)))))
#       seq.logll[l1, l2] <- seq.logll[l1, l2] + y[i]*log(pred) + 
#         (1 - y[i])*log(1 - pred)
#     }
#   }
# }
# lambda1 <- seq.lambda1[which(seq.logll==max(seq.logll), arr.ind=TRUE)[1]]
# lambda2 <- seq.lambda2[which(seq.logll==max(seq.logll), arr.ind=TRUE)[2]]
# heatmap(seq.logll, Rowv=NA, Colv=NA, labRow=NA, labCol=NA)
# points(lambda1, lambda2, pch=4)

# set.seed(123)
# test1 <- cv.glmnet(x, y, alpha=1, family="binomial", 
#                    standardize=FALSE)
# test2 <- fusedlasso(x, y, lambda1=lambda1, lambda2=lambda2, family="binomial", 
#                     graph=graph.list)
# test3 <- fusedlasso(x, y, lambda1=test1$lambda.min, lambda2=lambda2*10, 
#                     family="binomial", graph=graph.list)
# best <- as.matrix(cbind(coef(test1, "lambda.min"),
#                         rbind(test2$Intercept, test2$beta), 
#                         rbind(test3$Intercept, test3$beta)))
# colnames(best) <- c("normal", "fused1", "fused2")
# preds <- 1/(1 + exp(-cbind(1, xtest) %*% best))
# 
# roc(ytest, preds[, 1])$auc
# roc(ytest, preds[, 2])$auc
# roc(ytest, preds[, 3])$auc





### KEGGgraph
# download and import pathway
setwd("/Users/magnusmunch/Documents/OneDrive/PhD/netren/data")
library(KEGGgraph)
tmp <- "hsa05210.xml"
retrieveKGML("05210", organism="hsa", destfile=tmp, method="curl")
network <- parseKGML2Graph(tmp)

# convert to igraph object
library(igraph)
network2 <- igraph.from.graphNEL(network, name=TRUE, weight=TRUE, 
                                 unlist.attrs=TRUE)

# set symbols and name of the vertices
library("org.Hs.eg.db")
genenames <- unlist(mget(x=translateKEGGID2GeneID(nodes(network)),
                         envir=org.Hs.egGENENAME))
genesymbols <- unlist(mget(x=translateKEGGID2GeneID(nodes(network)),
                           envir=org.Hs.egSYMBOL))
network3 <- set_vertex_attr(network2, "name", value=genesymbols)
network3 <- set_vertex_attr(network3, "label", value=genenames)

V(network3)$name
V(network3)$label
plot(network3)