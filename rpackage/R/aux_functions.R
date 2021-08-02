# expit function
.expit <- function(x) {exp(x)/(1 + exp(x))}

# log likelihood function
.llcor <- function(gamma, dval, dvec, n, nk, p, ldetGamma) {
  ll <- ldetGamma + max(0, p - n + nk)*log(gamma) + 
    min(p, n - nk)*log((1 - gamma)/(n - nk - 1)) +
    sum(log(dval^2 + (n - nk - 1)*gamma/(1 - gamma))) +
    (n - nk - 1)/((1 - gamma)*(nk - 1))*
    sum(dvec/(c(dval^2, rep(0, max(0, p - n + nk))) + 
                (n - nk - 1)*gamma/(1 - gamma)))
  return(ll)
}

# cross validated log likelihood function
.cvllcor <- function(gamma, dval, dvec, n, nk, p, ldetGamma, x, Gamma, foldid) {
  nfolds <- length(nk)
  
  # check whether we can use fast calculation
  if(is.null(ldetGamma)) {
    cvll <- 0
    for(k in 1:length(nk)) {
      R <- (1 - gamma)*cor(x[foldid!=k, , drop=FALSE])
      diag(R) <- diag(R) + gamma*Gamma
      cvll <- cvll + 
        (determinant(R)$modulus + 
           sum(diag(cor(x[foldid==k, , drop=FALSE]) %*% solve(R))))*nk
    }
    cvll <- cvll/nfolds
  } else {
    cvll <- sum(sapply(1:nfolds, function(k) {
      .llcor(gamma, dval[[k]], dvec[[k]], n, nk[k], p, ldetGamma)})*nk)/nfolds  
  }
  return(cvll)
}

# calculate constant part of ELBO
.const <- function(hyper, n, d, m) {
  
  # auxiliary variables
  p <- length(hyper$kappa)
  
  # constant ELBO part
  const <- -(n*p - sum(m))*log(2*pi)/2 + (n*d + n*p + d*p + sum(m))/2 +
    sum(lgamma(n/2
               # + d/2
               + hyper$kappa)) - sum(lgamma((hyper$kappa))) -
    d*sum(digamma(n/2
                  # + d/2
                  + hyper$kappa))/2 + 
    d*sum(log(hyper$gamma[is.finite(hyper$gamma)]))/2 +
    sum(hyper$kappa*(1 + log(hyper$nu)))
  return(const)
  
}

# calculate constant part of ELBO in logistic regression
.logconst <- function(hyper, d, linconst) {
  
  # auxiliary variables
  py <- length(hyper$gamma$y)
  
  # constant ELBO part
  const <- linconst + py*(d + 1)/2 + 
    d*sum(log(hyper$gamma$y[is.finite(hyper$gamma$y)]))/2
  return(const)
  
}
# 
# param <- linparam 
# hyper <- linhyper
# const <- 0
# m <- mx
# calculate ELBO
.elbo <- function(param, hyper, const, m) {
  
  # auxiliary variables
  p <- length(param$Omega)
  n <- nrow(param$Phi)
  d <- ncol(param$Phi)
  tau <- (n/2 + d/2 + hyper$kappa)/(2*param$zeta)
  
  # parts that change (excluding Xi)
  elbo <- sum(m*log(param$chi))/2 - sum(m*tau*param$chi) - 
    sum(((n - d)/2 + hyper$kappa)*log(param$zeta)) - sum(hyper$nu*tau)/2 +
    sum(sapply(1:p, function(j) {
      as.numeric(determinant(param$Omega[[j]])$modulus)/2 -
        tau[j]*sum(diag(param$Omega[[j]]))/hyper$gamma[j] - 
        tau[j]*sum(param$Phi*(param$Phi %*% param$Omega[[j]]))})) + 
    2*sum((param$Phi %*% param$Mu)*t(t(param$Upsilon)*tau))  -
    sum(colSums(param$Mu^2)*tau/hyper$gamma) -
    sum(colSums((param$Phi %*% param$Mu)^2)*tau) - 
    sum(tau*colSums(param$Upsilon^2)) - sum(param$Phi^2)/2
  
  # including Xi
  if(is.list(param$Xi)) {
    elbo <- elbo - sum(sapply(param$Xi, function(xi) {
      sum(sapply(1:p, function(j) {tau[j]*sum(xi*param$Omega[[j]])})) -
        sum((t(param$Mu) %*% xi)*(t(param$Mu)*tau)) +
        as.numeric(determinant(xi)$modulus)/2 - sum(diag(xi))/2}))
  } else {
    elbo <- elbo -
      sum(sapply(1:p, function(j) {n*tau[j]*sum(param$Xi*param$Omega[[j]])})) - 
      n*sum((t(param$Mu) %*% param$Xi)*(t(param$Mu)*tau)) + 
      n*as.numeric(determinant(param$Xi)$modulus)/2 - n*sum(diag(param$Xi))/2  
  }
  
  # add constants and return
  elbo <- elbo + const
  return(elbo)
  
}

# calculate ELBO
.logelbo <- function(param, hyper, const, linelbo, y) {
  
  # auxiliary variables
  py <- length(param$Omega$y)
  n <- nrow(param$Phi)
  d <- ncol(param$Phi)
  mid <- lapply(1:py, function(jy) {which(is.na(y[, jy]))})
  m <- sapply(mid, length)
  
  # logistic parts that change
  elbo <- ifelse(any(m!=0), -sum(sapply(which(m!=0), function(jy) {
    ups <- param$Upsilon$y[mid[[jy]], jy]
    sum((1 - ups)*log(1 - ups) + ups*log(ups))})),
    0) +
    sum((param$Upsilon$y - 0.5)*
          (t(t(param$Phi %*% param$Mu$y[-1, , drop=FALSE]) + param$Mu$y[1, ]) - 
             param$delta)) - 
    sum(log(1 + exp(param$delta))) +
    sum(tanh(param$delta/2)*param$delta^2/(4*param$delta)) -
    sum(tanh(param$delta/2)*param$delta^2/(4*param$delta)*
          sapply(1:n, function(i) {
            aux <- rbind(c(1, param$Phi[i, ]), cbind(
              param$Phi[i, ], as.matrix(param$Phi[i, ]) %*% t(param$Phi[i, ]) + 
                param$Xi[[i]]))
            colSums(param$Mu$y*(aux %*% param$Mu$y)) -
              sapply(1:py, function(jy) {
                sum(param$Omega$y[[jy]]*aux)})})) +
    sum(sapply(param$Omega$y, function(om) {
      as.numeric(determinant(om)$modulus)}))/2 -
    sum(sapply(1:py, function(jy) {
      0.5*sum(diag(param$Omega$y[[jy]][-1, -1]))/hyper$gamma$y[jy]})) - 
    sum(colSums(param$Mu$y[-1, , drop=FALSE]^2)/hyper$gamma$y)/2
  
  # add constants, linear part, and return
  elbo <- elbo + const + linelbo
  return(elbo)
  
}

# calculate regression coefficient from EM
.embeta <- function(param) {
  beta <- param$B[, ncol(param$B)]
  B <- param$B[, -ncol(param$B)]
  psi <- param$psi[-length(param$psi)]
  Btilde <- t(t(B)/psi)
  mat <- Btilde %*% t(B)
  diag(mat) <- diag(mat) + 1
  beta <- colSums(beta*(solve(mat) %*% Btilde))
  return(beta)
}

# calculate regression coefficient from MCMC
.mcmcbeta <- function(param, burnin=0) {
  d <- dim(param$B)[1]
  p <- dim(param$B)[2]
  nsamples <- dim(param$B)[3]
  beta <- matrix(nrow=p - 1, ncol=nsamples - burnin)
  for(k in 1:(nsamples - burnin)) {
    mat <- param$B[, -p, burnin + k] %*% 
      (t(param$B[, -p, burnin + k])/param$psi[-p, burnin + k])
    diag(mat) <- diag(mat) + 1
    beta[, k] <- colSums(solve(mat) %*% t(t(param$B[, -p, burnin + k])/
                                            param$psi[-p, burnin + k])*
                           param$B[, p, burnin + k]) 
  }
  return(rowMeans(beta))
}

# calculate regression coefficient from VB
.vbbeta <- function(param, hyper) {
  
  # auxiliary variables
  n <- nrow(param$Phi)
  p <- length(param$zeta)
  d <- nrow(param$Mu)
  
  # pre steps
  mat2 <- t(param$Mu[, -p])*(n/2 + hyper$kappa[-p] - 1)/param$zeta[-p]
  mat1 <- param$Mu[, -p] %*% mat2
  diag(mat1) <- diag(mat1) + 1
  mat1 <- -mat2 %*% solve(mat1) %*% t(mat2)
  diag(mat1) <- diag(mat1) + (n/2 + hyper$kappa[-p] - 1)/param$zeta[-p]
  mat2 <- mat1 %*% t(param$Mu[, -p])
  mat3 <- param$Mu[, -p] %*% mat2
  m <- param$Mu[, p]
  vec1 <- colSums(param$Mu[, -p]*m)
  
  vec2 <- vec3 <- numeric(p - 1)
  for(j in 1:(p - 1)) {
    vec3 <- vec3 + mat1[j, j]*sum(mat1[, j]*vec1)*param$zeta[-p][j]^2/
      ((n/2 + hyper$kappa[-p][j] - 1)^2*(n/2 + hyper$kappa[-p][j] - 2))*
      mat1[, j]
    for(h in 1:d) {
      for(l in 1:d) {
        vec2 <- vec2 + 
          ((sum((mat2[j, l]*mat2[j, ] + mat1[j, j]*mat3[l, ])*m) - 
              mat1[j, j]*m[l])*mat2[, h] +
             (sum((mat2[j, h]*mat2[j, ] + mat1[j, j]*mat3[h, ])*m) - 
                mat1[j, j]*m[l])*mat2[, l] +
             (sum((2*mat3[h, l]*mat2[j, ] - 2*(l==h)*mat2[j, ] + 
                     mat2[j, h]*mat3[l, ])*m) - mat2[j, l]*m[h] - 
                mat2[j, h]*m[l])*mat1[, j])*param$Omega[[j]][h, l]
      }
    }
  }
  
  # calculate and return beta
  naive <- as.numeric(colSums(t(mat2)*m))
  correction <- as.numeric(vec2/2 + vec3)
  beta <- naive + correction
  return(beta)
}