p <- 20
n <- 10
nedges <- 20
A1 <- matrix(0, p/2, p/2)
A1[upper.tri(A1)] <- sample(rep(c(1, 0), times=c(nedges/2, (p/2)*(p/2-1)/2 - 
                                                   nedges/2)))
A1[lower.tri(A1)] <- t(A1)[lower.tri(A1)]
A2 <- matrix(0, p/2, p/2)
A2[upper.tri(A2)] <- sample(rep(c(1, 0), times=c(nedges/2, (p/2)*(p/2-1)/2 - 
                                                   nedges/2)))
A2[lower.tri(A2)] <- t(A2)[lower.tri(A2)]
A <- rbind(cbind(A1, matrix(0, nrow=p/2, ncol=p/2)),
           cbind(matrix(0, nrow=p/2, ncol=p/2), A2))

chi <- rnorm(sum(A[upper.tri(A)]), 1)
Achi <- A
Achi[upper.tri(Achi)][Achi[upper.tri(Achi)]==1] <- chi
Achi[lower.tri(Achi)] <- t(Achi)[lower.tri(Achi)]
D <- diag(rowSums(A))
L <- D - A
Dchi <- diag(rowSums(Achi))
Lchi <- Dchi - Achi

gammag <- c(0.5, 1, 2)
Gammag <- matrix(0, nrow=p, ncol=p)
Gammag[1:(p/2), 1:(p/2)] <- gammag[1]
Gammag[(p/2 + 1):p, 1:(p/2)] <- Gammag[1:(p/2), (p/2 + 1):p] <- gammag[2]
Gammag[(p/2 + 1):p, (p/2 + 1):p] <- gammag[3]

beta <- rnorm(p)
(t(beta) %*% diag(diag(Gammag)) %*% L %*% beta)

bsum <- 0
for(i in 2:nrow(A)) {
  for(j in 1:(i - 1)) {
    if(A[i, j]==1) {bsum <- bsum + Gammag[i, j]*(beta[i] - beta[j])^2}
  }
}
bsum

diag(Gammag)[i]

# calculates the auxiliary variables (not tested)
.aux.var <- function(double aold, arma::vec bold, arma::vec y, arma::mat x, 
             arma::rowvec ytx) {
  
  H <- 
}