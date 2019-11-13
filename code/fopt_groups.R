# objective function for penalty multiplier estimation with constraint:
# gammag/lambdag=u
fopt_groups1 <- function(par, sizes, sum1) {
  
  s <- par[1]
  loglambdag <- par[-1]
  
  magn <- sum((0.5*exp(loglambdag)*sum1 + 
                 (s - 0.5)*sizes)^2) + sum(sizes*loglambdag)^2
  return(magn)
  
}

# under assumption of proportional multipliers
fopt_groups2 <- function(par, sizes, sum1, u) {
  
  s1 <- par[1]
  s2 <- par[2]
  loglambdag <- par[-c(1, 2)]
  
  magn <- sqrt(sum((0.5*exp(loglambdag)*sum1 + 
                      (s1 - 0.5)*sizes[[1]] + s2*sizes[[2]])^2) + 
                 sum(sizes[[1]]*loglambdag)^2 + 
                 (sum(sizes[[2]])*log(u) + sum(sizes[[2]]*loglambdag))^2)
  return(magn)
  
}





