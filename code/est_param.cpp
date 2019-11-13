// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

/* This function estimates new variational parameters*/
// [[Rcpp::export]]
List est_param(arma::mat x, arma::vec kappa, arma::vec m, arma::mat L, int n, 
               int p, int G, arma::vec sizes, arma::vec ciold, arma::vec etaold, 
               arma::mat Zold, double lambda1, double lambda2, double gamma1, 
               double gamma2, arma::vec lambdag, arma::vec gammag) {
  
  // initialize c, eta, zeta, sigma and mu
  arma::vec ci(n);
  arma::vec eta(p);
  arma::mat zeta(p,p);
  arma::mat sigma(p,p);
  arma::vec mu(p);
  
  arma::vec ominv = 0.5*ciold/(m%tanh(ciold/2));
  arma::mat Ainv(p,p);
  arma::mat xAinv(n,p);
  arma::mat xAinvxt(n,n);
  arma::mat xsigma(n,p);
  arma::vec lambdagvec(p);
  arma::vec gammagvec(p);

  arma::vec starts = cumsum(sizes);
  starts.insert_rows(0,1,true);
  starts.shed_row(G);
  for(int g=0; g<G; g++) {
    arma::mat Ag(sizes(g),sizes(g));
    Ag = gamma2*gammag(g)*L.submat(starts(g),starts(g),starts(g)+sizes(g)-1,
                       starts(g)+sizes(g)-1);
    Ag = Ag + gamma1*sqrt(gamma2)*gammag(g)*Zold.submat(starts(g),starts(g),
                          starts(g)+sizes(g)-1,starts(g)+sizes(g)-1);
    Ag.diag() += lambda2*lambdag(g) + lambda2*sqrt(lambda1)*lambdag(g)/
      sqrt(etaold.subvec(starts(g),starts(g)+sizes(g)-1));
    arma::mat Aginv = Ag.i();
    
    Ainv.submat(starts(g),starts(g),starts(g)+sizes(g)-1,
                starts(g)+sizes(g)-1) = Aginv;
    xAinv.submat(0,starts(g),n-1,starts(g)+sizes(g)-1) = 
      x.cols(starts(g),starts(g)+sizes(g)-1)*Aginv;
    
    arma::vec clambdag(sizes(g));
    arma::vec cgammag(sizes(g));
    lambdagvec.subvec(starts(g),starts(g)+sizes(g)-1) = 
      clambdag.fill(arma::as_scalar(lambdag(g)));
    gammagvec.subvec(starts(g),starts(g)+sizes(g)-1) = 
      cgammag.fill(arma::as_scalar(gammag(g)));
  }
  
  xAinvxt = xAinv*x.t();
  xAinvxt.diag() += ominv;
  xAinvxt = xAinvxt.i();
  sigma = Ainv - xAinv.t()*xAinvxt*xAinv;
  xsigma = x*sigma;
  mu = xsigma.t()*kappa;
  ci = sqrt(sum(xsigma % x,1) + arma::square(x*mu));
  eta = lambda2*lambdagvec%(diagvec(sigma) + arma::square(mu));
  
  for(arma::uword i=0; i<zeta.n_rows; ++i) {
    for(arma::uword j=0; j<i; j++) {
      if(Zold(i,j)!=0.0) {
        zeta(i,j) = gamma2*(pow(mu(i),2.0) + pow(mu(j),2.0) - 2*mu(i)*mu(j) +
          sigma.diag()(i) + sigma.diag()(j) - 2*sigma(i,j));
      }
    }
  }
  
  zeta = zeta + zeta.t();
  zeta = zeta.each_col()%gammagvec;
  
  return List::create(Named("ci") = ci,
                      Named("eta") = eta,
                      Named("sigma") = sigma,
                      Named("mu") = mu,
                      Named("zeta") = zeta);

}


