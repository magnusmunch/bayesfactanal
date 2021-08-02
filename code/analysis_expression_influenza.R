### libraries
library(Biobase)
library(FMradio)
library(glmnet)
library(bayesfactanal)

### functions
# estimation functions
# source("code/factreg.R")
# source("code/bayesfactreg.R")

### load data
load("data/eset.filt416_expression_influenza.Rdata")

# #Get official gene symbols: TAKES SOME TIME!
# library("hgu133plus2.db")
# x <- hgu133plus2SYMBOL 
# symbolall <- Lkeys(x)
# for (i in 1:length(symbolall)) {
#   symbolall[i] <- get(symbolall[i], env=hgu133plus2SYMBOL)
# }

idu <- which(pData(eset.filt)$accnr=="GSE29614")
idl <- which(pData(eset.filt)$accnr=="GSE29617")

xu <- scale(t(exprs(eset.filt))[idu, ])
xl <- scale(t(exprs(eset.filt))[idl, ])
x <- rbind(xu, xl)

yu <- scale(pData(eset.filt)$logmaxtiter[idu])[, 1]
yl <- scale(pData(eset.filt)$logmaxtiter[idl])[, 1]

xy <- cbind(rbind(xu, xl), 
            c(rep(NA, length(yu)), yl))

### regular penalized methods
set.seed(2021)
time.ridge <- proc.time()[3]
fit.ridge <- cv.glmnet(xl, yl, alpha=0, intercept=FALSE, grouped=FALSE)
time.ridge <- proc.time()[3] - time.ridge
time.lasso <- proc.time()[3]
fit.lasso <- cv.glmnet(xl, yl, alpha=1, intercept=FALSE, grouped=FALSE)
time.lasso <- proc.time()[3] - time.lasso

### labeled only
# find correlation matrix
set.seed(2021)
time.factanal1 <- proc.time()[3]
fit.covmat1 <- bayesfactanal::regcor(x=xl, Gamma=rep(1, ncol(x)), nfolds=5, 
                                     lower=0, upper=1)
covmat1 <- fit.covmat1$cor  
est.gamma1 <- fit.covmat1$gamma

# find factors
# dimGB1 <- dimGB(covmat1)
# dimVAR1 <- dimVAR(covmat1, max(dimGB1), graph=TRUE)
est.factors1 <- sum(svd(xl/sqrt(nrow(xl) - 1), nu=0, nv=0)$d > 1)
est.factors1 <- max(est.factors1, 2)
# methods are in agreement

# fit factor model with varimax rotation
fit.factanal1 <- factanal(factors=est.factors1, covmat=covmat1, 
                          rotation="varimax")
scores1 <- facScore(xl, fit.factanal1$loadings, 
                    diag(fit.factanal1$uniquenesses))
fit.factanal1 <- c(factanal=list(fit.factanal1),
                   glmnet=list(
                     cv.glmnet(as.matrix(scores1), yl, alpha=0, intercept=FALSE,
                               grouped=FALSE)))
time.factanal1 <- proc.time()[3] - time.factanal1

# check determinacy of factor scores
# df1 <- facSMC(covmat1, fit.factanal1$loadings); df1

### unlabeled + labeled
# find correlation matrix
set.seed(2021)
time.factanal2 <- proc.time()[3]
fit.covmat2 <- bayesfactanal::regcor(x, Gamma=rep(1, ncol(x)), nfolds=5, 
                                     lower=0, upper=1)
covmat2 <- fit.covmat2$cor  
est.gamma2 <- fit.covmat2$gamma

# find factors
# dimGB2 <- dimGB(covmat2)
# dimVAR2 <- dimVAR(covmat2, max(dimGB2), graph=TRUE)
est.factors2 <- sum(svd(x/sqrt(nrow(x) - 1), nu=0, nv=0)$d > 1)
est.factors2 <- max(est.factors2, 2)
# methods are in agreement

# fit factor model with varimax rotation
fit.factanal2 <- factanal(factors=est.factors2, covmat=covmat2, 
                          rotation="varimax")
scores2 <- facScore(xl, fit.factanal2$loadings, 
                    diag(fit.factanal2$uniquenesses))
fit.factanal2 <- c(factanal=list(fit.factanal2),
                   glmnet=list(
                     cv.glmnet(as.matrix(scores2), yl, alpha=0, intercept=FALSE,
                               grouped=FALSE)))
time.factanal2 <- proc.time()[3] - time.factanal2

# check determinacy of factor scores
# df2 <- facSMC(covmat2, fit.factanal2$loadings); df2

### penalized and Bayesian factor regressions
time.emfactanal <- proc.time()[3]
fit.emfactanal <- tryCatch({
  emfactanal(x=xy, factors=est.factors2, gamma=est.gamma2, 
             Gamma=rep(1, ncol(xy)), rotation="none", start=NULL,
             control=list(maxit=1000, nstart=10, lower=0.005,
                          epsilon=sqrt(.Machine$double.eps)))},
  error=function(e) {list()})
if(length(fit.emfactanal)==0) {
  time.emfactanal <- NA
} else {
  time.emfactanal <- proc.time()[3] - time.emfactanal
}

hyper <- list(kappa=rep(9, ncol(x) + 1), nu=rep(4, ncol(x) + 1))
hyper <- c(hyper, list(gamma=(hyper$kappa - 1)/(hyper$nu*est.factors2)))
time.vbfactanal <- proc.time()[3]
fit.vbfactanal <- bayesfactanal(x=xy, factors=est.factors2, 
                               hyper=hyper, method="vb", rescale=TRUE,
                               control=list(maxit=1000, nsamples=1000, 
                                            epsilon=sqrt(.Machine$double.eps),
                                            trace=FALSE))
time.vbfactanal <- proc.time()[3] - time.vbfactanal

### assessing prediction error
# create folds
set.seed(2021)
n <- length(yl)
nfolds <- n
foldid <- sample(rep(1:nfolds, times=round(c(rep(
  n %/% nfolds + as.numeric((n %% nfolds)!=0), times=n %% nfolds),
  rep(n %/% nfolds, times=nfolds - n %% nfolds)))))

# set hyperparameters
hyper <- list(kappa=rep(9, ncol(x) + 1), nu=rep(4, ncol(x) + 1))
hyper <- c(hyper, list(gamma=(hyper$kappa - 1)/(hyper$nu*est.factors2)))

# loop over folds
methods <- c("ridge", "lasso", "FMradio1", "FMradio2", "emfactanal", 
             "vbfactanal")
pred <- matrix(NA, nrow=n, ncol=length(methods), dimnames=list(NULL, methods))
for(k in 1:nfolds) {
  
  cat("fold", k, "of", nfolds, "\r")
  
  # create training and test data
  xtrain <- scale(xl[foldid!=k, , drop=FALSE])
  xtest <- xl[foldid==k, , drop=FALSE]
  
  scores1train <- as.matrix(scores1[foldid!=k, , drop=FALSE])
  scores1test <- as.matrix(scores1[foldid==k, , drop=FALSE])
  
  scores2train <- as.matrix(scores2[foldid!=k, , drop=FALSE])
  scores2test <- as.matrix(scores2[foldid==k, , drop=FALSE])
  
  ytrain <- scale(yl[foldid!=k])
  
  xytrain <- cbind(rbind(xu, xl), 
                   c(rep(NA, length(yu)), 
                     replace(yl, which(foldid==k), NA)))
  
  # fit models
  cv.ridge <- cv.glmnet(xtrain, ytrain, alpha=0, intercept=FALSE, 
                        grouped=FALSE)
  cv.lasso <- cv.glmnet(xtrain, ytrain, alpha=1, intercept=FALSE, 
                        grouped=FALSE)
  cv.FMradio1 <- cv.glmnet(scores1train, ytrain, alpha=0, intercept=FALSE,
                           grouped=FALSE)
  cv.FMradio2 <- cv.glmnet(scores2train, ytrain, alpha=0, intercept=FALSE,
                           grouped=FALSE)
  
  
  # cv.emfactanal <- emfactanal(x=na.omit(xytrain), factors=est.factors2, 
  #                             gamma=est.gamma2, 
  #                             Gamma=rep(1, ncol(xytrain)), 
  #                             rotation="none", start=NULL,
  #                             control=list(maxit=1000, nstart=10, lower=0.005, 
  #                                          epsilon=sqrt(.Machine$double.eps)))
  
  cv.emfactanal <- tryCatch({
    emfactanal(x=xytrain, factors=est.factors2,
               gamma=est.gamma2,
               Gamma=rep(1, ncol(xytrain)),
               rotation="none", start=NULL,
               control=list(maxit=1000, nstart=10, lower=0.005,
                            epsilon=sqrt(.Machine$double.eps)))},
    error=function(e) {list()})
  
  cv.vbfactanal <- bayesfactanal(x=xytrain, factors=est.factors2, 
                                 hyper=hyper, method="vb", rescale=TRUE,
                                 control=list(maxit=1000, nsamples=1000, 
                                              epsilon=sqrt(.Machine$double.eps),
                                              trace=FALSE))
  
  # create emfactanal predictions
  if(length(cv.emfactanal)==0) {
    pred.emfactanal <- rep(NA, nrow(xtest))
  } else {
    pred.emfactanal <- predict.emfactreg(xtest, fit=cv.emfactanal)$pred
  }
  
  # make predictions
  pred[foldid==k, ] <- cbind(ridge=predict(cv.ridge, xtest, s="lambda.min"), 
                             lasso=predict(cv.lasso, xtest, s="lambda.min"), 
                             FMradio1=predict(
                               cv.FMradio1, scores1test, s="lambda.min"),
                             FMradio2=predict(
                               cv.FMradio2, scores2test, s="lambda.min"),
                             emfactanal=pred.emfactanal,
                             vbfactanal=predict.vbfactreg(
                               xtest, fit=cv.vbfactanal)$pred)
}

pcor <- c(apply(pred, 2, cor, y=yl), null=0)
pmse <- c(colMeans((pred - yl)^2), null=mean(yl^2))
time <- c(ridge=time.ridge, lasso=time.lasso, FMradio1=time.factanal1, 
          FMradio2=time.factanal2, emfactanal=time.emfactanal, 
          vbfactanal=time.vbfactanal)

save(pcor, pmse, pred, time,
     file="results/analysis_expression_influenza_filt416_res1.Rdata")

