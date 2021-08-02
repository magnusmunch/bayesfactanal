### libraries
library(CoRF)
library(bayesfactanal)
library(glmnet)
library(pROC)
library(GRridge)
library(FMradio)
library(freeknotsplines)

### load data
data(LNM_Example)

### data preparation
# select genes
fid <- which(na.omit(CoDataTrain$pvalsVUmc <= 0.01))

# create codata
CDT <- CoDataTrain[fid, ]
grridge2gren <- function(partit){
  vec <- array(dim=length(unlist(partit)))
  for(i in 1:length(partit)){
    vec[partit[[i]]] <- i
  }
  return(vec)
}

# pv2 <- (CDT$pvalsVUmc <= spl2@optknot[1]) + (CDT$pvalsVUmc <= spl2@optknot[2]) +
#   (CDT$pvalsVUmc <= spl2@optknot[3]) + (CDT$pvalsVUmc <= spl2@optknot[4]) + 1

# pvpart <- grridge2gren(CreatePartition(CDT$pvalsVUmc, ngroup=3, 
#                                        decreasing=FALSE))

corkn <- knots(ecdf(CDT$Corrs))
corspl <- freelsgen(corkn, 1:length(corkn), degree=1, numknot=2, seed=2021,
                    stream=0)
corpart1 <- (CDT$Corrs <= corspl@optknot[1]) +
  (CDT$Corrs <= corspl@optknot[2]) + 1
# + (CDT$Corrs <= corspl@optknot[3]) + (CDT$Corrs <= spl1@optknot[4])

corpart2 <- grridge2gren(CreatePartition(CDT$Corrs, ngroup=3, decreasing=TRUE))
# cbind(corpart2, CDT$Corrs)[order(-CDT$Corrs), ]

# create data
xu <- scale(ValidationData[, fid])
xl <- scale(TrainData[, fid])
x <- rbind(xu, xl)

yu <- as.numeric(levels(RespValidation))[RespValidation]
yl <- RespTrain
y <- c(rep(NA, length(yu)), yl)

xy <- cbind(rbind(xu, xl), 
            c(rep(NA, length(yu)), yl))

### fit models
# penalized models
set.seed(2021)
time.ridge <- proc.time()[3]
fit.ridge <- cv.glmnet(xl, yl, family="binomial", alpha=0, intercept=FALSE, 
                       grouped=FALSE)
time.ridge <- proc.time()[3] - time.ridge
time.lasso <- proc.time()[3]
fit.lasso <- cv.glmnet(xl, yl, family="binomial", alpha=1, intercept=FALSE, 
                       grouped=FALSE)
time.lasso <- proc.time()[3] - time.lasso

### two step factor analysis with labeled only
# find correlation matrix
set.seed(2021)
time.FMradio1 <- proc.time()[3]
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
fit.FMradio1 <- factanal(factors=est.factors1, covmat=covmat1, 
                          rotation="varimax")
scoresl1 <- facScore(xl, fit.FMradio1$loadings, 
                    diag(fit.FMradio1$uniquenesses))
scoresu1 <- facScore(xu, fit.FMradio1$loadings, 
                     diag(fit.FMradio1$uniquenesses))
fit.FMradio1 <- c(factanal=list(fit.FMradio1),
                   glmnet=list(
                     cv.glmnet(as.matrix(scoresl1), yl, family="binomial", 
                               alpha=0, intercept=TRUE)))
time.FMradio1 <- proc.time()[3] - time.FMradio1
# check determinacy of factor scores
# df1 <- facSMC(covmat1, fit.factanal1$factanal$loadings); df1

### two step factor analysis with unlabeled + labeled
# find correlation matrix
set.seed(2021)
time.FMradio2 <- proc.time()[3]
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
fit.FMradio2 <- factanal(factors=est.factors2, covmat=covmat2, 
                          rotation="varimax")
scoresl2 <- facScore(xl, fit.FMradio2$loadings, 
                     diag(fit.FMradio2$uniquenesses))
scoresu2 <- facScore(xu, fit.FMradio2$loadings, 
                     diag(fit.FMradio2$uniquenesses))
fit.FMradio2 <- c(factanal=list(fit.FMradio2),
                   glmnet=list(
                     cv.glmnet(as.matrix(scoresl2), yl, family="binomial", 
                               alpha=0, intercept=TRUE)))
time.FMradio2 <- proc.time()[3] - time.FMradio2

# check determinacy of factor scores
# df2 <- facSMC(covmat2, fit.factanal2$factanal$loadings); df2

### VB models
set.seed(2021)
hyper1 <- list(kappa=rep(9, ncol(x)), nu=rep(4, ncol(x)))
hyper1 <- c(hyper1, list(gamma=list(
  # x=(hyper1$kappa - 1)/(hyper1$nu*est.factors2),
  x=rep(1/est.factors2, p),
  y=Inf)))
hyper2 <- c(hyper1, 
            list(gammagroups=list(x=c(1, 1, 1), y=1)),
            list(gammatotal=list(
              x=1/est.factors2,
              y=Inf)))
                
time.vbfactanal <- proc.time()[3]
fit.vbfactanal <- logbayesfactanal(x, as.matrix(y), factors=est.factors2, 
                                   hyper=hyper1, groups=NULL, method="vb", 
                                   start=NULL, rescale=TRUE, 
                                   control=list(maxit=1000, epsilon=1e-4,
                                                trace=TRUE))
time.vbfactanal <- proc.time()[3] - time.vbfactanal

time.ebfactanal1 <- proc.time()[3]
fit.ebfactanal1 <- logbayesfactanal(x, as.matrix(y), factors=est.factors2, 
                                   hyper=hyper2, groups=list(x=corpart1, y=0), 
                                   method="eb", start=fit.vbfactanal$param, 
                                   rescale=TRUE, 
                                   constrained=list(x=TRUE, y=FALSE),
                                   control=list(maxit=1000, epsilon=1e-4,
                                                trace=TRUE))
time.ebfactanal1 <- proc.time()[3] - time.ebfactanal1

time.ebfactanal2 <- proc.time()[3]
fit.ebfactanal2 <- logbayesfactanal(x, as.matrix(y), factors=est.factors2, 
                                   hyper=hyper2, 
                                   groups=list(x=as.numeric(corpart2), y=0), 
                                   method="eb", start=fit.vbfactanal$param, 
                                   rescale=TRUE, 
                                   constrained=list(x=TRUE, y=FALSE),
                                   control=list(maxit=1000, epsilon=1e-4,
                                                trace=TRUE))
time.ebfactanal2 <- proc.time()[3] - time.ebfactanal2

# ### EB models
# hyper1 <- list(kappa=rep(9, ncol(x)), nu=rep(4, ncol(x)), 
#                gamma=list(x=rep(1/est.factors2, ncol(x)), y=Inf),
#                gammagroups=list(x=rep(1, 3), y=1), 
#                gammatotal=list(x=1/est.factors2, y=Inf))
# hyper2 <- list(kappa=rep(9, ncol(x)), nu=rep(4, ncol(x)), 
#                gamma=list(x=rep(1/est.factors2, ncol(x)), y=1/est.factors2),
#                gammagroups=list(x=rep(1, 3), y=1), 
#                gammatotal=list(x=1/est.factors2, y=1/est.factors2))
# fit.ebfactanal1 <- logbayesfactanal(x, as.matrix(y), factors=est.factors2, 
#                                     hyper=hyper1, 
#                                     groups=list(x=part$corr, y=0), method="eb", 
#                                     start=NULL, rescale=TRUE, 
#                                     constrained=list(x=TRUE, y=FALSE), 
#                                     control=list(maxit=1000, epsilon=1e-4,
#                                                  trace=TRUE))
# fit.ebfactanal2 <- logbayesfactanal(x, as.matrix(y), factors=est.factors2, 
#                                     hyper=hyper2, 
#                                     groups=list(x=part$corr, y=0), method="eb", 
#                                     start=NULL, rescale=TRUE, 
#                                     constrained=list(x=TRUE, y=FALSE), 
#                                     control=list(maxit=1000, epsilon=1e-4,
#                                                  trace=TRUE))
# fit.ebfactanal3 <- logbayesfactanal(x, as.matrix(y), factors=est.factors2, 
#                                     hyper=hyper1, 
#                                     groups=list(x=part$pv, y=0), method="eb", 
#                                     start=NULL, rescale=TRUE, 
#                                     constrained=list(x=TRUE, y=FALSE), 
#                                     control=list(maxit=1000, epsilon=1e-4,
#                                                  trace=TRUE))
# fit.ebfactanal4 <- logbayesfactanal(x, as.matrix(y), factors=est.factors2, 
#                                     hyper=hyper2, 
#                                     groups=list(x=part$pv, y=0), method="eb", 
#                                     start=NULL, rescale=TRUE, 
#                                     constrained=list(x=TRUE, y=FALSE), 
#                                     control=list(maxit=1000, epsilon=1e-4,
#                                                  trace=TRUE))

# save(fit.ridge, fit.lasso, fit.vbfactanal1, fit.vbfactanal2, fit.vbfactanal1, 
#      fit.vbfactanal2, fit.vbfactanal3, fit.vbfactanal4, fit.ebfactanal1, 
#      fit.ebfactanal2, fit.ebfactanal3, fit.ebfactanal4,
#      file="results/rnaseq_oral_cancer_metastasis_res1.Rdata")

# load(file="results/rnaseq_oral_cancer_metastasis_res1.Rdata")

pred <- cbind(ridge=predict(fit.ridge, xu, type="response", 
                            s="lambda.min")[, 1],
              lasso=predict(fit.lasso, xu, type="response", 
                            s="lambda.min")[, 1],
              FMradio1=predict(fit.FMradio1$glmnet, as.matrix(scoresu1), 
                                s="lambda.min", type="response")[, 1],
              FMradio2=predict(fit.FMradio2$glmnet, as.matrix(scoresu2), 
                               s="lambda.min", type="response")[, 1],
              vbfactanal=predict.logvbfactreg(xu, fit=fit.vbfactanal)$pred[, 1],
              ebfactanal1=predict.logvbfactreg(xu, fit=fit.ebfactanal1)$pred[, 1],
              ebfactanal2=predict.logvbfactreg(xu, fit=fit.ebfactanal2)$pred[, 1])

auc <- c(apply(pred, 2, function(pr) {pROC::auc(yu, pr)}))
briers <- 1 - colMeans((pred - yu)^2)/mean((mean(yu) - yu)^2)
time <- c(ridge=time.ridge, lasso=time.lasso, FMradio1=time.FMradio1, 
          FMradio2=time.FMradio2, vbfactanal=time.vbfactanal, 
          ebfactanal1=time.ebfactanal1, ebfactanal2=time.ebfactanal2)

load(file="results/rnaseq_oral_cancer_metastasis_res4.Rdata")
# save(fit.ridge, fit.lasso, fit.FMradio1, fit.FMradio2, fit.vbfactanal,
#      fit.ebfactanal1, fit.ebfactanal2, pred, auc, briers, time,
#      file="results/rnaseq_oral_cancer_metastasis_res4.Rdata")

