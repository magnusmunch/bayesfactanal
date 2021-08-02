#!/usr/bin/env Rscript

### libraries
library(TCGAbiolinks)
library(edgeR)
library(miRBaseConverter)
library(glmnet)
library(scout)  
library(netren)
library(beam)
library(pROC)

### parallelisation
parallel <- TRUE

### data
# study data
load("data/forMagnusN88.Rdata")

# TCGA data
load("data/TCGA-COAD.RData") #loads TCGA miRNA data

# micro RNAs
miRNA_TCGA <- X$miRNA
miRNA_study <- t(mirnormcen_resp)
colnames(miRNA_study) <- sapply(strsplit(colnames(miRNA_study), " "), "[[", 1)

# retrieve isomer miRNAs from TCGA
query_iso <- GDCquery(project="TCGA-COAD", 
                      data.category="Transcriptome Profiling", 
                      data.type="Isoform Expression Quantification",
                      sample.type="Primary Tumor")
GDCdownload(query_iso, directory="data")
tcga_iso <- GDCprepare(query=query_iso, summarizedExperiment=FALSE, 
                       directory="data")
tcga_iso$mirbase_accession <- sapply(strsplit(tcga_iso$miRNA_region, ","), tail, 
                                     n=1)
tcga_iso$miRNA_ID_full <- miRNA_AccessionToName(tcga_iso$mirbase_accession, 
                                                targetVersion="v22")$TargetName

# sum all precursors and versions that lead to the same mature miRNA
tcga_rpm <- aggregate(read_count ~ miRNA_ID_full + barcode, data=tcga_iso, sum)
tcga_mirna <- reshape(tcga_rpm, idvar="barcode", timevar="miRNA_ID_full", 
                      direction="wide")
colnames(tcga_mirna) <- c("barcode", unique(tcga_rpm$miRNA_ID_full))
subjectids <- sapply(strsplit(tcga_mirna$barcode, "-"), function(s) {
  paste0(s[1], "-", s[2], "-", s[3])})
tcga_mirna <- as.matrix(aggregate(tcga_mirna[, -1], list(subjectids), 
                                  mean)[, -1])
rownames(tcga_mirna) <- subjectids[!duplicated(subjectids)]
tcga_mirna[is.na(tcga_mirna)] <- 0
miRNA_TCGA <- tcga_mirna[, -which(colSums(tcga_mirna)==0)]

# normalise TCGA data
dge <- calcNormFactors(t(miRNA_TCGA), method="TMM")
tcga.prep <- scale(sqrt(miRNA_TCGA*dge))

# normalise study data
study.prep <- scale(miRNA_study)

# target vector
benefit <- as.numeric(resp) - 1

# unpenalized covariates
unpen <- model.matrix(~ adjth + thscheme + age + pcrcdiff, data=datfr)[, -1]

# save data
save(tcga.prep, study.prep, benefit, unpen,
     file="data/micrornaseq_colorectal_cancer_dat1.Rdata")

# load data
load(file="data/micrornaseq_colorectal_cancer_dat1.Rdata")

### analysis 1
# select common miRNAs
y <- ytest <- as.matrix(benefit)
x <- study.prep
u <- unpen
n <- nrow(x)
px <- ncol(x)
py <- ncol(y)
mid <- sort(c(sample(which(y==0), floor(sum(y==0)/2)),
              sample(which(y==1), floor(sum(y==1)/2))))
y[mid] <- NA


hyper <- list(kappa=rep(3, px), nu=rep(1, px), 
              gamma=list(x=rep(10, px), y=rep(Inf, py)))

est.factors <- max(sum(eigen(cor(x), only.values=TRUE)$values > 1), 2)

fit.vbfactanal1 <- logvbfactanal(x, y, factors=est.factors, hyper=hyper, 
                                 start=NULL, rescale=FALSE, 
                                 control=list(
                                   maxit=100,
                                   epsilon=1e-4))
fit.vbfactanal2 <- logvbfactanal(x, y, factors=20, hyper=hyper, 
                                 start=NULL, rescale=FALSE, 
                                 control=list(
                                   maxit=100,
                                   epsilon=1e-4))
fit.glmnet1 <- cv.glmnet(x[-mid, ], y[-mid, 1], family="binomial", alpha=0)
fit.glmnet2 <- cv.glmnet(x[-mid, ], y[-mid, 1], family="binomial", alpha=1)


fit.vbfactanal1$param$Mu$y
pred <- cbind(vbfactanal1=predict.logvbfactanal(fit.vbfactanal1)$y$pred[, 1],
              vbfactanal2=predict.logvbfactanal(fit.vbfactanal2)$y$pred[, 1],
              glmnet1=predict(fit.glmnet1, x, type="response", 
                              s="lambda.min")[, 1],
              glmnet2=predict(fit.glmnet2, x, type="response", 
                              s="lambda.min")[, 1])


auc <- apply(pred, 2, function(pr) {pROC::auc(ytest[mid, 1], pr[mid])})
pmse <- colMeans((pred[mid, ] - ytest[mid, 1])^2)

pairs(data.frame(true=ytest, pred)[mid, ])
cbind(ytest, 1*(pred > 0.5))[mid, ]

