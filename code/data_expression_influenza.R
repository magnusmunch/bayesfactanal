### libraries
library(affy)
library(GEOquery)
library(genefilter)

### functions
# eset creation function
read_eset <- function(accnr, dir="data", raw.dir="rawdata") {
  
  # load series and platform data from GEO
  if(!getGEOSuppFiles(accnr, baseDir="data", fetch_files=FALSE)$fname %in% 
     list.files(paste0("data/", accnr))) {
    getGEOSuppFiles(accnr, baseDir="data")  
  }
  if(!file.exists(paste0(dir, "/", accnr, "/", raw.dir))) {
    untar(paste(dir, "/",accnr,"/", accnr, "_RAW.tar",sep=""), 
          exdir=paste0(dir, "/", accnr, "/", raw.dir))  
  }
  cels <- list.celfiles(paste0("data/", accnr, "/", raw.dir))
  data <- ReadAffy(filenames=paste0("data/", accnr, "/", raw.dir, "/", cels))
  
  # read measurements times
  con <- gzcon(url(paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", 
                          substr(accnr, 1, 5), "nnn/", accnr, "/matrix/", accnr, 
                          "_series_matrix.txt.gz")))
  txt <- readLines(con)
  times <- read.table(textConnection(txt), skip=41, nrows=1)[1, -1]
  times <- sapply(strsplit(as.character(times), ": "), "[", 2)
  subjects <- read.table(textConnection(txt), skip=40, nrows=1)[, -1]
  subjects <- sapply(strsplit(as.character(subjects), ": "), "[", 2)
  
  # create expression data
  data <- rma(data)
  pData(data) <- cbind(pData(data), subject=as.numeric(subjects),
                       day=as.numeric(gsub("D", "", times)),
                       accnr=accnr)
  return(data)
}

# data set creation function
read_outcomes <- function(accnr) {
  
  # retrieve titer data
  con <- gzcon(url(paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", 
                          substr(accnr, 1, 5), "nnn/", accnr, "/matrix/", accnr, 
                          "_series_matrix.txt.gz")))
  txt <- readLines(con)
  titer <- read.table(textConnection(txt), skip=43, nrows=6)[, -1]
  subjects <- read.table(textConnection(txt), skip=40, nrows=1)[, -1]
  subjects <- sapply(strsplit(as.character(subjects), ": "), "[", 2)
  times <- read.table(textConnection(txt), skip=41, nrows=1)[1, -1]
  times <- sapply(strsplit(as.character(times), ": "), "[", 2)
  
  # remove duplicates (all three microarray times are duplicates of each other)
  titer <- titer[, !duplicated(subjects)]
  
  # combine into one data.frame
  data <- Reduce("rbind", lapply(1:ncol(titer), function(s) {
    t(sapply(strsplit(titer[, s], "hai titer \\(day |\\) - |: "), "[", -1))}))
  hai_factor <- factor(data[, 2], levels=unique(data[, 2]))
  hai_legend <- paste0(1:length(levels(hai_factor)), "='", levels(hai_factor), 
                       "'", collapse=", ")
  data <- data.frame(subject=as.numeric(rep(subjects[!duplicated(subjects)], 
                                            each=nrow(titer))), 
                     day_titer=as.numeric(data[, 1]), 
                     hai=as.numeric(hai_factor), titer=as.numeric(data[, 3]),
                     hai_legend=hai_legend)
  
  # change to wide format and add accnr
  data <- reshape(data, direction="wide", idvar=c("subject", "day_titer"), 
                  timevar="hai", v.names="titer")
  data$day <- 0 + (data$day_titer==28)*3
  data$accnr <- accnr
  return(data)     
  
}

### create data objects
# read data
eset_2007 <- read_eset(accnr="GSE29614")
titer_2007 <- read_outcomes(accnr="GSE29614")
eset_2008 <- read_eset(accnr="GSE29617")
titer_2008 <- read_outcomes(accnr="GSE29617")

# combine two batches
featureNames(eset_2008) <- gsub("PM_", "", featureNames(eset_2008), fixed=TRUE)
annotation(eset_2008) <- annotation(eset_2007)
eset <- combine(eset_2007, eset_2008)
titer <- rbind(titer_2007, titer_2008)

# combine expression and titer data
pData(eset) <- merge(pData(eset), titer, all.x=TRUE)

### transform data
# select and order subjects with measurements on days 0 and 3 
remove <- aggregate(day ~ subject + accnr, data=pData(eset),
                    function(x) {!(any(x==0) & any(x==3))})
remove <- remove[remove$day, -3]
sampleid <- which(apply(pData(eset)[, c("subject", "accnr")], 1, function(m) {
  !any(apply(remove, 1, function(s) {all(m==s)}))}))
eset <- eset[, sampleid]
eset <- eset[, order(pData(eset)$accnr, pData(eset)$subject, pData(eset)$day)]

# calculate max log difference score for titers
logmaxtiter <- aggregate(cbind(titer.1, titer.2, titer.3) ~ subject + accnr,
                         data=pData(eset),
                         FUN=function(x) {x[2]/x[1]})
logmaxtiter <- cbind(logmaxtiter[, c(1, 2)], 
                     logmaxtiter=apply(logmaxtiter[, -c(1, 2)], 1, max))
logmaxtiter <- cbind(logmaxtiter[, c(1, 2)],
                     logmaxtiter=log2(logmaxtiter[, 3]) - 
                       mean(log2(logmaxtiter[, 3])))

# take difference of microarray measurements day 3 an day 0
expr <- aggregate(t(exprs(eset)), 
                  by=list(pData(eset)$subject, pData(eset)$accnr), 
                  FUN=function(x) {x[2] - x[1]})[, -c(1, 2)]
eset.diff <- ExpressionSet(assayData=t(expr),
                           phenoData=AnnotatedDataFrame(logmaxtiter),
                           featureData=featureData(eset),
                           experimentData=experimentData(eset),
                           annotation=annotation(eset))

### filter features
# object to store filter log
filter.log <- list()

# filter features with missing values
var.func <- function(x) {1 - any(is.na(x))}
filter <- nsFilter(eset.diff, var.func=var.func, var.cutoff=0.5, 
                   filterByQuantile=FALSE, var.filter=TRUE, feature.exclude="", 
                   require.entrez=FALSE, require.GOBP=FALSE, require.GOCC=FALSE,
                   require.GOMF=FALSE, require.CytoBand=FALSE, 
                   remove.dupEntrez=FALSE)
eset.filtNA <- filter$eset
filter.log <- c(filter.log, list(filter$filter.log))

# calculate coefficient of variation (incl. non-parametric version) filter
psel <- 200
var.func <- function(x) {abs(mean(x))/sd(x)}
filter <- nsFilter(eset.filtNA, var.func=var.func, 
                   var.cutoff=1 - psel/nrow(exprs(eset.filt)), 
                   filterByQuantile=TRUE, var.filter=TRUE, feature.exclude="", 
                   require.entrez=FALSE, require.GOBP=FALSE, require.GOCC=FALSE,
                   require.GOMF=FALSE, require.CytoBand=FALSE, 
                   remove.dupEntrez=FALSE)
eset.filt <- filter$eset
save(eset.filt, file="data/eset.filt200_expression_influenza.Rdata")
filter.log <- c(filter.log, list(filter$filter.log))
# var.func <- function(x) {
#   abs(median(x))/abs(diff(quantile(x, probs=c(0.25, 0.75))))
# }
# filter <- nsFilter(eset, var.func=var.func,
#                    var.cutoff=1 - psel/nrow(exprs(eset.filt)), 
#                    filterByQuantile=TRUE, var.filter=TRUE, feature.exclude="", 
#                    require.entrez=FALSE, require.GOBP=FALSE, require.GOCC=FALSE,
#                    require.GOMF=FALSE, require.CytoBand=FALSE, 
#                    remove.dupEntrez=FALSE)

psel <- 416
var.func <- function(x) {abs(mean(x))/sd(x)}
filter <- nsFilter(eset.filtNA, var.func=var.func, 
                   var.cutoff=1 - psel/nrow(exprs(eset.filtNA)), 
                   filterByQuantile=TRUE, var.filter=TRUE, feature.exclude="", 
                   require.entrez=FALSE, require.GOBP=FALSE, require.GOCC=FALSE,
                   require.GOMF=FALSE, require.CytoBand=FALSE, 
                   remove.dupEntrez=FALSE)
eset.filt <- filter$eset
save(eset.filt, file="data/eset.filt416_expression_influenza.Rdata")
filter.log <- c(filter.log, list(filter$filter.log))

### save data
rm(list=ls())

