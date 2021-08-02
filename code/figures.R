# res <- res9
# methods=c("true", "ridge", "lasso", "FMradio",
#           "emfactanal", "vbfactanal",
#           "ebfactanal", "null")
# m=NULL; plotid=NULL; measure="median";
# legend=TRUE; col=NULL; lty=NULL; pch=NULL; labels=NULL;
# lwd=1; cex.lab=1; cex.axis=1; cex=1;
# cex.main=1; cex.sub=1
# EMSE and PMSE figures
simfig1 <- function(res, methods=c("true", "ridge", "lasso", "FMradio",
                                   "emfactanal", "vbfactanal", 
                                   "ebfactanal", "null"),
                    m=NULL, plotid=NULL, measure="median",
                    legend=TRUE, col=NULL, lty=NULL, pch=NULL, labels=NULL, 
                    lwd=1, cex.lab=1, cex.axis=1, cex=1,
                    cex.main=1, cex.sub=1) {
  
  # set some visual
  if(is.null(col)) {
    col <- brewer.pal(n=length(methods), name="Dark2")
  }
  if(is.null(lty)) {
    lty <- rep(1, length(methods))
  }
  if(is.null(pch)) {
    pch <- rep(1, length(methods))
  }
  if(is.null(labels)) {
    labels <- methods
  }
  
  # calculate number of fits
  lengths <- sapply(methods, function(mt) {
    sum(substr(colnames(res$emse), 1, nchar(mt))==mt)})
  if(is.null(m)) {
    m <- 1:max(lengths)
  }
  if(is.null(plotid)) {
    plotid <- 1:length(m)
  }
  
  # create plot tables
  temse <- lapply(methods, function(mt) {
    vals <- apply(res$emse[, substr(colnames(res$emse), 1, nchar(mt))==mt, 
                           drop=FALSE], 2, measure, na.rm=TRUE)
    cbind(m, rep(vals, length.out=length(m)))[plotid, , drop=FALSE]})
  tpmse <- lapply(methods, function(mt) {
    vals <- apply(res$pmse[, substr(colnames(res$pmse), 1, nchar(mt))==mt, 
                           drop=FALSE], 2, measure, na.rm=TRUE)
    cbind(m, rep(vals, length.out=length(m)))[plotid, , drop=FALSE]})
  tpcor <- lapply(methods, function(mt) {
    vals <- apply(res$pcor[, substr(colnames(res$pcor), 1, nchar(mt))==mt, 
                           drop=FALSE], 2, measure, na.rm=TRUE)
    cbind(m, rep(vals, length.out=length(m)))[plotid, , drop=FALSE]})
  
  # create plots
  opar <- par(no.readonly=TRUE)
  par(mar=opar$mar*c(1, 1.3, 1, 1))
  layout(matrix(c(1:3), nrow=1, ncol=3, byrow=TRUE))
  plot(temse[[1]], type="l", main="a)", xlab="Number of unlabeled", 
       ylab="EMSE", ylim=range(sapply(temse, function(x) {x[, 2]}), na.rm=TRUE), 
       col=col[1], lty=lty[1], lwd=lwd, cex.lab=cex.lab, cex.axis=cex.axis, 
       cex=cex, cex.main=cex.main, cex.sub=cex.sub)
  for(mt in 2:length(methods)) {
    if(length(na.omit(temse[[mt]][, 2]))==1) {
      points(na.omit(temse[[mt]]), col=col[mt], pch=pch[mt], cex.lab=cex.lab,
             cex.axis=cex.axis, cex=cex, cex.main=cex.main, cex.sub=cex.sub)
    } else {
      lines(temse[[mt]], col=col[mt], lty=lty[mt], lwd=lwd, cex.lab=cex.lab, 
            cex.axis=cex.axis, cex=cex, cex.main=cex.main, cex.sub=cex.sub)  
    }
  }
  plot(tpmse[[1]], type="l", main="b)", xlab="Number of unlabeled", 
       ylab="PMSE", ylim=range(sapply(tpmse, function(x) {x[, 2]}), na.rm=TRUE), 
       col=col[1], lty=lty[1], lwd=lwd, cex.lab=cex.lab, cex.axis=cex.axis, 
       cex=cex, cex.main=cex.main, cex.sub=cex.sub)
  for(mt in 2:length(methods)) {
    if(length(na.omit(tpmse[[mt]][, 2]))==1) {
      points(na.omit(tpmse[[mt]]), col=col[mt], pch=pch[mt], cex.lab=cex.lab,
             cex.axis=cex.axis, cex=cex, cex.main=cex.main, cex.sub=cex.sub)
    } else {
      lines(tpmse[[mt]], col=col[mt], lty=lty[mt], lwd=lwd, cex.lab=cex.lab, 
            cex.axis=cex.axis, cex=cex, cex.main=cex.main, cex.sub=cex.sub)
    }
  }
  plot(tpcor[[1]], type="l", main="c)", xlab="Number of unlabeled", 
       ylab=expression("corr(y,"~hat(y)~")"), 
       ylim=range(sapply(tpcor, function(x) {x[, 2]}), na.rm=TRUE), 
       col=col[1], lty=lty[1], lwd=lwd, cex.lab=cex.lab, cex.axis=cex.axis, 
       cex=cex, cex.main=cex.main, cex.sub=cex.sub)
  for(mt in 2:length(methods)) {
    if(length(na.omit(tpcor[[mt]][, 2]))==1) {
      points(na.omit(tpcor[[mt]]), col=col[mt], pch=pch[mt], cex.lab=cex.lab,
             cex.axis=cex.axis, cex=cex, cex.main=cex.main, cex.sub=cex.sub)
    } else {
      lines(tpcor[[mt]], col=col[mt], lty=lty[mt], lwd=lwd, cex.lab=cex.lab, 
            cex.axis=cex.axis, cex=cex, cex.main=cex.main, cex.sub=cex.sub)
    }
  }
  if(legend) {
    legend("bottomright", legend=labels, col=col, lty=lty, pch=pch, lwd=lwd, 
           cex=cex, bg="white", merge=TRUE)
  }
  par(opar)
}

# res <- res9
# methods=c("ebfactanal");
# m=NULL; measure="median";
# legend=TRUE; col=NULL; lty=NULL; labels=NULL; 
# lwd=1; cex.lab=1; cex.axis=1; cex=1;
# cex.main=1; cex.sub=1
# log multiplier figure
simfig2 <- function(res, methods=c("ebfactanal"),
                    m=NULL, measure="median",
                    legend=TRUE, col=NULL, lty=NULL, labels=NULL, 
                    lwd=1, cex.lab=1, cex.axis=1, cex=1,
                    cex.main=1, cex.sub=1, legend.pos="topright") {
  
  # find number of groups
  G <- length(unique(sapply(strsplit(colnames(res$gamma), ".", fixed=TRUE), "[", 
                            2)))
  
  # calculate number of fits
  lengths <- sapply(methods, function(mt) {
    sum(substr(colnames(res$gamma), 1, nchar(mt))==mt)})/G
  
  # set some visual
  if(is.null(col)) {
    col <- brewer.pal(n=max(3, G), name="Dark2")
  }
  if(is.null(lty)) {
    lty <- c(1:(length(methods) + 1)) + 1
  }
  if(is.null(labels)) {
    labels <- c(methods, paste0("group ", 1:G))
  }
  
  # set number of unlabeled
  if(is.null(m)) {
    m <- 1:max(lengths)
  }
  
  # create plot tables
  tlgamma <- lapply(methods, function(mt) {
    lapply(1:G, function(g) {
      vals <- apply(
        log(res$gamma[, substr(colnames(res$gamma), 1, nchar(mt))==mt & 
                        sapply(strsplit(colnames(res$gamma), ".", fixed=TRUE), 
                               "[", 2)==as.character(g), drop=FALSE]), 
        2, measure)
      cbind(m, vals)})})
  
  # create plots
  opar <- par(no.readonly=TRUE)
  par(mar=opar$mar*c(1, 1.5, 1, 1))
  for(mt in 1:length(methods)) {
    for(g in 1:G) {
      if(mt==1 & g==1) {
        plot(tlgamma[[mt]][[g]], type="l", xlab="Number of unlabeled", 
             ylab=expression(log~hat(gamma)[g]),
             ylim=range(sapply(tlgamma, function(x) {
               sapply(x, function(y) {y[, 2]})}), na.rm=TRUE), col=col[g], 
             lty=lty[mt], lwd=lwd, cex.lab=cex.lab, cex.axis=cex.axis, cex=cex,
             cex.main=cex.main, cex.sub=cex.sub)
      } else {
        lines(tlgamma[[mt]][[g]], col=col[g], lty=lty[mt], 
              lwd=lwd, cex.lab=cex.lab, cex.axis=cex.axis, cex=cex,
              cex.main=cex.main, cex.sub=cex.sub)
      }
    }
  }
  abline(a=0, b=0, lty=lty[length(methods) + 1])
  if(legend) {
    legend(legend.pos, legend=c("Bayes", labels), 
           col=c(rep(1, length(methods) + 1), col[1:G]), 
           lty=c(rev(lty), rep(0, G)), 
           fill=c(rep(0, length(methods) + 1), col[1:G]),
           border=c(rep(0, length(methods) + 1), rep("white", G)),
           lwd=lwd, cex=cex, bg="white", merge=TRUE)
  }
  par(opar)
}