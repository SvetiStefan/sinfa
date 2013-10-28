## sinfa.R: some experiments trying to reproduce the SINFA procedure from Wang and Bilger (1973)
## Copyright (c) 2013 by David A. van Leeuwen

## This program is free software: you can redistribute it and/or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, version 3 of the License.

##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.

##     You should have received a copy of the GNU General Public License
##     along with this program.  If not, see <http://www.gnu.org/licenses/>.

## version: 0.1 (this is not a git version)

## To understand what this does, read Wang an Bilger, JASA vol 54(5), 1973, and type
# source("sinfa.R")
# read.wang.bilger()
# x
# f
# s <- sinfa(x, f, 7)
# s
# summary(s)


read.wang.bilger <- function() {
  assign("x", read.csv("wang-bilger-table6.csv", row.names=1), envir=.GlobalEnv)
  assign("f", read.csv("wang-bilger-features.csv", row.names=1), envir=.GlobalEnv)
  cons <- colnames(x)
}

## cond can be multi-dimensional, i.e., an array of nphones x nfeatures
## cond can also be NULL, then this is just entropy()
cond.entropy <- function(x, feature, cond) {
   pm <- apply(x, 1, sum)/sum(x)
   f <- cbind(feature,cond)             # combine feature and condition, 1st is feature
   nph <- nrow(f)                       # number of phonemes
   nf <- ncol(f)
   for (i in 1:nf) 
     f[,i] <- tapply(numeric(nph), f[,i]) # normalize levels to 1..n, keep dim attrib
   dims <- apply(f, 2, max)             # total number of levels
   p <- array(dim=dims)
   for (i in 1:prod(dims)) {
     ind <- arrayInd(i,dims)
     sel <- rep(T, nph)                 # the phonemes in the current selection
     for (fi in seq(along=ind))
       sel <- sel & f[,fi]==ind[fi]     # select the right phones from this feature
     p[i] <- sum(pm[sel])
   }
   if (nf>1) {                     # condition actually specified...
     pcond <- apply(p, 2:nf, sum)  # marginal over feature or interest
     lp <- log(sweep(p, 2:nf, pcond, '/'), base=2)
   } else lp <- log(p, base=2)                 # same as entropy()
   lp[is.nan(lp) | lp==-Inf] <- 0       # lot of 0/0's out there
   return(-sum(p*lp))
 }

## multidimensional condditional mutual entropy.  cond may be NULL
cond.mutual <- function(x, feature, cond=NULL) {
  f <- cbind(feature,feature,cond) # combine feature and condition, 1-2 are feature
  nph <- nrow(f)                       # number of phonemes
  nf <- ncol(f)
  for (i in 1:nf) 
    f[,i] <- tapply(numeric(nph), f[,i]) # normalize levels to 1..n, keep dim attrib
  dims <- apply(f, 2, max)             # total number of levels
  p <- array(0, dim=dims)
  n <- sum(x)
  for (i in 1:prod(dims)) {
    ind <- arrayInd(i,dims)
    sel <- rep(T, nph)                 # the phonemes in the current selection
    for (fi in seq(from=3, length.out=nf-2)) 
      sel <- sel & f[,fi]==ind[fi]     # select the right phones from this feature
    y <- x[sel & f[,1]==ind[1], sel & f[,2]==ind[2]]
    if (prod(dim(y))>0) p[i] <- sum(y)/n
  }
  if (nf>2) {                           # condition was specified
    pz <- apply(p, 3:nf, sum)
    pxz <- apply(p, c(1, 3:nf), sum)
    pyz <- apply(p, 2:nf, sum)
    arg <- sweep(p, 3:nf, pz, '*')
    arg <- sweep(arg, c(1, 3:nf), pxz, '/')
    arg <- sweep(arg, 2:nf, pyz, '/')
  } else {
    px <- apply(p, 1, sum)
    py <- apply(p, 2, sum)
    arg <- sweep(p, 1, px, '/')
    arg <- sweep(arg, 2, py, '/')
  }
  lp <- log(arg, base=2)
  lp[is.nan(lp) | lp==-Inf] <- 0       # lot of 0/0's out there
  return(sum(p*lp))
}


# run a sinfa for confusion matrix x and feature data frame f
sinfa <- function(x, f, nmax=1) {
  stopifnot(all(row.names(x)==row.names(f)))
  f.names <- names(f)                   # current list of features to tabulate
  cond.f <- NULL
  cond.names <- NULL
##  cond.pasted.names <- "nothing"
  g <- function(fea) c(ent=cond.entropy(x, f[[fea]], cond.f),
                       mut=cond.mutual(x, f[[fea]], cond.f))
  n <- 1
  tot <- list()
  while(n<=nmax) {
    res <- data.frame(t(sapply(f.names, g)))
    if (sum(res$ent)==0) {
      warning("No more entropy left\n")
      break
    }
    no.info <- res$ent==0
    res$rel <- res$mut / res$ent
    res$rel[no.info] <- NA
#    if (verbose) {
#      cat(sprintf("Iteration %d: keeping %s constant\n", n, cond.names))
#      print(res)
#    }
    max.fea <- which(res$rel == max(res$rel,na.rm=T))
    tot[[n]] <- list(table=res[!no.info,], max=f.names[max.fea],
                     no.info=f.names[no.info])
    cond.names <- append(cond.names, f.names[max.fea])
    f.names <- f.names[-c(max.fea,no.info)]
    cond.f <- sapply(cond.names, function(x) f[[x]])
##    cond.pasted.names <- paste(cond.names, collapse=":")
##    if (verbose)
##      cat(sprintf("Max at feature '%s'\n\n", cond.names[length(cond.names)]))
    n <- n+1
  }
  class(tot) <- "sinfa"
  tot
}

## default way of printing the SINFA analysis
print.sinfa <- function(x) {
  stopifnot(is(x, "sinfa"))
  res <- NULL
  cond.pasted.names <- "nothing"
  cond.names <- NULL
  n <- 1
  for (iter in x) {
    cat(sprintf("Iteration %d: knowing %s\n", n, cond.pasted.names))
    print(iter$table)
    if (length(iter$no.info)) {
      cat(sprintf("No more information from '%s'\n", iter$no.info))
      cond.names <- append(cond.names, paste(iter$no.info, "=0", sep=""))
    }
    maxinfo <- paste(iter$max, collapse="/")
    cat(sprintf("Max relative information at feature '%s'\n\n", maxinfo))
    cond.names <- append(cond.names, maxinfo)
    cond.pasted.names <- paste(cond.names, collapse=":")
    n <- n+1
  }
}

## Summary of the analysis, just show the line of maximum information features.       
summary.sinfa <- function(x) {
  stopifnot(is(x, "sinfa"))
  res <- NULL
  for (iter in x) {
    m <- iter$max[1]
    line <- iter$table[m,]
    row.names(line) <- paste(iter$max, collapse="/")
    res <- rbind(res, line)
  }
  transform(res, cum.ent=cumsum(ent), cum.mut=cumsum(mut))
}

read.natalie <- function(av="ao", noise="noise") {
  global <- function(x, y) assign(x, y, env=.GlobalEnv)
  f <- read.table("data/features.txt")
  f <- data.frame(t(f[row.names(f)!="vcless",]))                 # transpose, remove vcless
  global("conditions", c("bao", "bat", "con", "hel", "hoo", "niq", "rub", "sur", "tap"))
  g <- function(x) {
    x <- read.csv(sprintf("data/cm_%s_%s - %s.csv", av, noise, toupper(x)), row.names=1)
    x[-nrow(x),-ncol(x)]           # ditch totals
  }
  x <- lapply(conditions, g)
  names(x) <- conditions
  global("x", x)
  global("f", f)
}

