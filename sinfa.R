## some experiments trying to reproduce the SINFA procedure from Wang and Bilger (1973)
## (c) 2013 David A. van Leeuwen

read.wang.bilger <- function() {
  assign("x", read.csv("wang-bilger-table6.csv", row.names=1), envir=.GlobalEnv)
  assign("f", read.csv("wang-bilger-features.csv", row.names=1), envir=.GlobalEnv)
  cons <- colnames(x)
}

read.natalie <- function() {
  global <- function(x, y) assign(x, y, env=.GlobalEnv)
  f <- read.table("data/features.txt")
  f <- data.frame(t(f))                 # transpose
  global("conditions", c("bao", "bat", "con", "hel", "hoo", "niq", "rub", "sur", "tap"))
  g <- function(x) {
    x <- read.csv(sprintf("data/cm_ao_noise - %s.csv", toupper(x)), row.names=1)
    x[-nrow(x),-ncol(x)]           # ditch totals
  }
  x <- lapply(conditions, g)
  names(x) <- conditions
  global("x", x)
  global("f", f)
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
run.sinfa <- function(x, f, nmax=1) {
  stopifnot(all(row.names(x)==row.names(f)))
  f.names <- names(f)                   # current list of features to tabulate
  cond.f <- NULL
  cond.names <- NULL
  cond.pasted.names <- "nothing"
  g <- function(fea) c(ent=cond.entropy(x, f[[fea]], cond.f),
                       mut=cond.mutual(x, f[[fea]], cond.f))
  n <- 1
  while(n<=nmax) {
    res <- data.frame(t(sapply(f.names, g)))
    res$rel <- res$mut / res$ent
    cat(sprintf("Iteration %d: keeping %s constant\n", n, cond.pasted.names))
    print(res)
    max.fea <- which.max(res$rel)
    cond.names <- append(cond.names, f.names[max.fea])
    f.names <- f.names[-max.fea]
    cond.f <- sapply(cond.names, function(x) f[[x]])
    cond.pasted.names <- paste(cond.names, collapse=":")
    cat(sprintf("Max at feature '%s'\n\n", cond.names[length(cond.names)]))
    n <- n+1
  }
}

## from here onwards old code
##
##


entropy <- function(x, feature, margin=1) {
  pm <- apply(x, margin, sum)/sum(x)
  p <- tapply(pm, feature, sum)
  lp <- log(p + (p==0), base=2)           # fix log
  return(-sum(p*lp))
}

cond.entropy1 <- function(x, feature, cond, margin=1) {
  pm <- apply(x, margin, sum)/sum(x)
  index.cond <- tapply(numeric(length(cond)), cond)
  nlevels.cond <- max(index.cond)
  index.feature <- tapply(numeric(length(feature)), feature)
  nlevels.feature <- max(index.feature)
  p <- matrix(nrow=nlevels.cond, ncol=nlevels.feature)
  for (i in 1:nlevels.cond) for (j in 1:nlevels.feature) {
    p[i,j] <- sum(pm[index.cond==i & index.feature==j])
  }
  pi <- tapply(pm, cond, sum)
  arg <- p / matrix(pi, nrow=nlevels.cond, ncol=nlevels.feature, byrow=F) 
  lp <- log(arg + (arg==0), base=2)
  return(-sum(p*lp))
}


mutual <- function(x, feature) {
  index <- tapply(numeric(length(feature)), feature)
  nlevels <- max(index)
  p <- matrix(nrow=nlevels, ncol=nlevels)
  n <- sum(x)
  for (i in 1:nlevels) for (j in 1:nlevels) {
    p[i,j] <- sum(x[index==i,index==j])/n
  }
  pi <- apply(p, 1, sum)
  pj <- apply(p, 2, sum)
  arg <- outer(pi, pj, '*')/p
  lp <- log(arg, base=2)
  lp[is.nan(lp) | lp==Inf] <- 0         # fix p->0 cases
  return(-sum(p*lp))
}

rel.transmission <- function(x, feature, margin=1) 
  return(mutual(x, feature) / entropy(x, feature, margin))

cond.mutual1 <- function(x, feature, cond) {
  index.cond <- tapply(numeric(length(cond)), cond)
  nlevels.cond <- max(index.cond)
  index.feature <- tapply(numeric(length(feature)), feature)
  nlevels.feature <- max(index.feature)
  p <- array(0, dim=c(nlevels.cond,nlevels.feature,nlevels.feature))
  n <- sum(x)
  for (i in 1:nlevels.cond) for (j in 1:nlevels.feature) for (k in 1:nlevels.feature) {
    y <- x[index.cond==i & index.feature==j, index.cond==i & index.feature==k]
    if (prod(dim(y))>0) 
      p[i,j,k] <- sum(y)/n
  }
  pi <- apply(p, 1, sum)
  pij <- apply(p, c(1,2), sum)
  pik <- apply(p, c(1,3), sum)
  arg <- sweep(p, 1, pi, '*')
  arg <- sweep(arg, c(1:2), pij, '/')
  arg <- sweep(arg, c(1,3), pik, '/')
  lp <- log(arg, base=2)
  lp[is.nan(lp) | lp==-Inf] <- 0                   # lot of 0/0's out there
  return (sum(p * lp))
}
  
# alternative computation, same result. 
cond.entropy2 <- function(x, feature, cond, margin=1) {
  pm <- apply(x, margin, sum)/sum(x)
  index.cond <- tapply(numeric(length(cond)), cond)
  nlevels.cond <- max(index.cond)
  index.feature <- tapply(numeric(length(feature)), feature)
  nlevels.feature <- max(index.feature)
  p <- matrix(nrow=nlevels.cond, ncol=nlevels.feature)
  for (i in 1:nlevels.cond) for (j in 1:nlevels.feature) {
    p[i,j] <- sum(pm[index.cond==i & index.feature==j]) / sum(pm[index.cond==i])
  }
  pi <- tapply(pm, cond, sum)
  lp <- log(p + (p==0), base=2)
  ent <- apply(p*lp, 1, sum)
  return(-sum(pi*ent))
}

# this function computes outer products varying only the first dimension,
# while keeping the remaining dimensions constant
partial.outer <- function(x, y, FUN) {
  nd <- length(dim(x))
  d.const <- dim(x)[2:nd]
  stopifnot(nd == length(dim(y)), all(d.const==dim(y)[2:nd]))
  nd.const <- prod(d.const)
  x <- array(x, dim=c(dim(x)[1], nd.const))
  y <- array(y, dim=c(dim(y)[1], nd.const))
  res <- array(NA, dim=c(dim(x)[1], dim(y)[1], nd.const))
#  print(sprintf("Output dim %d, %d, %d", dim(x)[1], dim(y)[1], d.const))
  for (i in 1:nd.const)
    res[,,i] <- outer(x[,i], y[,i], FUN)
  dim(res) <- c(dim(x)[1], dim(y)[1], d.const)
  res                
}

## features
#voice <- cons %in% c("b", "d", "g", "v", "dh", "z", "zj", "dzj")
#fric <- cons %in% c("f", "th", "s", "sh", "v", "dh", "z", "zj", "tsh", "dzj")
#sib <- cons %in% c("s", "sh", "z", "zj", "tsh", "dzj")
#str <- cons %in% c("f", "s", "sh", "v", "z", "zj", "tsh", "dzj")
#dur <- cons %in% c("s", "sh", "z", "zj")
#pl1 <- c(0, 1, 2, 0, 1, 2, 0, 1, 1, 2, 0, 1, 1, 2, 2, 2)
#high <- cons %in% c("k", "g", "sh", "zj", "tsh", "dzj")
