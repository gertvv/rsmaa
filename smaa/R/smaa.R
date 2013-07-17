smaa.values <- function(meas, pref) {
  N <- dim(meas)[1]
  m <- dim(meas)[2]
  n <- dim(meas)[3]
  stopifnot(identical(dim(pref), c(N, n)) || identical(length(pref), n))

  values <- t(.C("smaa_values",
    as.double(aperm(meas, c(2,3,1))), as.double(t(pref)),
    as.integer(N), as.integer(m), as.integer(n), is.vector(pref),
    v=matrix(0.0, nrow=m, ncol=N))$v)
  dimnames(values) <- dimnames(meas)[1:2]
  class(values) <- "smaa.values"
  values
}

plot.smaa.values <- function(x, ...) {
  boxplot(lapply(apply(x, 2, function(y) { list(y) }), function(z) { unlist(z) }),
    main="Alternative values")
}

smaa.ranks <- function(values) {
  N <- dim(values)[1]
  m <- dim(values)[2]
  ranks <- t(.C("smaa_ranks", as.double(t(values)),
    as.integer(N), as.integer(m),
    r=matrix(0L, nrow=m, ncol=N),
    NAOK=FALSE, DUP=FALSE)$r) + 1
  dimnames(ranks) <- dimnames(values)
  class(ranks) <- "smaa.ranks"
  ranks
}

plot.smaa.ranks <- function(x, ...) {
  plot(smaa.ra(x), ...)
}

print.smaa.ranks <- function(x, ...) {
  print(unclass(x), ...)
}

summary.smaa.ranks <- function(object, ...) {
  smaa.ra(object)
}

smaa.ra <- function(ranks) {
  N <- dim(ranks)[1]
  m <- dim(ranks)[2]

  ra <- t(apply(ranks, 2, tabulate, nbins=m) / N)
  attr(ra, "smaa.N") <- N
  class(ra) <- "smaa.ra"
  ra
}

plot.smaa.ra <- function(x, ...) {
  barplot(t(x), main="Rank acceptabilities", ...)
}

print.smaa.ra <- function(x, ...) {
  cat(paste("Rank acceptabilities (N = ", attr(x, "smaa.N"), " iterations): \n", sep=""))
  attr(x, "smaa.N") <- NULL
  print(unclass(x), ...)
}

smaa.cw <- function(ranks, pref) {
  N <- dim(ranks)[1]
  m <- dim(ranks)[2]
  n <- dim(pref)[2]
  stopifnot(identical(dim(pref), c(N, n)))

  cw <- t(apply(ranks, 2, function(r) {
    apply(pref[r == 1, ], 2, mean)
  }))
  attr(cw, "smaa.N") <- N
  class(cw) <- "smaa.cw"
  cw
}

plot.smaa.cw <- function(x, ...) {
  # FIXME: use layout() instead?
  par(mar=c(8.1, 4.1, 4.1, 8.1))
  plot(NA, xlim=c(1, ncol(x)), ylim=c(0, max(x)), xlab="", ylab="Weight", xaxt='n', bty='L',
    main='Central weights', ...)
  for (i in 1:nrow(x)) {
    lines(x[i,], pch=i, type="b")
  }
  axis(side=1, at=1:ncol(x), labels=colnames(x), las=2)

  legend("topright", inset=c(-0.25,0), legend=rownames(x), pch=(1:nrow(x)), xpd=TRUE)
}

print.smaa.cw <- function(x, ...) {
  cat(paste("Central weights (N = ", attr(x, "smaa.N"), " iterations): \n", sep=""))
  attr(x, "smaa.N") <- NULL
  print(unclass(x), ...)
}

smaa.entropy.ranking <- function(ranks, p0=1) {
  N <- dim(ranks)[1]

  p <- rle(sort(apply(ranks, 1, paste, collapse=".")))$lengths / N * p0
  -sum(p * log2(p))
}

smaa.entropy.choice <- function(ra, p0=1) {
  if (class(ra) == 'smaa.ranks') { ra <- smaa.ra(ra) }
  stopifnot(class(ra) == 'smaa.ra')

  p <- ra[, 1] * p0 # first-rank acceptabilities
  -sum(p * log2(p))
}

smaa.cf <- function(meas, cw) {
  N <- dim(meas)[1]
  m <- dim(meas)[2]
  n <- dim(meas)[3]
  stopifnot(identical(dim(cw), c(m, n)))

  cf <- diag(apply(cw, 1, function(w) {
    if (all(!is.na(w))) {
      smaa.ra(smaa.ranks(smaa.values(meas, w)))[,1]
    } else {
      NA
    }
  }))

  result <- list(cf=cf, cw=cw)
  attr(result, "smaa.N") <- N
  class(result) <- "smaa.cf"

  result
}

print.smaa.cf <- function(x, ...) {
  cat(paste("Central weights (N = ", attr(x$cw, "smaa.N"), " iterations) and\n",
    "  confidence factors (N = ", attr(x, "smaa.N"), " iterations): \n", sep=""))
  attr(x, "smaa.N") <- NULL
  print(cbind(x$cf, x$cw), ...)
}

plot.smaa.cf <- function(x, ...) {
  plot(x$cw)
}

smaa <- function(meas, pref) {
  N <- dim(meas)[1]
  m <- dim(meas)[2]
  n <- dim(meas)[3]
  stopifnot(identical(dim(pref), c(N, n)) || identical(length(pref), n))

  result <- .C("smaa", as.double(aperm(meas, c(2,3,1))), as.double(t(pref)),
    as.integer(N), as.integer(m), as.integer(n), is.vector(pref),
    h=matrix(0, nrow=m, ncol=m),
    cw=matrix(0, nrow=m, ncol=n),
    NAOK=FALSE, DUP=FALSE)

  # Introduce NAs where central weights are undefined
  cw <- t(apply(result$cw, 1, function(w) { if (sum(w) < 0.5) rep(NA, length(w)) else w }))

  rownames(cw) <- dimnames(meas)[[2]]
  colnames(cw) <- dimnames(meas)[[3]]
  class(cw) <- "smaa.cw"
  attr(cw, "smaa.N") <- N

  ra <- result$h / N
  rownames(ra) <- dimnames(meas)[[2]]
  class(ra) <- "smaa.ra"
  attr(ra, "smaa.N") <- N

  result <- list(cw=cw, ra=ra)
  class(result) <- "smaa.result"
  result
}

print.smaa.result <- function(x, ...) {
  print(x$ra)
  cat("\n")
  print(x$cw)
}

plot.smaa.result <- function(x, ...) {
  plot(x$ra, ...)
}
