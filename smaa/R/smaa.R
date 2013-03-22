smaa.values <- function(meas, pref) {
	N <- dim(meas)[1]
	m <- dim(meas)[2]
	n <- dim(meas)[3]
	stopifnot(identical(dim(pref), c(N, n)))

	values <- t(.C("smaa_values",
		as.double(aperm(meas, c(2,3,1))), as.double(t(pref)),
		as.integer(N), as.integer(m), as.integer(n),
		v=matrix(0.0, nrow=m, ncol=N))$v)
	dimnames(values) <- dimnames(meas)[1:2]
	values
}

smaa.ranks <- function(values) {
	N <- dim(values)[1]
	m <- dim(values)[2]
	ranks <- t(.C("smaa_ranks", as.double(t(values)),
		as.integer(N), as.integer(m),
		r=matrix(0L, nrow=m, ncol=N),
		NAOK=FALSE, DUP=FALSE)$r) + 1
	dimnames(ranks) <- dimnames(values)
	ranks
}

smaa.ra <- function(ranks) {
	N <- dim(ranks)[1]
	m <- dim(ranks)[2]

	t(apply(ranks, 2, tabulate, nbins=m) / N)
}

smaa.cw <- function(ranks, pref) {
	N <- dim(ranks)[1]
	m <- dim(ranks)[2]
	n <- dim(pref)[2]
	stopifnot(identical(dim(pref), c(N, n)))

	t(apply(ranks, 2, function(r) {
		apply(pref[r == 1, ], 2, mean)
	}))
}

smaa.entropy.ranking <- function(ranks, p0=1) {
	N <- dim(ranks)[1]

	p <- rle(sort(apply(ranks, 1, paste, collapse=".")))$lengths / N * p0
	-sum(p * log2(p))
}

smaa.entropy.choice <- function(ra, p0=1) {
	# FIXME: could implement this if I remembered whether row or columns are
	# the alternatives.
}

smaa <- function(meas, pref, m=dim(meas)[2], n=dim(meas)[3], N=dim(meas)[1], generate.values=FALSE) {
	stopifnot(identical(dim(meas), c(N, m, n)))
	stopifnot(identical(dim(pref), c(N, n)))

	t <- numeric(0)
	if (generate.values) {
		t <- matrix(0, nrow=m, ncol=N)
	}

	result <- .C("smaa", as.double(aperm(meas, c(2,3,1))), as.double(t(pref)),
		as.integer(N), as.integer(m), as.integer(n),
		as.integer(generate.values),
		h=matrix(0, nrow=m, ncol=m),
		cw=matrix(0, nrow=m, ncol=n),
		t=t,
		NAOK=FALSE, DUP=FALSE)

	# Introduce NAs where central weights are undefined
	result$cw <- t(apply(result$cw, 1, function(w) { if (sum(w) < 0.5) rep(NA, length(w)) else w }))

	rownames(result$cw) <- dimnames(meas)[[2]]
	colnames(result$cw) <- dimnames(meas)[[3]]

	rownames(result$h) <- dimnames(meas)[[2]]

	if (generate.values) {
		rownames(result$t) <- dimnames(meas)[[2]]
	}

	list(rankAcc=result$h/N, centralWeights=result$cw, values=t(result$t))
}
