smaa <- function(meas, pref, m=dim(meas)[2], n=dim(meas)[3], N=dim(meas)[1], generate.values=FALSE) {
	stopifnot(identical(dim(meas), c(N, m, n)))
	stopifnot(identical(dim(pref), c(N, n)))

	t <- NULL
	if (generate.values) {
		t <- matrix(0, nrow=m, ncol=N)
	}

	result <- .C("smaa", as.double(aperm(meas, c(2,3,1))), as.double(t(pref)),
		as.integer(N), as.integer(m), as.integer(n),
		h=matrix(0.0, nrow=m, ncol=m),
		wc=matrix(0, nrow=m, ncol=n),
		t=t,
		NAOK=FALSE, DUP=FALSE)

	list(rankAcc=result$h/N, centralWeights=result$wc, values=result$t)
}
