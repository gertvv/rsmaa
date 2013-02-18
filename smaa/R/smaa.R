# Calculate SMAA rank acceptability indices and central weights
# N: no. of iterations (default 10,000)
# m: no. of alternatives
# n: no. of criteria
# meas: criteria measurements: N.m.n array, where meas[i,,] is a matrix where
# the m alternatives are the rows, the n criteria the columns. The values must
# be *standardized* measurements (i.e. after application of the partial value
# function).
# pref: weights: N.n array, where pref[i,] is a (normalized) weight vector
smaa <- function(m, n, meas, pref, N=10000) {
	result <- .C("smaa", as.double(aperm(meas, c(2,3,1))), as.double(t(pref)),
		as.integer(N), as.integer(m), as.integer(n),
		h=matrix(0.0, nrow=m, ncol=m),
		wc=matrix(0, nrow=m, ncol=n),
		NAOK=FALSE, DUP=FALSE)

	list(rankAcc=result$h/N, centralWeights=result$wc)
}
