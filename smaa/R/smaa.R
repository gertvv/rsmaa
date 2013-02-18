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
	h <- matrix(0, nrow=m, ncol=m) # rank frequency table 
	wc <- matrix(0, nrow=m, ncol=n) # central weights
	for (k in 1:N) {
		w <- pref[k,]
		x <- meas[k,,]
		t <- x %*% w # calculate additive linear value function

		r <- rank(-t, ties.method='min') # rank alternatives
		for (i in 1:m) {
			h[i, r[i]] <- h[i, r[i]] + 1 # update rank counts
			if (r[i] == 1) {
				wc[i,] <- wc[i,] + w # update central weights
			}
		}
	}
	for (i in 1:m) {
		if (h[i, 1] > 0) {
			wc[i,] <- wc[i,] / h[i, 1] # normalize central weight
		}
	}
	list(rankAcc=h/N, centralWeights=wc)
}

#dyn.load("smaa.so")

smaa <- function(m, n, meas, pref, N=10000) {
	result <- .C("smaa", as.double(aperm(meas, c(2,3,1))), as.double(t(pref)),
		as.integer(N), as.integer(m), as.integer(n),
		h=matrix(0.0, nrow=m, ncol=m),
		wc=matrix(0, nrow=m, ncol=n),
		NAOK=FALSE, DUP=FALSE)

	list(rankAcc=result$h/N, centralWeights=result$wc)
}
