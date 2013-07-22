library(smaa)

# First alternative beats second on all criteria
# --> NA central weights
N <- 1E4; m <- 2; n <- 3
meas <- array(dim=c(N,m,n))
meas[,1,] <- 1
meas[,2,] <- 0
pref <- matrix(1/3, nrow=N, ncol=n)

result <- smaa(meas, pref)
stopifnot(all.equal(result$cw[1,], rep(1/3, 3)))
stopifnot(all(is.na(result$cw[2,])))

ranks <- smaa.ranks(smaa.values(meas, pref))
cw <- smaa.cw(ranks, pref)
stopifnot(all.equal(cw[1,], rep(1/3, 3)))
stopifnot(all(is.na(cw[2,])))


# First alternative beats second in a single iteration
# --> R gives us a vector instead of a matrix (gee, thanks!)
meas[1,1,] <- 0
meas[1,2,] <- 1

result <- smaa(meas, pref)
stopifnot(all.equal(result$cw[1,], rep(1/3, 3)))
stopifnot(all.equal(result$cw[2,], rep(1/3, 3)))

ranks <- smaa.ranks(smaa.values(meas, pref))
cw <- smaa.cw(ranks, pref)
stopifnot(all.equal(cw[1,], rep(1/3, 3)))
stopifnot(all.equal(cw[2,], rep(1/3, 3)))
