library(smaa)

meas <- dget('hansen-meas.txt')
pref <- dget('hansen-weights-ord.txt')

expect.ra <- dget('hansen-ra-ord.txt')
expect.cf <- dget('hansen-cf-ord.txt')

result <- smaa(meas, pref)
cf <- smaa.cf(meas, result$cw)

stopifnot(all.equal(result$ra, expect.ra))
stopifnot(all.equal(cf, expect.cf))
