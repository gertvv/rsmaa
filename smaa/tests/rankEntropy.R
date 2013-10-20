# Test the C implementation of ranking entropy (about 40x speed up)
library(smaa)

data <- dget('rankEntropy.txt')
N <- dim(data$ranks)[1]
m <- dim(data$ranks)[2]

counts <- .C("smaa_countRankings", as.integer(t(data$ranks)),
  as.integer(N), as.integer(m),
  counts=as.integer(rep(0, N)),
  NAOK=FALSE, DUP=FALSE)$counts

stopifnot(all.equal(counts[counts > 0], data$counts))
stopifnot(all.equal(smaa.entropy.ranking(data$ranks), data$entropy))
