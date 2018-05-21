\name{smaa.pwi}
\alias{smaa.pwi}
\title{
SMAA pair-wise winning indices
}
\description{
Calculate SMAA pair-wise winning indices from sampled rankings.
}
\usage{
smaa.pwi(ranks)
}
\arguments{
  \item{ranks}{An \eqn{N \times m}{N * m} array of sampled rankings, where \eqn{N}{N} is the number
of SMAA iterations and \eqn{m}{m} is the number of alternatives.}
}
\value{
An \eqn{m \times m}{m * m} matrix of pair-wise winning indices. The
index at (i, j) describes the share of samples for which alternative i
has a better (lower) rank than alternative j.
}
\author{
Tommi Tervonen
}
\seealso{
\code{\link{smaa.ranks}}
}
\examples{
N <- 1E4; m <- 2; n <- 3
meas <- dget(system.file("extdata/thrombo-meas.txt.gz", package="smaa"))
pref <- dget(system.file("extdata/thrombo-weights-nopref.txt.gz", package="smaa"))

# Calculate pair-wise winning indices
values <- smaa.values(meas, pref)
ranks <- smaa.ranks(values)
pwi <- smaa.pwi(ranks)
print(pwi)
}
