\name{smaa.cf}
\alias{smaa.cf}
\title{
SMAA confidence factors
}
\description{
Calculate SMAA confidence factors of the central weights.
}
\usage{
smaa.cf(meas, cw)
}
\arguments{
  \item{meas}{Criteria measurements. An \eqn{N \times m \times n}{N * m * n}
array, where \code{meas[i,,]} is a matrix where the m alternatives are the
rows and the n criteria the columns. The values must be standardized
measurements (i.e.  after application of the partial value function).}
  \item{cw}{An \eqn{m \times n}{m * n} matrix of central weights, where
each row corresponds to an alternative and each column to a criterion.}
}
\value{
An object of class \code{smaa.cf}, with the following elements:
  \item{cf}{A vector of confidence factors, one for each alternative.}
  \item{cw}{The central weights (see \code{\link{smaa.cw}}).}
The number of SMAA iterations is stored in \code{attr(x, "smaa.N")}.
}
\details{
The confidence factor for an alternative is its first-rank acceptability
under its central weight.
}
\author{
Gert van Valkenhoef
}
\seealso{
\code{\link{smaa.cw}}
}
\examples{
N <- 1E4; m <- 2; n <- 3
meas <- dget(system.file("extdata/thrombo-meas.txt.gz", package="smaa"))
pref <- dget(system.file("extdata/thrombo-weights-nopref.txt.gz", package="smaa"))

# Calculate central weights
values <- smaa.values(meas, pref)
ranks <- smaa.ranks(values)
cw <- smaa.cw(ranks, pref)
print(cw)
cf <- smaa.cf(meas, cw)
print(cf)
plot(cf)
}
