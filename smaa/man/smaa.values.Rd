\name{smaa.values}
\alias{smaa.values}
\title{
SMAA alternative values
}
\description{
Calculate the alternative values based on a set of samples from the criteria values distribution and a set of samples from the feasible weight space.
}
\usage{
smaa.values(meas, pref)
}
\arguments{
  \item{meas}{Criteria measurements. An \eqn{N \times m \times n}{N * m * n}
array, where \code{meas[i,,]} is a matrix where the m alternatives are the
rows and the n criteria the columns. The values must be standardized
measurements (i.e.  after application of the partial value function).
\code{\link{smaa.pvf}} provides a convenience method to standardize partial values.}
  \item{pref}{Weights. An \eqn{N \times n}{N * n} array, where \code{pref[i,]}
is a normalized weight vector.}
}
\value{
An \eqn{N \times m}{N * m} array of alternative values.
}
\author{
Gert van Valkenhoef
}
\seealso{
\code{\link{smaa.pvf}}
\code{\link{smaa.ranks}}
}
\examples{
N <- 1E4; m <- 2; n <- 3
meas <- dget(system.file("extdata/thrombo-meas.txt.gz", package="smaa"))

# Sample / read weights
\dontrun{
library(hitandrun)
pref <- simplex.sample(n, N)$samples}
pref <- dget(system.file("extdata/thrombo-weights-nopref.txt.gz", package="smaa"))

# Calculate alternative values
values <- smaa.values(meas, pref)
summary(values)
plot(values)

values.expected <- dget(system.file("extdata/thrombo-values-nopref.txt.gz", package="smaa"))
stopifnot(all.equal(values, values.expected))
}
