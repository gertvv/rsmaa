\name{smaa}
\alias{smaa}
\title{One-stage SMAA analysis}
\description{
Calculate SMAA decision indices based on a set of samples from the criteria values distribution and a set of samples from the feasible weight space.
}
\usage{
    smaa(meas, pref)
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
\details{
The one-stage method does not store the alternatives' values or the raw
rankings. Instead, only standard summary metrics are provided.
}
\value{
  \item{ra}{Rank acceptabilities (see \code{\link{smaa.ra}}).}
  \item{cw}{Central weights (see \code{\link{smaa.cw}}).}
}
\author{
Gert van Valkenhoef
}
\seealso{
\code{\link{smaa.pvf}}
}
\examples{
N <- 1E4; m <- 2; n <- 3
meas <- dget(system.file("extdata/thrombo-meas.txt.gz", package="smaa"))

# Sample / read weights
\dontrun{
library(hitandrun)
pref <- simplex.sample(n, N)$samples}
pref <- dget(system.file("extdata/thrombo-weights-nopref.txt.gz", package="smaa"))

# Calculate SMAA metrics
result <- smaa(meas, pref)
print(result)
plot(result)

result <- smaa(meas, c(0.5, 0.2, 0.3))
print(result)
}
\keyword{SMAA}
