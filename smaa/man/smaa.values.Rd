\name{smaa.values}
\alias{smaa.values}
\title{
SMAA alternative values
}
\description{
Calculate the alternative values based on a set of samples from the criteria values discribution and a set of samples from the feasible weight space.
}
\usage{
smaa.values(meas, pref)
}
\arguments{
  \item{meas}{Criteria measurements. An \eqn{N \times m \times n}{N * m * n}
array, where \code{meas[i,,]} is a matrix where the m alternatives are the
rows and the n criteria the columns. The values must be standardized
measurements (i.e.  after application of the partial value function).}
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
\code{\link{smaa.ranks}}
}
\examples{
}
