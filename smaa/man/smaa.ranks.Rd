\name{smaa.ranks}
\alias{smaa.ranks}
\title{
SMAA ranking
}
\description{
Calculate SMAA ranks based on the sampled alternatives' values.
}
\usage{
smaa.ranks(values)
}
\arguments{
  \item{values}{An \eqn{N \times m}{N * m} array of sampled alternative values, where \eqn{N}{N} is
the number of SMAA iterations and \eqn{m}{m} is the number of alternatives.}
}
\value{
An \eqn{N \times m}{N * m} array of ranks obtained by each alternative in each iteration.
}
\author{
Gert van Valkenhoef
}
\note{
The value is given class \code{smaa.ranks}, use \code{unclass(x)} to treat it as a regular matrix.
}
\seealso{
\code{\link{smaa.values}}
\code{\link{smaa.ra}}
\code{\link{smaa.cw}}
}
\examples{
N <- 1E4; m <- 2; n <- 3
meas <- dget(system.file("extdata/thrombo-meas.txt.gz", package="smaa"))
pref <- dget(system.file("extdata/thrombo-weights-nopref.txt.gz", package="smaa"))

# Calculate alternative ranks
values <- smaa.values(meas, pref)
ranks <- smaa.ranks(values)
summary(ranks)

ranks.expected <- dget(system.file("extdata/thrombo-ranks-nopref.txt.gz", package="smaa"))
stopifnot(all.equal(ranks, ranks.expected))
}
