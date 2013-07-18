\name{smaa.ra}
\alias{smaa.ra}
\title{
SMAA rank acceptabilities
}
\description{
Calculate SMAA rank acceptabilities from sampled rankings.
}
\usage{
smaa.ra(ranks)
}
\arguments{
  \item{ranks}{An \eqn{N \times m}{N * m} array of sampled rankings, where \eqn{N}{N} is the number
of SMAA iterations and \eqn{m}{m} is the number of alternatives.}
}
\value{
An \eqn{m \times m}{m * m} matrix of rank probabilities, where each row corresponds to an
alternative.
The number of SMAA iterations is stored in \code{attr(x, "smaa.N")}.
}
\note{
The value is given class \code{smaa.ra}, use \code{unclass(x)} to treat it as a regular matrix.
}
\author{
Gert van Valkenhoef
}
\seealso{
\code{\link{smaa.ranks}}
}
\examples{
N <- 1E4; m <- 2; n <- 3
meas <- dget(system.file("extdata/thrombo-meas.txt.gz", package="smaa"))
pref <- dget(system.file("extdata/thrombo-weights-nopref.txt.gz", package="smaa"))

# Calculate rank acceptabilities
values <- smaa.values(meas, pref)
ranks <- smaa.ranks(values)
ra <- smaa.ra(ranks)
print(ra)
plot(ra)
}
