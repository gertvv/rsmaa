\name{smaa.cw}
\alias{smaa.cw}
\title{
SMAA central weights
}
\description{
Calculate SMAA central weights from sampled rankings and the corresponding weights.
}
\usage{
smaa.cw(ranks, pref)
}
\arguments{
  \item{ranks}{An \eqn{N \times m}{N * m} array of sampled rankings, where \eqn{N}{N} is the number
of SMAA iterations and \eqn{m}{m} is the number of alternatives.}
  \item{pref}{An \eqn{N \times n}{N * n} array of sampled rankings, where \eqn{N}{N} is the number
of SMAA iterations and \eqn{n}{n} is the number of alternatives.}
}
\value{
An \eqn{m \times n}{m * n} matrix of central weights, where each row corresponds to an
alternative and each column to a criterion.
The number of SMAA iterations is stored in \code{attr(x, "smaa.N")}.
}
\note{
The value is given class \code{smaa.cw}, use \code{unclass(x)} to treat it as a regular matrix.
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

# Calculate central weights
values <- smaa.values(meas, pref)
ranks <- smaa.ranks(values)
cw <- smaa.cw(ranks, pref)
print(cw)
plot(cw)
}
