\name{smaa.entropy}
\alias{smaa.entropy}
\alias{smaa.entropy.ranking}
\alias{smaa.entropy.choice}
\title{
Decision entropy
}
\description{
Calculate decision entropy from the sampled SMAA rankings. For both ranking and choice problematics.
}
\usage{
smaa.entropy.ranking(ranks, p0 = 1)
smaa.entropy.choice(ra, p0 = 1)
}
\arguments{
  \item{ranks}{Object of class \code{smaa.ranks} containing sampled SMAA rankings.}
  \item{ra}{Object of class \code{smaa.ra} containing SMAA rank acceptabilities.
  Alternatively, an object of class \code{smaa.ranks} from which the rank acceptabilities will be
  calculated.}
  \item{p0}{Baseline probability for the entropy calculation.}
}
\details{
Calculates the entropy for the given problematic, quantifying either the uncertainty in the ranking
of the alternatives (where the outcome space \eqn{Y}{Y} consists of the \eqn{m!}{m!} possible
rankings) or in the choice of the best alternative (where the outcome space \eqn{Y}{Y} consists of
the \eqn{m}{m} alternatives). The entropy is given by:
\deqn{H(Y|W) = -\sum_{y \in Y} p_0 p(y|W) \log p_0 p(y|W)}{H(Y|W) = -sum [over y in Y] [ p0 p(y|W)
log p0 p(y|W) ]}
where \eqn{W}{W} is the space of feasible weights.
Since the SMAA analysis samples from the outcome space, the \eqn{p(y|W)}{p(y|W)} can be estimated
directly from the given sample.
}
\value{
The entropy (a single numeric value).
}
\references{
G. van Valkenhoef and T. Tervonen, \emph{Optimal weight constraint elicitation for additive
multi-attribute utility models}, presentation at EURO 2013, July 2013, Rome, Italy.
}
\author{
Gert van Valkenhoef
}
\note{
The number of samples needed to accurately estimate \eqn{H(Y|W)}{H(Y|W)} for the ranking problematic
is currently unknown.
}
\seealso{
\code{\link{smaa.ranks}}
\code{\link{smaa.ra}}
}
\examples{
N <- 1E4; m <- 2; n <- 3
meas <- dget(system.file("extdata/thrombo-meas.txt.gz", package="smaa"))
pref <- dget(system.file("extdata/thrombo-weights-nopref.txt.gz", package="smaa"))

# Calculate ranks
values <- smaa.values(meas, pref)
ranks <- smaa.ranks(values)

# Calculate ranking entropy
smaa.entropy.ranking(ranks)

# Calculate choice entropy
# (equal to ranking entropy because there are only two alternatives)
smaa.entropy.choice(ranks)
smaa.entropy.choice(smaa.ra(ranks))
}
