\name{smaa-package}
\alias{smaa-package}
\docType{package} 

\title{SMAA: Stochastic Multi-criteria Acceptability Analysis}
\description{
This R package implements the Stochastic Multi-criteria Acceptability Analysis (SMAA) family of
methods for stochastic Multiple Criteria Decision Analysis (MCDA). In SMAA methods, uncertainty in
criteria measurements and incomplete information on the weights are accounted for by Monte Carlo
integration of probability distributions.
}
\details{
The \code{\link{smaa}} function implements the core Monte Carlo integration method. It calculates
the SMAA decision metrics (rank acceptabilities and central weights) in one go.

Alternatively, the \code{\link{smaa.values}}, \code{\link{smaa.ranks}}, \code{\link{smaa.ra}}, and
\code{\link{smaa.cw}} perform the individual steps. Note that \code{\link{smaa}} is slightly more
efficient because it does not store the alternatives' values or rankings.

The \link[hitandrun]{hitandrun-package} is complementary to this package in that it provides methods
for sampling weights when incomplete preference information is available in the form of linear
constraints on the weight vector.
}
\author{Gert van Valkenhoef}
\references{
T. Tervonen and J.R. Figueira (2008),
\emph{A survey on stochastic multicriteria acceptability analysis methods},
Journal of Multi-Criteria Decision Analysis 15(1-2):1-14.
[\href{http://dx.doi.org/10.1002/mcda.407}{doi:10.1002/mcda.407}]

T. Tervonen, G. van Valkenhoef, N. Basturk, and D. Postmus (2012),
\emph{Hit-And-Run enables efficient weight generation for simulation-based multiple criteria decision analysis},
European Journal of Operational Research 224(3):552-559.
[\href{http://dx.doi.org/10.1016/j.ejor.2012.08.026}{doi:10.1016/j.ejor.2012.08.026}]
}
\examples{
N <- 1E4; m <- 2; n <- 3
meas <- dget(system.file("extdata/thrombo-meas.txt.gz", package="smaa"))
pref <- dget(system.file("extdata/thrombo-weights-nopref.txt.gz", package="smaa"))

# Calculate SMAA metrics (one-stage)
result <- smaa(meas, pref)
print(result)

# Calculate SMAA metrics (multi-stage)
values <- smaa.values(meas, pref)
summary(values)
ranks <- smaa.ranks(values)
smaa.ra(ranks)
smaa.entropy.ranking(ranks)
smaa.cw(ranks, pref)

# Calculate confidence factors
smaa.cf(meas, result$cw)
}
