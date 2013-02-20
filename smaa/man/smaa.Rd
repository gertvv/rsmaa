\name{smaa}
\alias{smaa}
\title{Core SMAA function}
\description{
Calculate SMAA decision indices based on a set of samples from the criteria values discribution and a set of samples from the feasible weight space.
}
\usage{
    smaa(meas, pref, m=dim(meas)[2], n=dim(meas)[3], N=dim(meas)[1], generate.values=FALSE)
}
\arguments{
    \item{meas}{Criteria measurements. An N * m * n array, where meas[i,,] is
    a matrix where the m alternatives are the rows and the n criteria the
    columns. The values must be standardized measurements (i.e. after
    application of the partial value function).}
    \item{pref}{Weights. An N * n array, where pref[i,] is a normalized
    weight vector.}
    \item{m}{Number of alternatives.}
    \item{n}{Number of criteria.}
    \item{N}{Number of iterations.}
	\item{generate.values}{If TRUE, the generated values / utilities for each of the alternatives are returned.}
}
\details{
TODO: some details here.
}
\value{
    \item{rankAcc}{Rank acceptabilities.}
    \item{centralWeigths}{Central weights.}
	\item{values}{Alternative values / utilities (if generate.values=TRUE)}
}
\examples{

}
\keyword{SMAA}
