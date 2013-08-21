\name{smaa.pvf}
\alias{smaa.pvf}
\title{
Compute piece-wise linear partial value functions
}
\description{
Given a set of reference levels and their values, compute a linearly
interpolated (piece-wise linear) partial value function.
}
\usage{
smaa.pvf(x, cutoffs, values, outOfBounds="error")
}
\arguments{
  \item{x}{Values to compute the PVF for, a numeric vector}
  \item{cutoffs}{Reference levels (ascending order)}
  \item{values}{Values of the reference levels}
  \item{outOfBounds}{What to do when some of the \code{x} are outside
the range of the given \code{cutoffs}. When \code{"error"}, throws an
error. When \code{"clip"}, clips to the value of the first or last
cutoff. When \code{"interpolate"}, interpolates according to the closest
interval.}
}
\value{
A numeric vector the same length as \code{x}.

The values are computed by linear interpolation between the values of
the two closest reference levels. This has been implemented in C for a
dramatic performance improvement.
}
\author{
Gert van Valkenhoef
}
\examples{
x <- c(50, 90, 100, 10, 40, 101, 120)
values <- smaa.pvf(x,
  cutoffs=c(50, 75, 90, 100),
  values=c(1, 0.8, 0.5, 0),
  outOfBounds="clip")
stopifnot(all.equal(values, c(1.0, 0.5, 0.0, 1.0, 1.0, 0.0, 0.0)))
}
