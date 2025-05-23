% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Heatmap_Cohen.R
\name{cohen_d}
\alias{cohen_d}
\title{Compute Cohen\'s d Effect Size}
\usage{
cohen_d(x, y)
}
\arguments{
\item{x}{A numeric vector representing the values for group 1.}

\item{y}{A numeric vector representing the values for group 2.}
}
\value{
A numeric value representing Cohen\'s d. Returns NA if either group has fewer than two observations
or if the pooled standard deviation is zero.
}
\description{
Computes the absolute Cohen\'s d effect size between two numeric vectors. This function returns
the absolute value of the difference in means divided by the pooled standard deviation.
}
\keyword{internal}
