% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmi510.R
\name{is_min}
\alias{is_min}
\title{Minimum Value Finder}
\usage{
is_min(x, na.rm = T)
}
\arguments{
\item{x}{An atomic vector of numbers}

\item{na.rm}{Whether to omit NAs or not (default is true)}
}
\value{
A logical vector with TRUE being the minimum value
}
\description{
Accepts an atomic vector x and returns a logical with TRUE where x equals its minimum value.
}
\examples{
test_vec = c(7, 6, 8, 5)
is_min(test_vec)
}
