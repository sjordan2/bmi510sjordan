% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmi510.R
\name{accuracy}
\alias{accuracy}
\title{Accuracy of a Model}
\usage{
accuracy(pred, truth)
}
\arguments{
\item{pred}{The predicted values from the model}

\item{truth}{The ground truth values}
}
\description{
Calculates the accuracy of a binary predictor using ground truth and prediction vectors
}
\examples{
ground_truth = c(1, 1, 0, 1, 0, 1, 0, 1, 1)
predicted_vals = c(1, 0, 1, 1, 0, 1, 0, 1, 0)
accuracy(predicted_vals, ground_truth)
}
