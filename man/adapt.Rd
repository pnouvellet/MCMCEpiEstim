% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adapt.R
\name{adapt}
\alias{adapt}
\title{Tuning}
\usage{
adapt(theta, s, it)
}
\arguments{
\item{theta}{matrix of posterior samples of posterior distribution. ncol: nb parameters, nrow: nb samples}

\item{s}{vector, proposal variances used to obtain posterior samples theta (size the same as ncol of theta)}

\item{it}{integer, number of theta sample}
}
\value{
a vector of corrected variances to decrease/increase acceptance
}
\description{
Tune the variance proposal.
Try to tune variance toward 20% acceptance.
}
