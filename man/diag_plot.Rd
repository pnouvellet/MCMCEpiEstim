% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diag_plot.R
\name{diag_plot}
\alias{diag_plot}
\title{diagnostic plots}
\usage{
diag_plot(I_NB, logged, max_x = 1, dist, Rt, k)
}
\arguments{
\item{logged}{plot y_axis on natural scale, i.e. FALSE, of log10 scale, logged = TRUE,}

\item{max_x}{plot up to x*max(x), default show the whole axis, i.e. max_x = 1}

\item{dist}{assumed offspring distribution: 'poisson' or 'nb'}

\item{Rt}{data.frame with t = time and vector of Rt's assumed or estimated}

\item{k}{relevant only when dist = 'nb', assumed or estimated overdispersion (single estimate)}
}
\description{
diagnostic plots
}
\details{
return corrected variances to decrease/increase acceptance
}
