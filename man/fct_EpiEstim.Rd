% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_EpiEstim.R
\name{fct_EpiEstim}
\alias{fct_EpiEstim}
\title{epiestim wrap}
\usage{
fct_EpiEstim(I0_t_import, I, t_window, mean_prior, std_prior, si, overlap)
}
\arguments{
\item{I0_t_import}{which of the initial incidence is imported}

\item{I}{dataframe of incidence. ncol: nb of location +1 (time), nrow: time}

\item{t_window}{integer, time window}

\item{mean_prior}{single real number, mean prior for Rts}

\item{std_prior}{single real number, stad deviation for prior for Rts}

\item{si}{serial distribution (as in EpiEstim include a 0 weighted SI on same day)}

\item{overlap}{TRUE/FALSE, whether using overlapping time window or not}
}
\description{
wrapper to call EpiEstim independently on multiple locations
}
\details{
res, list of 1 dataframe per location summarising incidence, time, overall infectivity,
             time window start and end date, estimated mean/std deviation/median and
              95% CI lower and upper quantile for Rts
}
