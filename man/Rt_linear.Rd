% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rt_linear.R
\name{Rt_linear}
\alias{Rt_linear}
\title{Rt_linear}
\usage{
Rt_linear(inp, n_week = 24)
}
\arguments{
\item{inp}{input, need to be a data.frame and include $R_max, $R_min,
$period, $step, the max/min of the reproduction and the period of the 
pattern and duration of the step-changes}

\item{n_week}{nb of weeks for the simulation, default of 24}
}
\value{
A data.frame of time ($t) and daily Rt values ($Rt), size: (n_week x 7) x 2
}
\description{
To specify stewise changes in Rt
}
\details{
Over the period of simulation (n_week), the function provide
Rt values following a periodic pattern with a number of step-wise changes, which
ultimately depend on the input.
}
