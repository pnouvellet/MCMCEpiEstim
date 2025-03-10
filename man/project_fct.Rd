% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_fct.R
\name{project_fct}
\alias{project_fct}
\title{Project incidence}
\usage{
project_fct(
  I0,
  Rt,
  n_loc,
  t_max,
  si,
  p,
  model = "poisson",
  over = NULL,
  k.seed = NULL,
  threshold_cumCase = NULL
)
}
\arguments{
\item{I0}{an incidence object, specifying the initial incidence. 
The initial incidence is here assumed the same in every location.}

\item{Rt}{A data.frame with time (\code{$t}) and associated daily Rts (instantaneous Rt, \code{$Rt}).}

\item{n_loc}{number of locations.}

\item{t_max}{number of simulated time steps.}

\item{si}{serial distribution (as in EpiEstim include a 0 weighted SI on same day).}

\item{p}{level of reporting, i.e. probability that a case simulated is actually observed.
p can take a single value and that would assume constant reporting across time and locations simulated.
Otherwise, p can be data.frame similar to Rt (\code{$t} for time and \code{$p} for the reporting). Such p data.frame
would be longer (1 more value) than Rt as it would be specifying in the reporting for imported cases seeding 
the simulation.
In one of our application, we first specify Rt and then then simulate reporting with a single time-step change 
in reporting:

t_change_p <- input$t_pi_change*input$step+1

Pi_t <- data.frame(t = c(Rt$t,tail(Rt$t,1)+1),
                   pi = c(rep(input$pi,t_change_p),
                         rep(input$pi2, nrow(Rt)+1-t_change_p )) )
 input$pi, input$pi2, input$t_pi_change specifying the 2 reporting probability and last time (in weeks) 
 when pi1 applies respectively.}

\item{model}{takes 'poisson' or 'negbin', i.e. offspring distribution assumed. 
By default, 'poisson' is assumed.}

\item{over}{a real number quantifying the level of overdispersion in the offspring
distribution. The overdispersion is defined in the context of a negative binomial distribution 
as: if X follows a negative binomial distribution with mean \mu and overdispersion \delta, then the 
$$ E[X] = \mu $$ and $$ Var[X] = \mu + frac{\mu^2}{\delta} $$.}

\item{k.seed}{an optional argument to set the seed of the simulations}

\item{threshold_cumCase}{an optional argument to remove location with a cumulative number of 
cases (observed) below threshold_cumCase. if some locations are removed, the function keep simulating new ones
until the the right number of simulations (n_loc) have reached the required cumulative number of cases.}
}
\value{
A list containing 2 data.frames of time (\code{$t}) and daily incidence values at each location 
(\code{$simX}, with X an integer for each location). 
The first data.frame, \code{$I_true}, contains the 'true' incidence, while the second, \code{$I_obs}, contains the 'observed'
incidence once under-reporting is accounted for.
}
\description{
Wrapper for projection with/without overdispersion and with under-reporting.
}
\details{
Simulate incidence using the package 'projection', with a time-varying Rt, The possibility
to specify the presence of overdispersion and the presence of under-reporting of cases.
}
\examples{
# a serial interal
si <- c(0,1)
# simualte an Rt
input <- data.frame(R_min = 0.95,
                    R_max = 1.1,
                    k = 1,
                    period = 24, 
                    step = 4,
                    I0 = 100,
                    pi = 0.5)
Rt <- Rt_linear(inp = input, n_week = input$period)
plot(Rt$t, Rt$Rt,type='l')
# set initial incidence
I0 <- incidence::as.incidence(x = input$I0, dates = 1, interval = 1)

# simulate for 10 locations
res <- MCMCEpiEstim::project_fct(I0 = I0,
                                 Rt = Rt,
                                 n_loc = 10,
                                 t_max = nrow(Rt),
                                 si = si,
                                 p = input$pi,
                                 model = 'negbin',
                                 over = input$k)

layout(matrix(1:2,nrow = 2))
matplot((res$I_true[,-1]), xlab = '', ylab = 'cases', main = 'true cases')
matplot((res$I_obs[,-1]), xlab = '', ylab = 'cases', main = 'observed cases')


}
