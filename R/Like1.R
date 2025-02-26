#' Compute the likelihood  for Poisson or NB offspring distributions
#'
#' 
#' This is an internal function.
#' it compute the likelihood for given parameter set and the data.
#' Run the MCMC to sample the posterior of Rts (and optionally overdispersion) at each location.
#' this function is called internally by  MCMC_iter.R, adapt_tuning.R, and MCMC_Full.R, which is itself called in the 
#' fct_MCMC_EpiEstim.R.
#' 
#' @param theta0 list of 2 vectors (for Rt and overdispersion) of inital values for parameters (set in MCMC_iter.R)
#'                  
#' @param data_long data.frame of incidence and overall infectivities by locations in long format 
#' (i.e. see  fct_MCMC_EpiEstim.R)
#' 
#' @param t_window integer, time window during which Rt is assumed constant.#' 
#' 
#' @param n_loc number of locations
#' 
#' @param n_tw number of time windows (set in fct_MCMC_EpiEstim.R)
#'  
#' @param param_agg TRUE or FALSE if Rts estimates are aggregated by location
#'
#' @param overdispersion TRUE or FALSE if overdispersion is assumed (Poisson vs. NB)
#' 
#' @param p_reps reporting probability, either a single real (if constant) or a vector 
#' with a value of reporting for each day. (set in fct_MCMC_EpiEstim.R)
#'
#' @return  logL log likelihood, a row vector of loglikelihoods. when param_agg is FALSE, 
#' logL has size number of independent time windows times number of location,
#' when param_agg is FALSE, logL has size number of independent time windows.
#' 
#' 
#' @export
#' 


# # NB likelihood
# # assume delta = overdispersion in offspring distribution:
# # overdispersion at population level: omicron = OverInfectivity_{t} * delta
# # so Var(X) = mu + mu^2/omicro , with mu = R_{t} * OverInfectivity_{t}
# logL <- dnbinom(x = data_long$Inc_lk, mu = Rts_lk_0*data_long$Oi_lk, 
#                 size = data_long$Oi_lk * theta0_over , log = TRUE) 

Like1 <- function(theta, data_long, t_window, n_loc, n_tw, param_agg = FALSE, overdispersion, p_reps  ){
  
  R <- theta$Rts
 
  # place new value in datalong to account for time window smoothing
  Rts_lk <- R[data_long$Rt]
  # get mean new number of cases
  lambda <- Rts_lk*data_long$Oi_lk
  
  if(overdispersion){
      # over <- theta$Over # 
      # varnb <- lambda*(1+Rts_lk/theta$Over).  # no under-reporting
      varnb <- lambda*(1+(1-p_reps)* Rts_lk +p_reps*Rts_lk/theta$Over) + 
        (1-p_reps)* Rts_lk*(1+ Rts_lk +p_reps*Rts_lk/theta$Over)
    # get the log-likelihood (minus constant bits)
    logL_ind <- dnbinom(x = data_long$Inc_lk, mu = lambda, 
                        size = lambda^2/(varnb-lambda) , log = TRUE)
  }else{
    if(mean(p_reps)==1){
      logL_ind <- dpois(x = data_long$Inc_lk, lambda = lambda, log = TRUE)
    }else{
      varnb <- lambda*(1+(1-p_reps)* Rts_lk ) + 
        (1-p_reps)* Rts_lk*(1+ Rts_lk )
      # get the log-likelihood (minus constant bits)
      logL_ind <- dnbinom(x = data_long$Inc_lk, mu = lambda, 
                          size = lambda^2/(varnb-lambda) , log = TRUE)
    }
  }
  
  
  # aggregate within time windows equivalent but faster than: logL_s <- aggregate(logL_ind_s, by = list(data_long$Rt),sum)[,-1]
  logL <- colSums(matrix(logL_ind, nrow = t_window, ncol = n_loc*n_tw, byrow = FALSE), na.rm = TRUE )
  if (param_agg){
    logL <- colSums(matrix(logL, nrow = n_loc, ncol = n_tw, byrow = TRUE), na.rm = TRUE )
  }
  return(logL)
}