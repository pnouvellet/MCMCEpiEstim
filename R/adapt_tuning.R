#' Adaptive tuning
#'
#' This is an internal function.
#' Run the MCMC to sample the posterior of Rts (and optionally overdispersion) at each location;
#' evaluate and tune the proposal variances toward achieving 20% acceptance.
#' This function is called internally by MCMC_Full.R, which is itself called in the 
#' fct_MCMC_EpiEstim.R.
#' 
#' 
#' @param repli_adapt number of time the variance of the proposal is tuned (10 tends to be ok)
#' 
#' @param within_iter iterations for evaluation of the acceptance rate with new proposal variances, when tuning
#' 
#' @param theta0 list of 2 vectors (for Rt and overdispersion) of inital values for parameters (set in fct_MCMC_EpiEstim.R)
#' 
#' @param sigma original variance to start with. List of 2 vectors (for Rt and overdispersion) 
#' of variances of proposal distributions (log-normal). set by MCMC_full.R, inherited format from fct_MCMC_EpiEstim.R
#'
#' @param data_long data.frame of incidence and overall infectivities by locations in long format 
#' (i.e. see  fct_MCMC_EpiEstim.R)
#' 
#' @param n_loc number of locations 
#' 
#' @param n_tw number of time windows (set in fct_MCMC_EpiEstim.R)
#' 
#' @param t_window integer, time window during which Rt is assumed constant.
#' 
#' @param prior prior for Rt parameter, assume gamma distributed Rt prior with \code{$shape} and \code{$scale}
#' 
#' @param overdispersion TRUE or FALSE if overdispersion is assumed (Poisson vs. NB)
#'              
#' @param param_agg TRUE or FALSE if Rts estimates are aggregated by location
#' 
#' @param p_reps reporting probability, either a single real (if constant) or a vector 
#' with a value of reporting for each day. (set in fct_MCMC_EpiEstim.R)
#' 
#' @param mean_k_prior real, assuming k prior as an exponential distribution, mean_k_prior is the mean of the prior distribution
#' 
#' @param k_upper_limit TRUE or FALSE if k estimates should be bounded to 1,000 (upper limit for k)
#' 
#' 
#' @return res list of 2 lists: theta0: parameter values (posterior samples) at the last iterations;
#'                       sigma: new variances for the proposal distributions.
#' @export
#' 


## adaptative tuning
# tune proposal and give good initial values to start 'proper' MCMC
# repli: number of time the variance of the proposal is tuned (10 tend to be ok)
# within_iter: iteration for evaluate the accpetance with new proposal variances
# sigma: original variance to start with
# others: as in MCMC_iter
adapt_tuning <- function(repli, within_iter, theta0, sigma, 
                         data_long, n_loc, n_tw, t_window, prior, 
                         overdispersion, param_agg = FALSE , p_reps,
                         mean_k_prior, k_upper_limit ){
  
  new_sigma <- sigma
  for (i in 1:repli){
    # run MCMC with wthin_iter iterations
    # res <- MCMC_iter(iter = within_iter,
    #                  theta0 = theta0,
    #                  s = new_sigma)
    res <- MCMC_iter(iter = within_iter, theta0 = theta0, s = new_sigma, 
                     data_long = data_long, n_loc = n_loc, n_tw = n_tw,
                     t_window = t_window, prior = prior, overdispersion = overdispersion, 
                     param_agg, p_reps, mean_k_prior, k_upper_limit )
    
    # colSums(diff(res$theta)!=0)/(within_iter-1) 
    # tune the variance according to accpetance
    new_sigma$Rts <- adapt(theta = res$theta_R,
                           s = new_sigma$Rts,
                           it = within_iter)
    # new starting value for the parameters for the next round of MCMC
    theta0$Rts <- res$theta_R[within_iter,]
    
    # overdispersion
    if(overdispersion){
      new_sigma$Over <- adapt(theta = res$theta_over,
                             s = new_sigma$Over,
                             it = within_iter)
      # new starting value for the parameters for the next round of MCMC
      theta0$Over <- res$theta_over[within_iter]
    }
    

    # plot(res$theta[,2],ylim=c(0,50))
  }
  res <- list(theta0 = theta0,
              sigma = new_sigma)
}


