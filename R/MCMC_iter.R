#' Iterate the MCMC inference
#'
#' This is an internal function.
#' Run the MCMC to sample the posterior of Rts (and optionally overdispersion) at each location.
#' this function is called internally by MCMC_Full.R, which is itself called in the 
#' fct_MCMC_EpiEstim.R.
#' 
#' @param iter integer, the number of iterations for the MCMC
#'
#' @param theta0 list of 2 vectors (for Rt and overdispersion) of inital values for parameters (set in fct_MCMC_EpiEstim.R)
#'
#' @param s list of 2 vectors (for Rt and overdispersion) of variances of proposal distributions (log-normal). (set in fct_MCMC_EpiEstim.R) 
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
#' @return  res a list containing matrices of: Rts posterior samples (ncol = number of Rts estimated, nrow = number of iteration), 
#'               overdispersions posterior samples (if estimated), and the associated log-likelihoods.
#' @export
#' 

MCMC_iter <- function(iter,theta0,s, data_long, n_loc, n_tw, t_window, prior, 
                      overdispersion, param_agg = FALSE, p_reps, mean_k_prior, k_upper_limit){
  
  # parameters
  n_param <- data.frame(Rt = length(theta0$Rts),
                        Over = length(theta0$Over))
 
  #storage matrices
  L <- matrix(NA, nrow = iter, ncol = sum(n_param))
  Rs <- matrix(NA, nrow = iter, ncol = n_param$Rt)
  Overs <- matrix(NA, nrow = iter, ncol = n_param$Over)
  
  logL_0 <- Like1(theta = theta0, data_long = data_long, t_window = t_window,
                  n_loc = n_loc, n_tw = n_tw, param_agg, overdispersion = overdispersion, p_reps )
  
  # fill first row of storage
  Rs[1,] <- theta0$Rts
  L[1,1:(n_param[[1]])] <- logL_0
  if(overdispersion){
    L[1,n_param[[1]]+1] <- sum(logL_0)
    Overs[1] <- theta0$Over
  }

  for (i in 2:iter){
    # propose new Rt
    theta_s <- theta0
    theta_s$Rts <- theta_s$Rts*exp(s$Rts*rnorm(n = n_param[[1]], mean = 0, sd = 1) )
    
    # get the log-likelihood (minus constant bits)
    logL_s <- Like1(theta = theta_s, data_long = data_long, t_window = t_window,
                    n_loc = n_loc, n_tw = n_tw, param_agg, overdispersion = overdispersion, p_reps  )
    
    # correct log-likelihood for gamma prior Rt
    corr_prior <- ( (prior$shape-1)*(log(theta_s$Rts) - log(theta0$Rts)) + (theta0$Rts - theta_s$Rts)/prior$scale )
    # get ratio of likelihood corrected for priors and proposal
    r <- exp(logL_s-logL_0)*theta_s$Rts/theta0$Rts * exp(corr_prior)  
    
    # find which Rt-associated time window to accept
    f <- which(runif(n = n_param[[1]], min = 0, max = 1) <= r)
    
    # if accepted replace Rt and likelihood
    if(length(f)>0){
      theta0$Rts[f] <- theta_s$Rts[f] 
      logL_0[f] <- logL_s[f]
    }
    
    # Overdispersion
    if(overdispersion){
      # propose new overdispersion
      theta_s <- theta0
      theta_s$Over <- theta_s$Over*exp(s$Over*rnorm(n = 1, mean = 0, sd = 1) )
      if(k_upper_limit){
        if(theta_s$Over>1e3) theta_s$Over <- theta0$Over
      }
      
      # get the log-likelihood (minus constant bits)
      logL_s <- Like1(theta = theta_s, data_long = data_long, 
                      t_window = t_window, n_loc = n_loc, n_tw = n_tw, param_agg, overdispersion = overdispersion, p_reps )
      
      # correct log-likelihood for prior of overdisp. (1/k~exp(mu_v)) assume mu_v=1/100: Poisson-like
      # corr_prior <- 2*log(theta0$Over/theta_s$Over) - 1 *(1/theta_s$Over - 1/theta0$Over)
      # exponential prior
      corr_prior <- (theta0$Over-theta_s$Over)/mean_k_prior
      
      
      # corr_prior <- 1
      
      # get ratio of likelihood corrected for priors and proposal
      r <- exp(sum(logL_s) - sum(logL_0) + corr_prior)*theta_s$Over/theta0$Over 
      
      # if accepted replace Rt and likelihood
      if(runif(n = 1, min = 0, max = 1) <= r){
        theta0$Over <- theta_s$Over 
        logL_0 <- logL_s
      }
      
    }
    
    # save the accepted new parameters estimates
    Rs[i,] <- theta0$Rts
    L[i,1:(n_param[[1]])] <- logL_0
    if(overdispersion){
      L[i,n_param[[1]]+1] <- sum(logL_0)
      Overs[i] <- theta0$Over
    }
    
  }
  
  ####
  
  if(overdispersion){
    res <- list(theta_R = Rs, theta_over = Overs, logL = L)
  }else{
    res <- list(theta_R = Rs, logL = L)
  }
  return(res)
  
}
