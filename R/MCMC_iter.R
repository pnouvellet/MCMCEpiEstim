#' MCMC iterate
#'
#' run the MCMC to sample posterior of of Rts (and overdispersion) at each location
#' 
#' @param iter integer, the number of iteration for the MCMC
#'
#' @param theta0 vector of inital parameters
#'
#' @param s variance of proposal distributions (log-normal) 
#' 
#' @param data_long dataframe of incidence and overall infectivities
#' 
#' @param n_loc number of locations
#' 
#' @param n_tw number of time windows 
#' 
#' @param t_window time windows
#' 
#' @param prior prior for parameters
#' 
#' @param overdispersion TRUE or FALSE if overdispersion is assumed (Poisson vs. NB)
#'              
#' @param param_agg TRUE or FALSE if Rt/overdispersion estimate are aggregated by location
#' 
#' 
#' @details  res a list containing matrices of: Rts posterior samples (before/after thinning), 
#'               overdispersions posterior samples (before/after thinning).
#' @export
#' 

MCMC_iter <- function(iter,theta0,s, data_long, n_loc, n_tw, t_window, prior, overdispersion, param_agg = FALSE ){
  
  # parameters
  n_param <- data.frame(Rt = length(theta0$Rts),
                        Over = length(theta0$Over))
  if(overdispersion){
    Like1 <- Like1NBsp 
  }else{
    Like1 <- Like1Poisson
  }
  #storage matrices
  L <- matrix(NA, nrow = iter, ncol = sum(n_param))
  Rs <- matrix(NA, nrow = iter, ncol = n_param$Rt)
  Overs <- matrix(NA, nrow = iter, ncol = n_param$Over)
  
  logL_0 <- Like1(theta = theta0, data_long = data_long, t_window = t_window, n_loc = n_loc, n_tw = n_tw, param_agg )
  
  # fill first raw of storage
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
    logL_s <- Like1(theta = theta_s, data_long = data_long, t_window = t_window, n_loc = n_loc, n_tw = n_tw, param_agg )
    
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
      if(theta_s$Over>1e3) theta_s$Over <- theta0$Over 
      
      # get the log-likelihood (minus constant bits)
      logL_s <- Like1(theta = theta_s, data_long = data_long, t_window = t_window, n_loc = n_loc, n_tw = n_tw, param_agg )
      
      # get ratio of likelihood corrected for priors and proposal
      r <- exp(sum(logL_s)-sum(logL_0))*theta_s$Over/theta0$Over  
      
      # if accepted replace Rt and likelihood
      if(runif(n = 1, min = 0, max = 1) <= r){
        theta0$Over <- theta_s$Over 
        logL_0 <- logL_s
      }
      
    }
    
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
