#' MCMC iterate
#'
#' run the MCMC to sample posterior of R and initial coniditions at each location 
#' FYI: this is called internally by adapt_tuning
#' 
#' @param incidence the incidence for the time window during which we assume Rt to be constant.  
#'           I is a dataframe, first column are dates then incidence for all locations
#'           nb of row is the size of time widows, dates must be sequential
#' 
#' @param N_geo integer of  numbers of locations
#'                   
#' @param iter integer, the number of iteration for the MCMC
#'
#' @param theta0 vector of inital parameters, here taken from the last MCMC iteration after tuning (save some burn-in)
#'
#' @param s variance of proposal distributions (log-normal) - tuned previously
#' 
#' @param SI Serial interval distribution (see SI_gamma_dist_EpiEstim)
#' 
#' @param mu0: initial conidtions to guaranty that if R=1, then we predict the number of cases in the future will stablise at the mean number of cases observed in the time window
#'              mu0 is also used as the mean of the (exponential) prior for intial conditions estimated
#' 
#' @details  res a list containing 2 matrices: theta: matrix of posterior samples
#'                      and logL: matrix of associated log-likelihood
#' @export
#' 

MCMC_iter <- function(iter,theta0,s, data_long, n_sim, n_tw, t_window, prior, overdispersion, param_agg = FALSE ){
  
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
  
  logL_0 <- Like1(theta = theta0, data_long = data_long, t_window = t_window, n_sim = n_sim, n_tw = n_tw, param_agg )
  
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
    logL_s <- Like1(theta = theta_s, data_long = data_long, t_window = t_window, n_sim = n_sim, n_tw = n_tw, param_agg )
    
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
      logL_s <- Like1(theta = theta_s, data_long = data_long, t_window = t_window, n_sim = n_sim, n_tw = n_tw, param_agg )
      
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
