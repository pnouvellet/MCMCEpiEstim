#' full MCMC procedure
#'
#' Run the MCMC to sample the posterior of Rts (and optionally overdispersion) at each location.
#' This function is called internally by the function fct_MCMC_EpiEstim.R.
#' This function automatically tune the proposal variances targeting a 20% aceeptance rate
#' The function also run some diagnotics (covergence checks) and estimate the DIC of model.
#' 
#' 
#' @param iter integer, the number of iterations for the MCMC
#'
#' @param theta0 list of 2 vectors (for Rt and overdispersion) of inital values for parameters (set in fct_MCMC_EpiEstim.R)
#'
#' @param s list of 2 vectors (for Rt and overdispersion) of variances of proposal distributions (log-normal). (set in fct_MCMC_EpiEstim.R) 
#'  
#' @param repli_adapt number of time the variance of the proposal is tuned (10 tends to be ok)
#' 
#' @param within_iter iterations for evaluation of the acceptance rate with new proposal variances, when tuning
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
#' @param thin thinning of posterior samples
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
#' @return  res a list containing matrices of: 
#'               Rts (\code{$theta_R}) posterior samples (ncol = number of Rts estimated, nrow = number of iteration), 
#'               overdispersions (\code{$theta_over}) posterior samples (if estimated), 
#'               and the associated log-likelihoods (\code{$logL}).
#'               The thinned versions of the above (\code{$theta_R_thinned}, \code{$theta_over_thinned}, \code{$logL_thinned})
#'               DIC for the model (\code{$DIC}), with the first number being the DIC and the second being the effective number of parameters.
#'               The Gelman and Rubin's convergence diagnostic, from the gelman.diag function (\code{$GRD}).
#'               The effective sample size for estimating the mean, from the effectiveSize function (\code{$ESS}).
#'               
#'               The function evaluate convergence , based on
#'               the upper 95%CI of potential scale reduction factors >1.1 (\code{$GRD_converged} being TRUE or FALSE).
#'               the function outputs a warning if evidence of convergence is not reached
#' 
#' @importFrom coda gelman.diag 
#' @importFrom coda as.mcmc.list 
#' @importFrom coda as.mcmc 
#' @importFrom coda effectiveSize
#' 
#' @export
#' 

MCMC_full <- function(iter, theta0, s, repli_adapt, within_iter, data_long,
                      n_loc, n_tw, t_window, prior, overdispersion, thin, 
                      param_agg = FALSE, p_reps, mean_k_prior, k_upper_limit = TRUE,
                      run_diagnostics){
  
  rep <- repli_adapt*within_iter
  # # initialise likelihood
  # if(overdispersion){
  #   Like1 <- Like1NBsp 
  # }else{
  #   Like1 <- Like1Poisson
  # }
  # 
  res0 <- adapt_tuning(repli = repli_adapt,
                       within_iter = within_iter,
                       theta0 = theta0,
                       sigma = s,
                       data_long = data_long, n_loc = n_loc, n_tw = n_tw, 
                       t_window = t_window, prior = prior, overdispersion = overdispersion,
                       param_agg, p_reps, mean_k_prior, k_upper_limit)
  # adaptive tuning bit: we run an MCMC with rep/10 iterations, then
  # adjust the proposal variance to reach 0.2
  # do again using parameter value from the last iteration of the previous MCMC
  # repeat 10 times
  # from experience, this is enough to tunes proposal variances well, but worth checking
  # see below for final acceptance rate output 
  # see Rscript/MCMC_Rt_2018.R for full function
  
  print('halfway!')             # message halfway through (effectively, including tuning, we do 2xrep iterations)
  
  # res <- MCMC_iter(iter = iter,
  #                  theta0 = res0$theta0,
  #                  s = res0$sigma)
  # browser()
  res <- MCMC_iter(iter = rep, 
                   theta0 = res0$theta0, 
                   s = res0$sigma, 
                   data_long = data_long, n_loc = n_loc, n_tw = n_tw, 
                   t_window = t_window, prior = prior, 
                   overdispersion = overdispersion, 
                   param_agg, p_reps, mean_k_prior, k_upper_limit)
  
  # thinning
  res$theta_R_thinned <- res$theta_R[seq(1, rep, by = thin),]
  res$logL_thinned <- res$logL[seq(1, rep, by = thin),]
  if(overdispersion){
    res$theta_over_thinned <- res$theta_over[seq(1, rep, by = thin),]
  }
  
  # run the MCMC to sample posterior of R and initial coniditions at each location
  # FYI: this is called internally by adapt_tuning
  # see Rscript/MCMC_Rt_2018.R for full function
  # needs:
  # I: the incidence for the time window during which we assume Rt to be constant. 
  # N_geo: the number of locations
  # iter: the number of iteration for the MCMC
  # theta0: inital parameters, here taken from the last MCMC iteration after tuning (save some burn-in)
  # s: variance of proposal distribution (log-normal)
  # SI: the serial interval(use SI_gamma_dist_EpiEstim to define), or need to be a list with vector dist for the daily distribution and SItrunc: an integer for the threshold of serial interval, if SItrunc=40, then dist is 41 element long to include day 0
  # mu0: initial conidtions to guaranty that if R=1, then we predict the number of cases in the future will stablise at the mean number of cases observed in the time window
  # mu0 is also used as the mean of the (exponential) prior for intial conditions estimated
  #Calculate L
  
  # DIC
  # parameters
  n_param <- data.frame(Rt = length(theta0$Rts),
                        Over = length(theta0$Over))
  
  theta_hat <- list(Rts = apply(res$theta_R_thinned,2,median), Over = NULL)
  if (overdispersion){
    theta_hat$Over <- median(res$theta_over_thinned)
  }
  L <- sum(Like1(theta = theta_hat, data_long = data_long, t_window = t_window, 
                 n_loc = n_loc, n_tw = n_tw, param_agg, overdispersion = overdispersion, p_reps ))
  
  ll_med = median(rowSums(res$logL_thinned[,1:(n_param[[1]])]))
  P = 2 * (L - ll_med)
  
  #Calculate DIC
  res$DIC = c(-2 * (L - P),P)
  
  # to check convergence parameter by parameter 
  pars <- names(res)[grep("theta_", names(res))] 
  thinned <- names(res)[grep("thinned", names(res))] 
  non_thinned_pars <- pars[!(pars %in% thinned)]
  
  
  if(run_diagnostics){
    ## get all parameter chains into a single matrix to feed into Gelman Rubin
    tmp_par <- res[[non_thinned_pars[1]]]
    if(length(non_thinned_pars) > 1) {
      for(k in seq(2, length(non_thinned_pars))) {
        tmp_par <- cbind(tmp_par, res[[non_thinned_pars[k]]])  
      }
    }
    spl1 <- tmp_par[seq_len(floor(nrow(tmp_par) / 2)), ]
    spl2 <- tmp_par[seq(ceiling(nrow(tmp_par) / 2) + 1, nrow(tmp_par)), ]
    res$GRD <- gelman.diag(as.mcmc.list(list(as.mcmc(spl1), as.mcmc(spl2))))
    res$ESS <- effectiveSize(tmp_par)
    
    ess_names <- unlist(sapply(non_thinned_pars, function(i)
    {
      if(ncol(res[[i]]) > 1) {
        ret <- paste(i, seq_len(ncol(res[[i]])), sep = "_")
      }else{
        ret <- i
      }
      ret
    }))
    
    names(res$ESS) <- ess_names
    
    # Is any of the potential scale reduction factors >1.1
    # (looking at the upper CI)?
    # If so this would suggest that the MCMC has not converged well.
    if (any(res$GRD$psrf[, "Upper C.I."] > 1.1)) {
      warning("The Gelman-Rubin algorithm suggests the MCMC may not have \n
      converged within the number of iterations (MCMC.burnin + n1) specified. \n
      You could visualise the MCMC chain & decide whether to rerun for longer.\n")
      res$GRD_converged <- FALSE
    } else {
      cat("\nGelman-Rubin MCMC convergence diagnostic was successful.")
      res$GRD_converged <- TRUE
    }
  }
  return(res)
}

