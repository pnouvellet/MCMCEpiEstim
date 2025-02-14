#' adaptive tuning
#'
#' tune proposal and give good initial values to start 'proper' MCMC
#'  try to tune variance toward 20% acceptance
#' 
#' @param repli number of time the variance of the proposal is tuned (10 tends to be ok)
#' 
#' @param within_iter iterations for evaluation of the acceptance with new proposal variances
#'                   
#' @param sigma original variance to start with
#'
#' @param others same as in MCMC_iter function
#'
#' 
#' @details res list of 2 vectors: theta0: parameters value (posterior samples) at the last iterations
#'                       sigma: new variances for the proposal distributions
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


