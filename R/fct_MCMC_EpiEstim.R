#' wrapper for MCMC EpiEstim
#'
#' run the MCMC to sample posterior of of Rts (and overdispersion) at each location, require
#' running first the 'fct_EpiEstim' (need to clean up input here so no duplications)
#' 
#' 
#' 
#' @param I0 the incidence for the time window during which we assume Rt to be constant.  
#'           I is a dataframe, first column are dates then incidence for all locations
#'           nb of row is the size of time widows, dates must be sequential
#' 
#' @param I integer of  numbers of locations
#'                   
#' @param t_window integer, the number of iteration for the MCMC
#'
#' @param mean_prior mean prior for Rt
#'
#' @param std_prior std deviation prior for Rt
#' 
#' @param res_EpiEstim starting at Rts at estimate from epiestim
#' 
#' @param overdispersion: FALSE/TRUE, TRUE: single overdispersion across locations
#'              
#' @param rep iterations for MCMC, rep for the final sampling, rep/10 10 times for optimisation of proposal
#' 
#' @param thin thinning of posterior sample
#' 
#' @param param_agg iterations for evaluate the accpetance with new proposal variances
#' 
#' @param Rt0_epiEstim iterations for evaluate the accpetance with new proposal variances
#' 
#' @details  res a list containing 2 matrices: theta: matrix of posterior samples
#'                      and logL: matrix of associated log-likelihood
#' @export
#' 

fct_MCMC_EpiEstim <- function(I0, I, t_window,
                              mean_prior, std_prior,
                              res_EpiEstim, overdispersion = FALSE, 
                              rep, thin = 10, param_agg = FALSE, Rt0_epiEstim = TRUE ){
  
  prior <- epitrix::gamma_mucv2shapescale(mu = mean_prior, cv = std_prior/mean_prior)
  t_max <- nrow(I)
  n_loc <- ncol(I)-1
  
  # time windows
  t_start <- seq(I0$timespan+1, t_max-t_window+1,by = 1)        
  t_end <- t_start + t_window - 1     
  n_tw <- length(t_start)
  
  # parameters & proposal log_normal standard deviation
  if (param_agg){
    if(Rt0_epiEstim){
      temp <- apply(matrix(unlist(lapply(res_EpiEstim, "[", ,'Mean(R)')), nrow = t_max, ncol = n_loc, byrow = FALSE),
                    1,mean,na.rm=TRUE)
      Rts_0 <- temp[(t_window+I0$timespan):length(temp)]
    }else{
      Rts_0 <- rep(1, n_tw)
    }
    s_Rt <- rep(0.1, n_tw)
  }else{
    if(Rt0_epiEstim){
      temp <- matrix(unlist(lapply(res_EpiEstim, "[", ,'Mean(R)')), nrow = t_max, ncol = n_loc, byrow = FALSE)
      Rts_0 <- c(temp[(t_window+I0$timespan):nrow(temp),])
    }else{
      Rts_0 <- rep(1, n_loc*n_tw)
    }
    s_Rt <- rep(0.1, n_loc*n_tw)
  }
  if(overdispersion == TRUE){
    Over_0 <- 100
    s_Over <- .1
  }else{
    Over_0 <- NULL
    s_Over <- NULL
  }
  Theta_0 <- list(Rts = Rts_0, Over = Over_0)
  s_0 <- list(Rts = s_Rt, Over = s_Over)
  
  # precompute I and infectivity matrices
  Inc <- as.matrix(I[,-1])
  Oi <- matrix(unlist(lapply(res_EpiEstim, "[", ,'Oi')), nrow = t_max, ncol = n_loc, byrow = FALSE)
  Oi[Oi==0] <- NA
  idx_inc <- c(apply(cbind(t_start,t_end), 1,f1_idx_inc))
  
  if(param_agg){
    data_long <- data.frame(Inc_lk = c(Inc[idx_inc,]),
                            Oi_lk = c(Oi[idx_inc,]),
                            sim = c(matrix(1:n_loc, nrow = n_tw*t_window, ncol = n_loc, byrow = TRUE)),
                            Rt = rep(seq(1,n_tw), each = t_window ) )
  }else{
    data_long <- data.frame(Inc_lk = c(Inc[idx_inc,]),
                            Oi_lk = c(Oi[idx_inc,]),
                            sim = c(matrix(1:n_loc, nrow = n_tw*t_window, ncol = n_loc, byrow = TRUE)),
                            Rt = rep(seq(1,n_tw*n_loc), each = t_window ) )
  }
  # res <- MCMC_iter(iter = rep, theta0 = Theta_0, s = s_0, data_long = data_long,
  #                  n_loc = n_loc, n_tw = n_tw, t_window = t_window,
  #                  prior = prior,overdispersion = overdispersion, param_agg )

  res <- MCMC_full(iter = rep, theta0 = Theta_0, s = s_0,
                   repli_adapt = 10, within_iter = rep/10,
                   data_long = data_long, n_loc = n_loc, n_tw = n_tw, 
                   t_window = t_window, prior = prior, 
                   overdispersion = overdispersion, thin = thin, param_agg )
  
  
  return(res)
}