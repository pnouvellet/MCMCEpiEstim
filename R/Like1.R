#' get likelihood  for Poisson or NB offspring distribution
#'
#' get likelihood  
#' internal to the MCMC
#' 
#' @param theta: parameters values   
#'                  
#' @param data_long dataframe of incidence and overall infectivities
#' 
#' @param t_window time windows
#' 
#' @param n_loc number of locations
#' 
#' @param n_tw number of time windows
#' 
#' @param param_agg TRUE or FALSE if Rt/overdispersion estimate are aggregated by location
#'
#' @details  L log likelihood
#' @export
#' 

Like1Poisson <- function(theta, data_long, t_window, n_loc, n_tw, param_agg = FALSE  ){
  
  R <- theta$Rts
  # place new value in datalong to account for time window smoothing
  Rts_lk <- R[data_long$Rt]
  # get mean new number of cases
  lambda <- Rts_lk*data_long$Oi_lk
  # get the log-likelihood (minus constant bits)
  logL_ind <- dpois(x = data_long$Inc_lk, lambda = lambda, log = TRUE)
  # logL_ind <- data_long$Inc_lk *log(lambda) - lambda 
  
  # aggregate within time windows equivalent but faster than: logL_s <- aggregate(logL_ind_s, by = list(data_long$Rt),sum)[,-1]
  logL <- colSums(matrix(logL_ind, nrow = t_window, ncol = n_loc*n_tw, byrow = FALSE), na.rm = TRUE )
  if (param_agg){
    logL <- colSums(matrix(logL, nrow = n_loc, ncol = n_tw, byrow = TRUE), na.rm = TRUE )
  }
  
  return(logL)
}

# # NB likelihood
# # assume delta = overdispersion in offspring distribution:
# # overdispersion at population level: omicron = OverInfectivity_{t} * delta
# # so Var(X) = mu + mu^2/omicro , with mu = R_{t} * OverInfectivity_{t}
# logL <- dnbinom(x = data_long$Inc_lk, mu = Rts_lk_0*data_long$Oi_lk, 
#                 size = data_long$Oi_lk * theta0_over , log = TRUE) 

Like1NBsp <- function(theta, data_long, t_window, n_loc, n_tw, param_agg = FALSE  ){
  
  R <- theta$Rts
  over <- theta$Over # !!!!!!!need to be modified in light of new results!!!!
  # place new value in datalong to account for time window smoothing
  Rts_lk <- R[data_long$Rt]
  # get mean new number of cases
  lambda <- Rts_lk*data_long$Oi_lk
  
  # get the log-likelihood (minus constant bits)
  logL_ind <- dnbinom(x = data_long$Inc_lk, mu = lambda, 
                      size = data_long$Oi_lk * over , log = TRUE)
  
  
  # aggregate within time windows equivalent but faster than: logL_s <- aggregate(logL_ind_s, by = list(data_long$Rt),sum)[,-1]
  logL <- colSums(matrix(logL_ind, nrow = t_window, ncol = n_loc*n_tw, byrow = FALSE), na.rm = TRUE )
  if (param_agg){
    logL <- colSums(matrix(logL, nrow = n_loc, ncol = n_tw, byrow = TRUE), na.rm = TRUE )
  }
  return(logL)
}