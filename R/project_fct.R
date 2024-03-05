#' project
#'
#' wrapper for projection with/without overdispersion 
#' 
#' @param I0 matrix of initial incidence per location (can be more than one day)
#' 
#' @param Rt dataframe with t (time) and associated Rts (instantaneous Rt)
#'                   
#' @param n_loc number of locations
#' 
#' @param t_max number of simulated time steps
#'
#' @param si serial distribution (as in EpiEstim include a 0 weighted SI on same day)
#'
#' @param model take 'poisson' or 'NB', i.e. offspring distribution assumed.
#'
#' @param over value of overdispersion in the offspring distribution (here no under-reporting)
#'
#' 
#' @details I simulated incidence, inclusive of starting values
#' 
#' @export
#' 
#' 

project_fct <- function(I0, Rt, n_loc, t_max, si, p,
                        model = 'poisson',over = NULL){
  
  I <- as.data.frame(matrix(NA,ncol = n_loc+1, nrow = t_max+1))
  names(I) <- c('t',paste0('sim',1:n_loc))
  I$t <- 1:(t_max+1)
  
  # project
  if (model == 'poisson'){
    temp <- as.data.frame(project(x = I0,
                                  R = Rt$Rt, 
                                  si = si[-1], n_sim = n_loc, time_change = 1:(t_max-1-I0$timespan+1),
                                  n_days = t_max-I0$timespan+1, 
                                  R_fix_within = TRUE, 
                                  model = 'poisson',instantaneous_R = TRUE))
  }else if (model == 'negbin'){
    temp <- as.data.frame(project(x = I0,
                                  R = Rt$Rt, 
                                  si = si[-1], n_sim = n_loc, time_change = 1:(t_max-1-I0$timespan+1),
                                  n_days = t_max-I0$timespan+1, 
                                  R_fix_within = TRUE, 
                                  model = 'negbin',instantaneous_R = TRUE, size = over))
  }
  
  
  I[1:I0$timespan,] <- cbind(I0$dates,matrix(data = I0$counts,
                                             nrow = nrow(I0),
                                             ncol = n_loc,byrow = FALSE))
  
  I[(I0$timespan+1):(nrow(I)),] <- temp
  
  I_obs <- I
  
  if(length(p)>1){
    p <- p$pi
  }
  
  temp <- matrix(rbinom(n = (t_max+nrow(I0))*n_loc,
                        size = as.matrix(I[,-1]),
                        prob = matrix(p,nrow = t_max+nrow(I0),ncol = n_loc,byrow = FALSE) ),
                 nrow = t_max+nrow(I0),ncol = n_loc,byrow = FALSE)
  
  I_obs[,-1] <- matrix(rbinom(n = nrow(I)*(ncol(I)-1),
                              size = as.matrix(I[,-1]),prob = p),
                       nrow = nrow(I),ncol = (ncol(I)-1),byrow = FALSE)
  
  return(list(I_true = I, I_obs = I_obs))
  
}
