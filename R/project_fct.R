#' project
#'
#' wrapper for projection
#' 
#' @param theta matrix, samples of posterior distribution. ncol: nb parameters, nrow: nb samples
#' 
#' @param s vector, proposal variances used to obtain posterior samples theta
#'                   
#' @param it integer, nro of theta
#'
#' 
#' @details lambda incidence weighted by serial interval
#' @export
#' 
#' 

project_fct <- function(I0, Rt, n_sim, t_max, si, model = 'poisson',size = NULL){
  
  I <- as.data.frame(matrix(NA,ncol = n_sim+1, nrow = t_max))
  names(I) <- c('t',paste0('sim',1:n_sim))
  I$t <- 1:t_max
  
  # project
  if (model == 'poisson'){
    temp <- as.data.frame(project(x = I0,
                                  R = Rt$Rt, 
                                  si = si[-1], n_sim = n_sim, time_change = 2:(t_max-0-I0$timespan),
                                  n_days = t_max-I0$timespan, 
                                  R_fix_within = TRUE, 
                                  model = 'poisson',instantaneous_R = TRUE))
  }else if (model == 'negbin'){
    temp <- as.data.frame(project(x = I0,
                                  R = Rt$Rt, 
                                  si = si[-1], n_sim = n_sim, time_change = 2:(t_max-0-I0$timespan),
                                  n_days = t_max-I0$timespan, 
                                  R_fix_within = TRUE, 
                                  model = 'negbin',instantaneous_R = TRUE, size = size))
  }
  
  
  I[1:I0$timespan,] <- cbind(I0$dates,matrix(data = I0$counts,
                                             nrow = nrow(I0),
                                             ncol = n_sim,byrow = FALSE))
  
  I[(I0$timespan+1):(nrow(I)),] <- temp
  
  return(I)
  
}