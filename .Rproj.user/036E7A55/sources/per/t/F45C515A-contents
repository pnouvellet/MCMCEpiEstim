#' simulate
#'
#' simulate sinusoidal Rt
#'  try to tune variance toward 20% acceptance
#' 
#' @param a matrix, samples of posterior distribution. ncol: nb parameters, nrow: nb samples
#' 
#' @param period vector, proposal variances used to obtain posterior samples theta
#'                   
#' @param m integer, nro of theta
#'
#' @param t_max vector, proposal variances used to obtain posterior samples theta
#'                   
#' @param t_start integer, nro of theta
#'
#' 
#' @details lambda incidence weighted by serial interval
#' 
#' @export
#' 
#' 

# 
sim_Rt_sin <- function(a, period, m, t_max, t_start){

  Rt <- data.frame(t = seq(t_start,t_max))
  Rt$Rt <- a*sin(2*pi*Rt$t/period)+m
  
  return(Rt)
}
