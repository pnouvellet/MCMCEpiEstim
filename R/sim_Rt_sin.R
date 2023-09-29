#' simulate
#'
#' simulate sinusoidal Rt
#' 
#' @param a amplitude of the wave
#' 
#' @param period period of Rts oscillations
#'                   
#' @param m overall mean of Rt
#'
#' @param t_start starting time for simulation of Rt
#'
#' @param t_max end time
#'                   
#' 
#' @details return dataframe with Rt and time
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
