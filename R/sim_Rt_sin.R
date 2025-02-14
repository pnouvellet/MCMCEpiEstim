#' Simulate Rt
#'
#' simulate a sinusoidal Rt over time
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
#' @return return data frame with 2 columns for time (\code{t}) and the Rt (\code{$Rt})
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
