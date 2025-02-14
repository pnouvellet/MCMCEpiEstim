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
#' @return Return data frame with 2 columns for time (\code{t}) and the Rt (\code{$Rt})
#' 
#' @export
#' 
#' @examples
#' 
#' Rt <- sim_Rt_sin(a = 3, period = 20, m = 1, t_start = 1, t_max = 100)
#' plot(Rt$t, Rt$Rt, type='l')
#' 
#' 

# 
sim_Rt_sin <- function(a, period, m, t_start, t_max){
  # make dataframe with time
  Rt <- data.frame(t = seq(t_start,t_max))
  # add Rt as a sinusoidal function
  Rt$Rt <- a*sin(2*pi*Rt$t/period)+m
  
  return(Rt)
}
