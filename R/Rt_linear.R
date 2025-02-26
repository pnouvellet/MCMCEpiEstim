#' Rt_linear
#'
#' Function to simulate a periodic and linearly changing Rt. The function allows 
#' stepwise changes of Rt.
#' 
#' @param inp input, need to be a data.frame and include \code{$R_max, $R_min,
#' $period, $step}: the max/min of the reproduction, the period of the 
#' pattern and duration of the step-changes (during 'step' days the Rt is constant).
#' Note: the ratio of period divided
#' by step should be an even integer to ensure that the Rt function is symmetrical,
#' while reaching its maximum and its minimum.
#'
#' @param n_week number of weeks for the simulation, default of 24 weeks.
#' Note: the number of week divide by the period * 7 should should be an integer, to
#' ensure we have complete cycles.
#' 
#' @return Over the period of simulation (\code{n_week}), the function provide
#' Rt values following a periodic pattern with a number of step-wise changes, which
#' depends on the input.
#'  
#' @return A data.frame of time (\code{$t}) and daily Rt values (\code{$Rt}), 
#' size: (\code{n_week} x 7) x 2
#' 
#' @export
#' 
#' @examples
#' # with linear change in Rt (step = 1)
#' inp <- data.frame(R_max=.5, R_min=2,period = 16, step = 1)
#' Rt <- Rt_linear(inp = inp, n_week = inp$period)
#' plot(Rt$t, Rt$Rt, type='l')
#' 
#' # with linear step changes in Rt (step = 4)
#' inp <- data.frame(R_max=.5, R_min=2,period = 16, step = 4)
#' Rt <- Rt_linear(inp = inp, n_week = inp$period)
#' plot(Rt$t, Rt$Rt, type='l')
#' 

Rt_linear <- function(inp,n_week = 24){
  
  # ratio of period by step
  ratio <- inp$period/inp$step
  
  # check satisfy conditions
  if((ratio %% 1)!=0) { stop("ratio of period by step should be integer")}
  if(ratio %% 2) { stop("ratio of period by step should be even")}
  
  # make half of the sequence, increasing phase
  temp <- seq(inp$R_max, inp$R_min,length.out = ratio/2+1)
  # inverse the sequence, for the decreasing phase (without the last step)
  temp <- c(temp,rev(temp[-1])[-1])
  
  
  n_rep_cycle <- n_week/(inp$period/7)
  if((n_rep_cycle %% 1)!=0) { stop("nb cycle repetitions should be integer (with 24 weeks simulation")}
  
  
  Rt <- data.frame(t = 1:(n_week*7), 
                   Rt = rep(rep(temp, each = inp$step), n_rep_cycle))
  
  
  return(Rt)
  
}
