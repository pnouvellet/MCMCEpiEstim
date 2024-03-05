#' Rt_linear
#'
#' To specify stewise changes in Rt
#' 
#' @param inp input, need to be a data.frame and include $R_max, $R_min,
#' $period, $step, the max/min of the reproduction and the period of the 
#' pattern and duration of the step-changes
#'
#' @param n_week nb of weeks for the simulation, default of 24
#' 
#' @details Over the period of simulation (n_week), the function provide
#' Rt values following a periodic pattern with a number of step-wise changes, which
#' ultimately depend on the input.
#'  
#' @return A data.frame of time ($t) and daily Rt values ($Rt), size: (n_week x 7) x 2
#' 
#' @export
#' 
#' 

Rt_linear <- function(inp,n_week = 24){
  
  ratio <- inp$period/inp$step
  
  if((ratio %% 1)!=0) { stop("ratio of period by step should be integer")}
  if(ratio %% 2) { stop("ratio of period by step should be even")}
  
  temp <- seq(inp$R_max, inp$R_min,length.out = ratio/2+1)
  temp <- c(temp,rev(temp[-1])[-1])
  
  
  n_rep_cycle <- n_week/(inp$period/7)
  if((n_rep_cycle %% 1)!=0) { stop("nb cycle repetitions should be integer (with 24 weeks simulation")}
  
  
  Rt <- data.frame(t = 1:(n_week*7), 
                   Rt = rep(rep(temp, each = inp$step), n_rep_cycle))
  
  
  return(Rt)
  
}
