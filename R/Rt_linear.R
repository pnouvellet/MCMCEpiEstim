#' project
#'
#' wrapper for projection with/without overdispersion 
#' 
#' @param inp see others
#'
#' @details n_week nb of weeks for the simulation, default of 24
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
