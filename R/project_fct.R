#' project
#'
#' wrapper for projection with/without overdispersion and with under-reporting
#' 
#' @param I0 an incidence object, specifying the initial incidence. 
#' The initial incidence is here assumed the same in every location.
#' 
#' @param Rt A data.frame with time ($t) and associated daily Rts (instantaneous Rt, $Rt).
#'                   
#' @param n_loc number of locations.
#' 
#' @param t_max number of simulated time steps.
#'
#' @param p level of reporting, i.e. probability that a case simulated is actually observed.
#' p can take a single value and that would assume constant reporting across time and locations simulated.
#' Otherwise, p can be data.frame similar to Rt ($t for time and $p for the reporting). Such p data.frame
#' would be longer (1 more value) than Rt as it would be specifying in the reporting for imported cases seeding 
#' the simulation.
#' In one of our application, we first specify Rt and then then simulate reporting with a single time-step change 
#' in reporting:
#' 
#' t_change_p <- input$t_pi_change*input$step+1
#' 
#' Pi_t <- data.frame(t = c(Rt$t,tail(Rt$t,1)+1),
#'                    pi = c(rep(input$pi,t_change_p),
#'                          rep(input$pi2, nrow(Rt)+1-t_change_p )) )
#'  input$pi, input$pi2, input$t_pi_change specifying the 2 reporting probability and last time (in weeks) 
#'  when pi1 applies respectively.
#'
#' @param si serial distribution (as in EpiEstim include a 0 weighted SI on same day).
#'
#' @param model takes 'poisson' or 'negbin', i.e. offspring distribution assumed. 
#' By default, 'poisson' is assumed.
#'
#' @param over a real number quantifying the level of overdispersion in the offspring
#'  distribution. The overdispersion is defined in the context of a negative binomial distribution 
#'  as: if X follows a negative binomial distribution with mean \mu and overdispersion \delta, then the 
#'  $$ E[X] = \mu $$ and $$ Var[X] = \mu + frac{\mu^2}{\delta} $$.
#'
#' 
#' @details Simulate incidence using the package 'projection', with a time-varying Rt, The possibility
#' to specify the presence of overdispersion and the presence of under-reporting of cases.
#' 
#' 
#' @return A list containing 2 data.frames of time ($t) and daily incidence values at each location 
#' ($simX, with X an integer for each location). 
#' The first data.frame, I_true, contains the 'true' incidence, while the second, I_obs, contains the 'observed'
#' incidence once under-reporting is accounted for.
#' 
#' @export
#' 
#' 

project_fct <- function(I0, Rt, n_loc, t_max, si, p,
                        model = 'poisson',over = NULL, k.seed = NULL, 
                        threshold_cumCase = NULL){
  
  if(!is.null(k.seed)){
    set.seed(k.seed)
  }
  
  I_obs_final <- I_true_final <- as.data.frame(matrix(1:(t_max+1),
                                                      ncol = 1, nrow = t_max+1))
  n_loc_above <- 0
  
  while(n_loc_above < n_loc){
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
    
    I_obs[,-1] <- temp
    
    if(is.null(threshold_cumCase)){
      f_CumCase_above <- 1:n_loc
    }else{
      f_CumCase_above <- which(colSums(I_obs[-1,-1]) > threshold_cumCase)
    }
    
    I_obs <- I_obs[,f_CumCase_above]
    I <- I[,f_CumCase_above]
    
    I_obs_final <- cbind(I_obs_final, I_obs[,-1])
    I_true_final <- cbind(I_true_final, I[,-1])
    
    n_loc_above <- ncol(I_obs_final)
  }
    
  I_obs_final <-  I_obs_final[,1:(n_loc+1)]
  I_true_final <- I_true_final[,1:(n_loc+1)]
  names(I_obs_final) <- c('t',paste('sim',1:n_loc))
  names(I_true_final) <- c('t',paste('sim',1:n_loc))
  
  return(list(I_true = I_true_final, I_obs = I_obs_final))
  
}
