#' Project incidence
#'
#' Wrapper for projection with/without overdispersion and with under-reporting.
#' 
#' @param I0 an incidence object, specifying the initial incidence. 
#' The initial incidence is here assumed the same in every location.
#' 
#' @param Rt A data.frame with time (\code{$t}) and associated daily Rts (instantaneous Rt, \code{$Rt}).
#'                   
#' @param n_loc number of locations.
#' 
#' @param t_max number of simulated time steps.
#'
#' @param p level of reporting, i.e. probability that a case simulated is actually observed.
#' p can take a single value and that would assume constant reporting across time and locations simulated.
#' Otherwise, p can be data.frame similar to Rt (\code{$t} for time and \code{$p} for the reporting). Such p data.frame
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
#' @param k.seed an optional argument to set the seed of the simulations
#'
#' @param threshold_cumCase an optional argument to remove location with a cumulative number of 
#' cases (observed) below threshold_cumCase. if some locations are removed, the function keep simulating new ones
#' until the the right number of simulations (n_loc) have reached the required cumulative number of cases.
#' 
#' @details Simulate incidence using the package 'projection', with a time-varying Rt, The possibility
#' to specify the presence of overdispersion and the presence of under-reporting of cases.
#' 
#' 
#' @return A list containing 2 data.frames of time (\code{$t}) and daily incidence values at each location 
#' (\code{$simX}, with X an integer for each location). 
#' The first data.frame, \code{$I_true}, contains the 'true' incidence, while the second, \code{$I_obs}, contains the 'observed'
#' incidence once under-reporting is accounted for.
#' 
#' @import projections
#' @export  
#' 
#' @examples
#' # a serial interal
#' si <- c(0,1)
#' # simualte an Rt
#' input <- data.frame(R_min = 0.95,
#'                     R_max = 1.1,
#'                     k = 1,
#'                     period = 24, 
#'                     step = 4,
#'                     I0 = 100,
#'                     pi = 0.5)
#' Rt <- Rt_linear(inp = input, n_week = input$period)
#' plot(Rt$t, Rt$Rt,type='l')
#' # set initial incidence
#' I0 <- incidence::as.incidence(x = input$I0, dates = 1, interval = 1)
#' 
#' # simulate for 10 locations
#' res <- MCMCEpiEstim::project_fct(I0 = I0,
#'                                  Rt = Rt,
#'                                  n_loc = 10,
#'                                  t_max = nrow(Rt),
#'                                  si = si,
#'                                  p = input$pi,
#'                                  model = 'negbin',
#'                                  over = input$k)
#' 
#' layout(matrix(1:2,nrow = 2))
#' matplot((res$I_true[,-1]), xlab = '', ylab = 'cases', main = 'true cases')
#' matplot((res$I_obs[,-1]), xlab = '', ylab = 'cases', main = 'observed cases')
#' 
#' 

project_fct <- function(I0, Rt, n_loc, t_max, si, p,
                        model = 'poisson',over = NULL, k.seed = NULL, 
                        threshold_cumCase = NULL){
  # check seed
  if(!is.null(k.seed)){
    set.seed(k.seed)
  }
  
  # allocate data.frame for results
  I_obs_final <- I_true_final <- as.data.frame(matrix(1:(t_max+1),
                                                      ncol = 1, nrow = t_max+1))
  # initialise number of location
  n_loc_above <- 0
  
  # loop process until enough locations reach the require number of total cases 
  # (see threshold_cumCase option above)
  while(n_loc_above < n_loc + 1 ){
    # initialise data.frame
    I <- as.data.frame(matrix(NA,ncol = n_loc+1, nrow = t_max+1))
    names(I) <- c('t',paste0('sim',1:n_loc))
    I$t <- 1:(t_max+1)
    
    # project, see projections R-package
    if (model == 'poisson'){
      temp <- as.data.frame(projections::project(x = I0,
                                    R = Rt$Rt, 
                                    si = si[-1], n_sim = n_loc, time_change = 1:(t_max-1-I0$timespan+1),
                                    n_days = t_max-I0$timespan+1, 
                                    R_fix_within = TRUE, 
                                    model = 'poisson',instantaneous_R = TRUE))
    }else if (model == 'negbin'){
      temp <- as.data.frame(projections::project(x = I0,
                                    R = Rt$Rt, 
                                    si = si[-1], n_sim = n_loc, time_change = 1:(t_max-1-I0$timespan+1),
                                    n_days = t_max-I0$timespan+1, 
                                    R_fix_within = TRUE, 
                                    model = 'negbin',instantaneous_R = TRUE, size = over))
    }
    
    # bind simulated incidence with that used to initialise the simulation (e.g. imported cases)
    I[1:I0$timespan,] <- cbind(I0$dates,matrix(data = I0$counts,
                                               nrow = nrow(I0),
                                               ncol = n_loc,byrow = FALSE))
    
    I[(I0$timespan+1):(nrow(I)),] <- temp
    
    I_obs <- I
    
    # extract the reporting (either single value throughout, or probability per day)
    if(length(p)>1){
      p_value <- p$pi
    } else{
      p_value <- p
    }
    
    # assume binomial sampling of simulated incidence with reporting
    temp <- matrix(rbinom(n = (t_max+nrow(I0))*n_loc,
                          size = as.matrix(I[,-1]),
                          prob = matrix(p_value,nrow = t_max+nrow(I0),ncol = n_loc,byrow = FALSE) ),
                   nrow = t_max+nrow(I0),ncol = n_loc,byrow = FALSE)
    
    I_obs[,-1] <- temp
    
    # search for locations where total simulated incidence is below the required threshold # (see threshold_cumCase option above)
    if(is.null(threshold_cumCase)){
      f_CumCase_above <- 1:n_loc
    }else{
      f_CumCase_above <- which(colSums(I_obs[-1,-1]) > threshold_cumCase)
    }
    
    # remove locations not meeting the criteria set
    I_obs[,-1] <- I_obs[,f_CumCase_above+1]
    I[,-1] <- I[,f_CumCase_above+1]
    
    I_obs_final <- cbind(I_obs_final, I_obs[,-1])
    I_true_final <- cbind(I_true_final, I[,-1])
    
    # update the number of 'correctly' simulated locations
    n_loc_above <- ncol(I_obs_final)
  }
    
  # store results
  I_obs_final <-  I_obs_final[,1:(n_loc+1)]
  I_true_final <- I_true_final[,1:(n_loc+1)]
  names(I_obs_final) <- c('t',paste('sim',1:n_loc))
  names(I_true_final) <- c('t',paste('sim',1:n_loc))
  
  return(list(I_true = I_true_final, I_obs = I_obs_final))
  
}
