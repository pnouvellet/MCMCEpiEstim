#' epiestim wrap
#'
#' wrapper to call EpiEstim independently on multiple locations
#' 
#' @param I0 matrix, samples of posterior distribution. ncol: nb parameters, nrow: nb samples
#' 
#' @param I matrix, samples of posterior distribution. ncol: nb parameters, nrow: nb samples
#' 
#' @param t_window matrix, samples of posterior distribution. ncol: nb parameters, nrow: nb samples
#' 
#' @param mean_prior single real number, mean prior for Rts
#' 
#' @param std_prior single real number, stad deviation for prior for Rts
#' 
#' @param si serial distribution (as in EpiEstim include a 0 weighted SI on same day)
#'                   
#' @param overlap TRUE/FALSE, whether using overlapping time window or not
#'
#' 
#' @details res, list of 1 dataframe per location summarising incidence, time, overall infectivity,
#'              time window start and end date, estimated mean/std deviation/median and
#'               95% CI lower and upper quantile for Rts
#' @export
#' 
#' 
# 

fct_EpiEstim <- function(I0 , I , t_window, 
                         mean_prior, std_prior,
                         si, overlap){
  
  t_max <- nrow(I)
  n_sim <- ncol(I)-1
  # time windows
  t_start <- seq(I0$timespan+1, t_max-t_window+1,by = 1)
  if (overlap == FALSE) {
    t_start <- t_start[seq(1,length(t_start),by = t_window)]
  }else{
    
  }
  t_end <- t_start + t_window - 1      
  
  # defines settings for estimate_R 
  config <- make_config(list(si_distr = si,
                             t_start = t_start,
                             t_end = t_end,
                             mean_prior = mean_prior,
                             std_prior = std_prior))
  
  res <- list()
  for (i in 1:n_sim){
    
    d_incidence <- data.frame(t = I[,1],
                              incidence = I[,i+1])
    
    # overall Infectivity
    d_incidence$Oi = overall_infectivity(incid = d_incidence$incidence,
                                         si_distr = si)
    d_incidence$Oi[1:I0$timespan] <- NA
    f_0 <- which(d_incidence$Oi == 0)
    
    I_corr <- data.frame(local = d_incidence$incidence,
                    imported = 0)
    # correct initial case as imported
    I_corr$imported[1:I0$timespan] <- d_incidence$incidence[1:I0$timespan]
    I_corr$local[1:I0$timespan] <- 0
    # correction for days with 0 overall infectivity
    I_corr$imported[f_0] <- d_incidence$incidence[f_0]
    I_corr$local[f_0] <- 0
    
    # estimation
    res_non_parametric_si <- estimate_R(I_corr, 
                                        method = "non_parametric_si",
                                        config = config)
    
    Repi <- res_non_parametric_si$R[,c("t_start","t_end",
                                       "Mean(R)","Std(R)",
                                       "Quantile.0.025(R)","Median(R)","Quantile.0.975(R)")]
    
    d_incidence[,c("t_start","t_end",
                   "Mean(R)","Std(R)",
                   "Quantile.0.025(R)","Median(R)","Quantile.0.975(R)")] <- NA
    
    d_incidence[Repi$t_end,c("t_start","t_end",
                             "Mean(R)","Std(R)",
                             "Quantile.0.025(R)","Median(R)","Quantile.0.975(R)")] <- Repi
    # f <- which(d_incidence$incidence>0 & d_incidence$Oi==0)
    # d_incidence[f,c("Mean(R)","Std(R)",
    #                 "Quantile.0.025(R)","Median(R)","Quantile.0.975(R)")] <- NA
    
    res[[i]] <- d_incidence
  }
  
  return(res)
  # plot(d_incidence$t -(t_window-1)/2,d_incidence$`Mean(R)`)
  # lines(Rt$t,Rt$Rt,col='red')
}
