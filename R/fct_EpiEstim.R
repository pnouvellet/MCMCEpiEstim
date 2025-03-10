#' EpiEstim wrapper
#'
#' Wrapper to call EpiEstim independently on multiple locations.
#' 
#' @param I0_t_import which of the initial incidence is imported.
#' 
#' @param I dataframe of incidence. ncol: nb of location +1 (time), nrow: time.
#' 
#' @param t_window integer, time window during which Rt is assumed constant.
#' 
#' @param mean_prior single real number, mean prior for Rts.
#' 
#' @param std_prior single real number, std deviation for prior for Rts.
#' 
#' @param si serial distribution (as in EpiEstim include a 0 weighted SI on same day).
#'                   
#' @param overlap TRUE/FALSE, whether using overlapping time window or not.
#'
#' @details Use EpiEstim independently on each location to estimate Rts without
#'  accounting for overdispersion, nor under-reporting.
#' 
#' @return A list composed of 1 data.frame per location summarizing daily time, 
#' daily incidence (full incidence time series, including the imported cases for seeding), 
#' daily overall infectivity,the start and end of the time windows associated with 
#' the estimated mean/std deviation/median and
#' 95% CI lower and upper quantile for Rts 
#' (i.e. Mean(R)", "Std(R)", "Quantile.0.025(R)", "Median(R)", "Quantile.0.975(R)").
#' If no overlapping windows is used, the entries for t_start/t_end and estimated Rt characteristics are only given 
#' when 't' is the end of the time window.
#' 
#' @import EpiEstim
#' @export
#' 
#' @examples
#' 
#' I <- data.frame(t = seq(1:50))
#' I$loc1 <- 2*exp(.05*I$t)
#' I$loc2 <- 2*exp(.1*I$t)
#' matplot(I[,-1], xlim = c(0,50), ylim = c(0,200))
#'
#' 
#' res <- fct_EpiEstim(I0_t_import = 1 , I = I , t_window = 5, 
#'                                 mean_prior = 1, std_prior = 1,
#'                                 si = c(0,1), overlap = FALSE)
#' 
#' a <- res[[1]]
#' Hmisc::errbar(x =a$t_start, 
#'               y = a$`Median(R)`,
#'               yplus = a$`Quantile.0.975(R)`,
#'               yminus = a$`Quantile.0.025(R)`,
#'               xlab = 'time', ylab = 'Rt')
#' 
#' a <- res[[2]]
#' Hmisc::errbar(x =a$t_start, 
#'               y = a$`Median(R)`,
#'               yplus = a$`Quantile.0.975(R)`,
#'               yminus = a$`Quantile.0.025(R)`,
#'               xlab = 'time', ylab = 'Rt')
#'


fct_EpiEstim <- function(I0_t_import , I , t_window, 
                         mean_prior, std_prior,
                         si, overlap){
  
  input_import <- FALSE
  if (!is.data.frame(I)){
    input_import <- TRUE
    I_import <- I$import
    I_local <- I$local
    I <- I_import
    I[,-1] <- I[,-1] + I_local[,-1]
  }
  
  t_max <- nrow(I)
  n_sim <- ncol(I)-1
  # time windows
  t_start <- seq(I0_t_import+1, t_max-t_window+1,by = 1)
  if (overlap == FALSE) {
    t_start <- t_start[seq(1,length(t_start),by = t_window)]
  }else{
    
  }
  t_end <- t_start + t_window - 1      
  
  # defines settings for estimate_R 
  config <- EpiEstim::make_config(list(si_distr = si,
                             t_start = t_start,
                             t_end = t_end,
                             mean_prior = mean_prior,
                             std_prior = std_prior))
  
  res <- list()
  for (i in 1:n_sim){
    
    d_incidence <- data.frame(t = I[,1],
                              incidence = I[,i+1])
    # overall Infectivity
    f_incidence_nonNA <- which(!is.na(d_incidence$incidence))
    # replace NAs by 0
    d_incidence$incidence[-f_incidence_nonNA] <- 0
    d_incidence$Oi <- 0
    
    d_incidence$Oi[f_incidence_nonNA] <- EpiEstim::overall_infectivity(incid = d_incidence$incidence[f_incidence_nonNA], 
                                                                       si_distr = si)
    # remove infectivity estimate when we know case are imported
    d_incidence$Oi[1:I0_t_import] <- NA
    
    I_corr <- data.frame(local = d_incidence$incidence,
                         imported = 0)
    if(input_import == FALSE){
      # correct initial case as imported
      I_corr$imported[1:I0_t_import] <- d_incidence$incidence[1:I0_t_import]
      I_corr$local[1:I0_t_import] <- 0
    }else{
      I_corr$imported <- I_import[,i+1]
      I_corr$local <- I_local[,i+1]
      
      # then account for importation of simulation
      I_corr$imported[1:I0_t_import] <- d_incidence$incidence[1:I0_t_import]
      I_corr$local[1:I0_t_import] <- 0
    }
    
    # correction for days with 0 overall infectivity

      f_0 <- which(d_incidence$Oi == 0)
      I_corr$imported[f_0] <- d_incidence$incidence[f_0]
      I_corr$local[f_0] <- 0

    
    # estimation
    res_non_parametric_si <- estimate_R(I_corr, 
                                        method = "non_parametric_si",
                                        config = config)
    # output with key moments/quantiles
    Repi <- res_non_parametric_si$R[,c("t_start","t_end",
                                       "Mean(R)","Std(R)",
                                       "Quantile.0.025(R)","Median(R)","Quantile.0.975(R)")]
    
    d_incidence[,c("t_start","t_end",
                   "Mean(R)","Std(R)",
                   "Quantile.0.025(R)","Median(R)","Quantile.0.975(R)")] <- NA
    
    d_incidence[Repi$t_end,c("t_start","t_end",
                             "Mean(R)","Std(R)",
                             "Quantile.0.025(R)","Median(R)","Quantile.0.975(R)")] <- Repi
    
    res[[i]] <- d_incidence
  }
  
  return(res)
}
