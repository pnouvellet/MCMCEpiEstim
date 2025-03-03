#' Truncate incidence
#'
#' Remove incidence (i.e. replace with NA) after a given time or once a certain incidence threshold is reached.
#' 
#' @param I0_t_import integer, which of the initial incidence is imported. 
#' 
#' @param I data.frame of incidence. ncol: number of locations +1 (\code{$time}), nrow: time
#'
#' @param t_window integer, time window during which Rt is assumed constant.
#' 
#' @param t_truncate integer, truncate incidence after \code{t_truncate* t_window days + I0_t_import}.
#'
#' @param incidence_truncate integer, set incidence at NA for all time windows after which
#' 'incidence_truncate' daily cases are reached.
#'
#' @return A data.frame similar to the input incidence but with truncation applied.
#'
#'
#' @export 
#' 
#' 
#' @examples
#' 
#' I <- data.frame(t = seq(1:100))
#' I$loc1 <- 2*exp(.04*I$t)
#' I$loc2 <- 2*exp(.08*I$t)
#' matplot(I[,-1], xlim = c(0,100), ylim = c(0,220))
#' 
#' I_trunc1 <- truncate_incidence(I0_t_import = 1, I = I,t_window = 5,t_truncate = 10)
#' matplot(I_trunc1[,-1], xlim = c(0,100), ylim = c(0,220))
#' 
#' I_trunc2 <- truncate_incidence(I0_t_import = 1, I = I,t_window = 5,incidence_truncate = 100)
#' matplot(I_trunc2[,-1], xlim = c(0,100), ylim = c(0,220))
#'

truncate_incidence <- function(I0_t_import, I, t_window, 
                              t_truncate = NULL, incidence_truncate = NULL){
  
  # check if a parameter is inputed for the incidence threshold (threshold on number of cases)
  if (!is.null(incidence_truncate)){
    if(ncol(I)==2){
      # check if incidence reach the threshold
      f <- which(I[,2]>incidence_truncate)[1]
    } else {
      # check if incidence reach the threshold, for multiple locations
      f <- apply(I[,-1],2,function(x) which(x>incidence_truncate)[1])
    }
    # identify the begining of the next time-window, after which the threshold is reached
    f <- ceiling((f-I0_t_import)/7)*7 + I0_t_import+1
    # replace incidence after the threshold is reached by NA
    for(i in 2:nrow(I)){
      if(!is.na(f[i-1])){
        I[f[i-1]:nrow(I),i] <- NA
      }
    }
  }
  # check if a parameter is inputed for t_truncate (threshold on time, after t_tuncate, set incidence to NA)
  if (!is.null(t_truncate)){
    I <- I[1:(t_truncate*t_window+I0_t_import),]
  }
  
  return(I)
}