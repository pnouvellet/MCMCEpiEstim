#' truncate incidence
#'
#' Remove incidence (i.e. replace with NA) after a given time or once a certain incidence threshold is reached.
#' 
#' 
#' @param t_truncate integer, truncate incidence after t_truncate* t_window days 
#'
#' @param incidence_truncate integer, set incidence at NA for all time windows after which
#' 'incidence_truncate' daily cases are reached
#'
#' @return A data.frame similar to the input incidence but with truncation applied.
#'
#'
#' @export 
#' 
#' 

truncate_incidence <- function(I0_t_import, I, t_window, 
                              t_truncate = NULL, incidence_truncate = NULL){
  
  # check incidence/time to truncate
  if (!is.null(incidence_truncate)){
    f <- apply(I[,-1],2,function(x) which(x>incidence_truncate)[1])
    f <- ceiling((f-I0_t_import)/7)*7 + I0_t_import+1
    for(i in 2:nrow(I)){
      if(!is.na(f[i-1])){
        I[f[i-1]:nrow(I),i] <- NA
      }
    }
  }
  
  if (!is.null(t_truncate)){
    I <- I[1:(t_truncate*t_window+I0_t_import),]
  }
  
 
  
  return(I)
}