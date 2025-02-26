#' Tuning
#'
#' Tune the variance proposal.
#' Try to tune variance toward 20% acceptance.
#' 
#' @param theta matrix of posterior samples of posterior distribution. ncol: nb parameters, nrow: nb samples
#' 
#' @param s vector, proposal variances used to obtain posterior samples theta (size the same as ncol of theta)
#'                   
#' @param it integer, number of theta sample
#'
#' 
#' @return a vector of corrected variances to decrease/increase acceptance
#' 
#' @export
#' 

adapt <- function(theta,s,it){
  Acc <- colSums(diff(theta)!=0)/(it-1)   # current acceptance rate
  s_out <- Acc/.2*s                       # reduce variance if acceptance was too low, or increase variance if accpetance was too high (work well!)
  s_out[which(Acc==0)] <- s[which(Acc==0)]/2   # if never accepted half the variance of the proposal 
  return(s_out)
}

