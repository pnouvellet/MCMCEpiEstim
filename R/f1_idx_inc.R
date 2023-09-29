#' useful for indices
#' 
#' @param x vector of 2 index integers
#' 
#' 
#' @details return corrected variances to decrease/increase acceptance
#' 
#' @export
#' 

f1_idx_inc <- function(x){
  x[1]:x[2]
}
