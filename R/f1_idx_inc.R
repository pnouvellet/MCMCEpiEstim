#' useful for indices
#' 
#' @param x vector of 2 index integers
#' 
#' 
#' @return a sequence of indices betwen the first and second number 
#' 
#' @export
#' 

f1_idx_inc <- function(x){
  x[1]:x[2]
}

f2_idx_inc <- function(x){
  rep(x[2],length(x[1]:x[2]))
}
