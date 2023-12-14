#' diagnostic plots 
#'
#' diagnostic plots 
#' 
#' @param logged plot y_axis on natural scale, i.e. FALSE, of log10 scale, logged = TRUE,
#' 
#' @param max_x plot up to x*max(x), default show the whole axis, i.e. max_x = 1
#'                   
#' @param dist assumed offspring distribution: 'poisson' or 'nb'
#' 
#' @param Rt data.frame with t = time and vector of Rt's assumed or estimated
#'
#' @param k relevant only when dist = 'nb', assumed or estimated overdispersion (single estimate)
#'
#' 
#' @details return corrected variances to decrease/increase acceptance
#' 
#' @export
#' 

diag_plot <- function(I_NB, logged, max_x=1, dist, Rt, k){
  
  # get overall infertivity from EpiEstim
  # E_NB <- matrix(unlist(lapply(res_EpiEstim, "[", ,"Oi")), 
  #                nrow = t_max, ncol = n_location, byrow = FALSE)
  if(ncol(Rt) == 2){
    Rts <- matrix(Rt$Rt,nrow = nrow(Rt),ncol = ncol(I_NB)-1,byrow = FALSE)
  }else{
    Rts <- as.matrix(Rt[,-1])
  }
  
  E_NB <- I_NB[1:(nrow(I_NB)-1),-1] * Rts
  E_NB <- as.matrix(rbind(NA,E_NB))
  
  check <- data.frame(time = rep(as.matrix(I_NB[-1,1]),ncol(I_NB)-1),
                      location = rep(paste('loc_',1:(ncol(I_NB)-1)), each = nrow(I_NB)-1),
                      Rt = c(Rts),
                      Obs = c(as.matrix(I_NB[-1,-1])),
                      Exp = c(E_NB[-1,]),
                      residual = NA, 
                      lim_up = NA,
                      lim_down = NA)
  
  check$residual <- check$Obs - check$Exp
    
  # 95%CI associated with particular observation
  if (dist == 'poisson'){
    check$lim_up <- qpois(p = 0.975,lambda = check$Exp,lower.tail = TRUE) - check$Exp
    check$lim_down <- qpois(p = 0.025,lambda = check$Exp,lower.tail = TRUE) - check$Exp
  }else if(dist == 'nb'){
    check$var <- check$Exp * (1 + check$Rt/k)
    check$lim_up <- qnbinom(p = 0.975,mu = check$Exp, size = check$Exp^2 / (check$var - check$Exp), lower.tail = TRUE) - check$Exp
    check$lim_down <- qnbinom(p = 0.025,mu = check$Exp, size = check$Exp^2 / (check$var - check$Exp),lower.tail = TRUE) - check$Exp
  }
  
  #remove all entries where expected case number == 0
  check[which(check$Exp==0),] <- NA
  check$log10Exp <- log10(check$Exp)
  
  # return true when observation lies within 95%CI
  check$out <- ((check$residual<=check$lim_up )+ (check$residual >=check$lim_down))==2
  

  
  # expected mean - 95%CI
  x <- seq(1,max(check$Exp,na.rm = TRUE), length.out = 1e3) # range of expected incidence
  if(dist == 'poisson'){
    y1 <- qpois(p = 0.975,lambda = x,lower.tail = TRUE)
    y2 <- qpois(p = 0.025,lambda = x,lower.tail = TRUE)
  }else if(dist == 'nb'){
    # # with nb
    var_max <- x*(1+max(check$Rt)/k)
    var_min <- x*(1+min(check$Rt)/k)
    ynb1_max <- qnbinom(p = 0.975,mu = x, size = x^2/(var_max-x),lower.tail = TRUE)
    ynb2_max <- qnbinom(p = 0.025,mu = x, size = x^2/(var_max-x),lower.tail = TRUE)
    ynb1_min <- qnbinom(p = 0.975,mu = x, size = x^2/(var_min-x),lower.tail = TRUE)
    ynb2_min <- qnbinom(p = 0.025,mu = x, size = x^2/(var_min-x),lower.tail = TRUE)
    
  }
  
  max_x <- max(c(check$Obs,check$Exp),na.rm = TRUE) * max_x
  min_x <- min(check$Exp,na.rm = TRUE)
  if(logged){
    Xobs <- check$log10Exp
    Xci <- log10(x)
    xlim <- log10(c(min_x,max_x) )
  }else{
    Xobs <- check$Obs
    Xci <- (x)
    xlim <- c(0,max_x) 
  }
  
  plot(Xobs,c(check$resid), xlim = xlim)
  abline(h = 0,col = 'red',lty = 1)
  
  if(dist=='poisson'){
    lines(Xci,y1-x,col = 'red',lty = 1)   
    lines(Xci,y2-x,col = 'red',lty = 1)   
  }else if(dist == 'nb'){
    lines(Xci,ynb1_max-x,col = 'red',lty = 1)   
    lines(Xci,ynb2_max-x,col = 'red',lty = 1)   
    lines(Xci,ynb1_min-x,col = 'blue3',lty = 1)   
    lines(Xci,ynb2_min-x,col = 'blue3',lty = 1)   
  }
  
  # confidence interval inclusion
  if (logged){
    check$int_I <- cut(check$log10Exp,breaks = seq(log10(min_x)*1.1,log10(max_x),length.out=5))
  }else{
    check$int_I <- cut(check$Exp,breaks = seq(0,max(check$Exp,na.rm = TRUE),length.out=5))
  }
  check$int_R <- cut(check$Rt,breaks = seq(min(check$Rt,na.rm = TRUE)*.9,max(check$Rt,na.rm = TRUE),length.out=5))
  check$int_t <- cut(check$time,breaks = seq(min(check$time,na.rm = TRUE)*.9,max(check$time,na.rm = TRUE),length.out=5))
  
  
  
  res_table <- list()
  # I
  tt <- table(check$int_I,check$out)
  temp <- tt[,2]/rowSums(tt)
  res_table$I <- data.frame( interval_I = c('all',names(temp)),
                    n = c(sum(rowSums(tt)),rowSums(tt)),
                    out = c(sum(tt[,1]),tt[,1]),
                    incl_p = c(sum(check$out,na.rm = TRUE)/sum(!is.na(check$out)),
                               temp))
  
  # Rt
  tt <- table(check$int_R,check$out)
  temp <- tt[,2]/rowSums(tt)
  res_table$R <- data.frame( interval_R = c('all',names(temp)),
                             n = c(sum(rowSums(tt)),rowSums(tt)),
                             out = c(sum(tt[,1]),tt[,1]),
                             incl_p = c(sum(check$out,na.rm = TRUE)/sum(!is.na(check$out)),
                                        temp))
  
  # t
  tt <- table(check$int_t,check$out)
  temp <- tt[,2]/rowSums(tt)
  res_table$t <- data.frame( interval_t = c('all',names(temp)),
                             n = c(sum(rowSums(tt)),rowSums(tt)),
                             out = c(sum(tt[,1]),tt[,1]),
                             incl_p = c(sum(check$out,na.rm = TRUE)/sum(!is.na(check$out)),
                                        temp))
  
  # location
  tt <- table(check$location,check$out)
  temp <- tt[,2]/rowSums(tt)
  res_table$loc <- data.frame( interval_loc = c('all',names(temp)),
                             n = c(sum(rowSums(tt)),rowSums(tt)),
                             out = c(sum(tt[,1]),tt[,1]),
                             incl_p = c(sum(check$out,na.rm = TRUE)/sum(!is.na(check$out)),
                                        temp))
  
  return(res_table)
}

