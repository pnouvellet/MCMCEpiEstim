#' diagnostic plots 
#'
#' diagnostic plots 
#' 
#' @param logged plot y_axis on natural scale, i.e. FALSE, or log10 scale, logged = TRUE,
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

diag_plot <- function(res, logged, max_x=1, dist , k, res_MCMC=NULL){
  
  check <- c()
  for (i in 1:length(res)){
    tmp <- res[[i]]
    f <- which(!is.na(tmp$t_start))
    
    idx_inc <- c(apply(cbind(tmp$t_start[f],tmp$t_end[f]),1,f1_idx_inc))
    idx_Rt <- c(apply(cbind(tmp$t_start[f],tmp$t_end[f]),1,f2_idx_inc))
    
    if(!is.null(res_MCMC)){
      if(res_MCMC$input_MCMC$param_agg){
        tmp$`Mean(R)`[f] <- apply(res_MCMC$theta_R_thinned,2,mean)
      }else{
        n <- length(f)
        tmp$`Mean(R)`[f] <- apply(res_MCMC$theta_R_thinned[,((i-1)*n+1):(n*i)],2,mean)
      }
    }
    check0 <- data.frame(time = tmp$t[idx_inc],
                        location = paste('loc_',i),
                        Rt = tmp$`Mean(R)`[idx_Rt],
                        Obs = tmp$incidence[idx_inc],
                        Exp = tmp$Oi[idx_inc] * tmp$`Mean(R)`[idx_Rt])

    check <- rbind(check,check0)
  }
  
  check$residual <- check$Obs - check$Exp
    
  # 95%CI associated with particular observation
  if (dist == 'poisson'){
    if(!is.null(res_MCMC)){
      p_reps <- res_MCMC$input_MCMC$p_reps
    }else{
      p_reps <- 1
    }
    if(p_reps ==1){
      check$lim_up <- qpois(p = 0.975,lambda = check$Exp,lower.tail = TRUE) - check$Exp
      check$lim_down <- qpois(p = 0.025,lambda = check$Exp,lower.tail = TRUE) - check$Exp
    }else{
      check$var <- check$Exp * (1 + (1- p_reps) * check$Rt) +
        (1- p_reps) * check$Rt * (1 + check$Rt)
      
      check$lim_up <- qnbinom(p = 0.975,mu = check$Exp, size = check$Exp^2 / (check$var - check$Exp), lower.tail = TRUE) - check$Exp
      check$lim_down <- qnbinom(p = 0.025,mu = check$Exp, size = check$Exp^2 / (check$var - check$Exp),lower.tail = TRUE) - check$Exp
    }
  }else if(dist == 'nb'){
    p_reps <- res_MCMC$input_MCMC$p_reps
    # check$var <- check$Exp * (1 + check$Rt/k)
    check$var <- check$Exp * (1 + (1- p_reps) * check$Rt  + p_reps * check$Rt/k) +
      (1- p_reps) * check$Rt * (1 + check$Rt + p_reps * check$Rt/k)
    
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
    var_max <- x*(1+max(check$Rt,na.rm = TRUE)/k)
    var_min <- x*(1+min(check$Rt,na.rm = TRUE)/k)
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
  check$int_R <- cut(check$Rt,breaks = seq(min(check$Rt,na.rm = TRUE),max(check$Rt,na.rm = TRUE),length.out=5))
  check$int_t <- cut(check$time,breaks = seq(min(check$time,na.rm = TRUE),max(check$time,na.rm = TRUE),length.out=5))
  
  
  
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
  
  return(list(res_table = res_table, check = check))
}

