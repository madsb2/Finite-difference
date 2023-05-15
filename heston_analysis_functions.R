########################################
############# Heston ###################
########################################

heston <- function(param, grid){
  startTime <- Sys.time()
  
  s_vector <- grid$s_vector
  v_vector <- grid$v_vector
  t_vector <- grid$t_vector
  
  J <- length(s_vector)
  M <- length(v_vector)
  N <- length(t_vector)
  
  K <- param$K
  r <- param$r
  eta <- param$eta
  rho <- param$rho
  kappa <- param$kappa
  sigma <- param$sigma
  
  U <- matrix(NA,M,J)
  for(j in 2:J){
    for(m in 1:M){
      U[m,j] <- Heston.Fourier(spot = s_vector[j],timetoexp = t_vector[N],strike = K,r = r,
                               divyield = 0,V = v_vector[m],theta = eta, kappa = kappa, epsilon = sigma,rho = rho)
    }
  }
  U[,1] <- 0
  
  endTime <- Sys.time()
  run_time <- difftime(endTime, startTime, units = 'mins')
  return(list("prices"=U,"run_time"=run_time))
}

################################################
############### Heston prices ##################
################################################

an_heston_prices <- function(cases, gridparam, m1, m2, unif, adi, case, wd){
  
  startTime <- Sys.time()
  setwd(wd)
  
  param <- cases[[case]]
  if(adi == 'MCS'){
    param$theta <- 1/3
  }else{
    param$theta <- 1/2
  }
  gridparam$Tmax <- param$Tmax
  
  for(i in 1:length(m1)){
    print(paste('STATUS:',(i-1)/length(m1)))
    
    #Update grid
    gridparam$J <- m1[i]+1
    gridparam$M <- m2[i]+1
    grid <- grid2(param, gridparam, unif)
    
    #Compute A matrix
    HES <- heston(param, grid)
    row.names(HES$prices) <- grid$v_vector; colnames(HES$prices) <- grid$s_vector
    saveRDS(HES, file = paste('HES_case',case,'_J',gridparam$J,'_M',gridparam$M,'_N',gridparam$N,'_unif',unif, sep = ""))
    
    if(adi == 'DO'){
      ADI <- Douglass(param, grid)
    }else if(adi == 'CS'){
      ADI <- Craig_sneyd(param, grid)
    }else if(adi == 'MCS'){
      ADI <- Modified_Craig_sneyd(param, grid)
    }
    else{print('Wrong ADI input')}
    
    row.names(ADI$prices) <- grid$v_vector; colnames(ADI$prices) <- grid$s_vector
    saveRDS(ADI, file = paste(adi,'_case',case,'_J',gridparam$J,'_M',gridparam$M,'_N',gridparam$N,'_unif',unif, sep = ""))
  }
  endTime <- Sys.time()
  run_time <- difftime(endTime, startTime, units = 'mins')
  return(print(run_time))
}
################################################
############### Heston errors ##################
################################################
an_heston_errors <- function(cases, gridparam, m1, m2, unif, adi, case, wd){
  setwd(wd)
  
  param <- cases[[case]]
  gridparam$Tmax <- param$Tmax
  
  #Create set of s and v
  mat_glob <- mat <- matrix(NA, length(m1), 3)
  
  for(i in 1:(length(m1))){
    #Update grid
    gridparam$J <- m1[i]+1
    gridparam$M <- m2[i]+1
    grid <- grid2(param, gridparam, unif = unif)
    
    set <- sv_set(param, grid)
    
    #Read files
    fd_an  <- readRDS(paste('HES_case',case,'_J',gridparam$J,'_M',gridparam$M,'_N',gridparam$N,'_unif',unif, sep = ""))
    fd_adi <- readRDS(paste(adi,'_case',case,'_J',gridparam$J,'_M',gridparam$M,'_N',gridparam$N,'_unif',unif, sep = ""))
    
    mat[i,1] <- err_l2(fd_an$prices[set$v,set$s], fd_adi$prices[set$v,set$s])
    mat[i,2] <- err_inf(fd_an$prices[set$v,set$s], fd_adi$prices[set$v,set$s])
    
    mat_glob[i,1] <- err_l2(fd_an$prices, fd_adi$prices)
    mat_glob[i,2] <- err_inf(fd_an$prices, fd_adi$prices)
    
    mat[i,3] <- fd_adi$run_time
  }
  
  return(list(region = mat, global = mat_glob))
}

################################################
############### Heston prices ##################
################################################

an_heston_dt <- function(cases, gridparam, n, unif, adi, case, wd){
  
  startTime <- Sys.time()
  setwd(wd)
  
  param <- cases[[case]]
  if(adi == 'MCS'){
    param$theta <- 1/3
  }else{
    param$theta <- 1/2
  }
  gridparam$Tmax <- param$Tmax
  
  for(i in 1:length(n)){
    print(paste('STATUS:',(i-1)/length(n)))
    
    #Update grid
    gridparam$N <- n[i]+1
    grid <- grid2(param, gridparam, unif)
    
    #Compute A matrix
    HES <- heston(param, grid)
    row.names(HES$prices) <- grid$v_vector; colnames(HES$prices) <- grid$s_vector
    saveRDS(HES, file = paste('HES_case',case,'_J',gridparam$J,'_M',gridparam$M,'_N',gridparam$N,'_unif',unif, sep = ""))
    
    if(adi == 'DO'){
      ADI <- Douglass(param, grid)
    }else if(adi == 'CS'){
      ADI <- Craig_sneyd(param, grid)
    }else if(adi == 'MCS'){
      ADI <- Modified_Craig_sneyd(param, grid)
    }else{print('Wrong ADI input')}
    
    row.names(ADI$prices) <- grid$v_vector; colnames(ADI$prices) <- grid$s_vector
    saveRDS(ADI, file = paste(adi,'_case',case,'_J',gridparam$J,'_M',gridparam$M,'_N',gridparam$N,'_unif',unif, sep = ""))
  }
  endTime <- Sys.time()
  run_time <- difftime(endTime, startTime, units = 'mins')
  return(print(run_time))
}

################################################
############## Heston errors dt ################
################################################
an_heston_errors_dt <- function(cases, gridparam, n , unif, adi, case, wd){
  setwd(wd)
  
  param <- cases[[case]]
  gridparam$Tmax <- param$Tmax
  
  #Create set of s and v
  mat_glob <- mat <- matrix(NA, length(n), 3)
  
  for(i in 1:(length(n))){
    #Update grid
    gridparam$N <- n[i]+1
    grid <- grid2(param, gridparam, unif = unif)
    
    set <- sv_set(param, grid)
    
    #Read files
    fd_an  <- readRDS(paste('HES_case',case,'_J',gridparam$J,'_M',gridparam$M,'_N',gridparam$N,'_unif',unif, sep = ""))
    fd_adi <- readRDS(paste(adi,'_case',case,'_J',gridparam$J,'_M',gridparam$M,'_N',gridparam$N,'_unif',unif, sep = ""))
    
    mat[i,1] <- err_l2(fd_an$prices[set$v,set$s], fd_adi$prices[set$v,set$s])
    mat[i,2] <- err_inf(fd_an$prices[set$v,set$s], fd_adi$prices[set$v,set$s])
    
    mat_glob[i,1] <- err_l2(fd_an$prices, fd_adi$prices)
    mat_glob[i,2] <- err_inf(fd_an$prices, fd_adi$prices)
    
    mat[i,3] <- fd_adi$run_time
  }
  
  return(list(region = mat, global = mat_glob))
}

#######################################
############## Price 0 ################
#######################################

an_heston_price0_sol <- function(param, gridparam, s0, v0){

  N <- gridparam$N
  Tmax <- gridparam$Tmax
  
  K <- param$K
  r <- param$r
  eta <- param$eta
  rho <- param$rho
  kappa <- param$kappa
  sigma <- param$sigma
  
  price0 <- Heston.Fourier(spot = s0,timetoexp = Tmax,strike = K,r = r,
                         divyield = 0,V = v0,theta = eta, kappa = kappa, epsilon = sigma,rho = rho)
  return(price0)
}

an_heston_price0_fd <- function(grid, fd, s0, v0){
  
  s_vector <- grid$s_vector
  v_vector <- grid$v_vector
  
  #estimate specific price
  i_sL <- which.min(abs(s_vector - s0))-1; i_sU <- i_sL+1
  i_vL <- which.min(abs(v_vector - v0))-1; i_vU <- i_vL+1
  
  #paste(s_vector[i_sL], s_vector[i_sL+1])
  
  sL <- s_vector[i_sL]; sU <- s_vector[i_sU]
  vL <- v_vector[i_vL]; vU <- v_vector[i_vU]
  
  f1 <- (vU-v0)/(vU-vL)*fd[i_vL,i_sL] + (v0-vL)/(vU-vL)*fd[i_vU,i_sL]
  f2 <- (vU-v0)/(vU-vL)*fd[i_vL,i_sU] + (v0-vL)/(vU-vL)*fd[i_vU,i_sU]
  
  price0 <- (sU-s0)/(sU-sL)*f1 + (s0-sL)/(sU-sL)*f2
  
  return(price0)
}


