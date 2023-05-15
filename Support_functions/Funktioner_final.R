##############################################################
#############  Thomas algorithm ##############################
##############################################################

#a1: Lower diagonal
#b1: Central diagonal
#c1: Upper diagonal
#d1: b-vector

thomas <- function(a1,b1,c1,d1){
  N <- length(c1)
  
  c2 <- rep(0,N)
  d2 <- rep(0,N)
  temp <- rep(0,N)
  a1[1] <- 0 #Redundant
  c1[N] <- 0 #Redundant
  
  c2[1] <- c1[1]/b1[1]
  
  for(n in 2:(N-1)){
    c2[n] <- c1[n]/(b1[n]-a1[n]*c2[n-1])
  }
  d2[1] <- d1[1]/b1[1]
  for(n in 2:N){
    d2[n] <- (d1[n]-a1[n]*d2[n-1])/(b1[n]-a1[n]*c2[n-1])
  }
  
  temp[n] <- d2[n]
  for(n in (N-1):1){
    temp[n] <- d2[n]-c2[n]*temp[n+1]
  }
  
  return(temp)
}
#################################################################
########################## Call #################################
#################################################################

bs <- function(S,K,r,sigma,tau,option='call'){
  if(is.na(log(S/K))==TRUE){print(S)}##
  d1 = 1/(sigma*sqrt(tau))*(log(S/K)+(r+0.5*sigma^2)*tau)
  if(is.na(d1) == TRUE){d1 <- 0} #If log(1) we get NaN
  d2 = d1-sigma*sqrt(tau)
  if(option=='call'){
    res <- S*pnorm(d1) - exp(-r*tau)*K*pnorm(d2)
  }else{
    res <- exp(-r*tau)*K*pnorm(-d2) - S*pnorm(-d1)
  }
  return(res)
}

#Price matrix
bs_call <- function(param,grid,option='call'){
  x_vector <- grid$s_vector
  t_vector <- grid$t_vector
  
  K <- param$K
  r <- param$r
  sigma <- param$sigma

  J = length(x_vector)
  N = length(t_vector)
  tau = t_vector[N]-t_vector
  
  f <- matrix(NaN,J,N)
  for (n in 1:N){ #Time
    for (j in 1:J){ #S mesh
      f[j,n] <- bs(x_vector[j],K,r,sigma,tau[n],option)
    }
  }
  return(f)
}


#################################################################
######################### Cut-off ###############################
#################################################################

#Define set
sv_set <- function(param, grid){
  K <- param$K
  
  s_vector <- grid$s_vector
  
  s_l <- K*(1-0.5) #Lower bound
  s_u <- K*(1+0.5) #Upper bound
  
  s_l <- min(which(s_vector > s_l)) #Index for lower bound
  s_u <- max(which(s_vector < s_u)) #Index for upper bound

  s <- s_l:s_u
  
  if(exists("v_vector", grid)){
    v_vector <- grid$v_vector
    
    v_l <- 1 #Index for lower bound
    v_u <-max(which(v_vector < 1)) #Index for upper bound
    v <- v_l:v_u
  }
  return(list(s = s, v = v))
}


#################################################################
##################### Error functions ###########################
#################################################################

err_l2 <- function(est, an){
  N <- length(est)
  err <- sqrt(sum((est-an)^2))/N
  return(err)
}
err_inf <- function(est, an){
  err <- max(abs(est-an))
  return(err)
}
err_mse <- function(est, an){
  N <- length(est)
  err <- sum((est-an)^2)/N
  return(err)
}

