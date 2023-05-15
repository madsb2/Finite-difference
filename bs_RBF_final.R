#Second derivative
gamma2 <- function(i, pos, ds,c){
  c <- c*ds[i]
  l <- ds[i+1]/ds[i]
  h <- ds[i] 
  if(pos==-1){
    temp <- 2/((1+l)*(h^2))
    temp2 <- (1+0.5*l*(l^2+l)*(h^2/(c^2)))
    return(temp*temp2)
  }else if(pos==0){
    return (-2/(ds[i]^2*l)*(1-(l^2-4*l+1)/(2*l) * ds[i]^2/(c^2)))
  }else if(pos==1){
    return(2/(ds[i]^2*l*(1+l)) * (1 + (3*l-1)/2 * ds[i]^2/c^2))
  }else{
    return(NA)
  }
}

#First derivative
beta2 <- function(i, pos, ds,c){
  l <- ds[i+1]/ds[i]
  c <- c*ds[i]
  if(pos==-1){
    return(-l/(ds[i]*(1+l))*(1+l*0.5*ds[i]^2/c^2))
  }else if(pos==0){
    return((l-1)/(ds[i]*l)*(1+l*0.5*ds[i]^2/c^2))
  }else if(pos==1){
    return(1/(ds[i]*l*(1+l))*(1+0.5*l*ds[i]^2/c^2))
    
  }else{
    return(NA)
  }
}

#RBF_FD approximation of Black-Scholes w. thomas
bs_rbf <- function(c, theta, param, grid, s0=param$k){
  
  startTime <- Sys.time()
  
  #Grid
  s_vector <- grid$s_vector
  t_vector <- grid$t_vector
  ds <- diff(s_vector)
  dt <- t_vector[2]-t_vector[1]
  
  #Parameters
  r <- param$r
  sigma <- param$sigma
  K <- param$K
  N <- length(s_vector)
  M <- length(t_vector)
  
  a1 <-a2 <- a3 <- b1 <- b2 <- b3 <- rep(0,N-2)
  
  #####
  #Setting up pay-off matrix with boundary conditions. The lower bound is implicitly code in as the floor is 0
  
  boundary <- rep(0,N-2)
  prices <- pmax(s_vector[2:(N-1)]-K,0) #Initial condition
  
  #Create A B C (Coefficients for FD)
  A <- sigma^2*0.5*s_vector[2:(N-1)]^2*gamma2(1:(N-2),-1,ds,c) + r*s_vector[2:(N-1)]*beta2(1:(N-2),-1,ds,c)
  B <- sigma^2*0.5*s_vector[2:(N-1)]^2*gamma2(1:(N-2),0,ds, c) + r*s_vector[2:(N-1)]*beta2(1:(N-2),0,ds, c)-r
  C <- sigma^2*0.5*s_vector[2:(N-1)]^2*gamma2(1:(N-2),1,ds, c) + r*s_vector[2:(N-1)]*beta2(1:(N-2),1,ds, c)
  
  #Diagonal
  B0 <- 1+(1-theta)*dt*B
  B1 <- -1 + theta*dt*B
  #Upper diagonal
  C0 <-(1-theta)*dt*C #We discard last element in thomas
  C1 <- theta*dt*C
  #Lower diagonal
  A0 <-(1-theta)*dt*A #We discard first element in thomas
  A1 <- theta*dt*A
  
  d <- rep(NA, N-2)
  #Boundary conditions
  for(i in 1:(M-1)){
    
    boundary[N-2] <- (s_vector[N]-K*exp(-r*t_vector[i])) *  (C0[N-2]+C1[N-2])#dt*C[(N-2)]
    
    d[1]      <- B0[1]*prices[1] + C0[1]*prices[2]
    d[N-2]    <- A0[N-2]*prices[N-3] + B0[N-2]*prices[N-2]
    d[2:(N-3)]<- A0[2:(N-3)]*prices[1:(N-4)] + B0[2:(N-3)]*prices[2:(N-3)] + C0[2:(N-3)]*prices[3:(N-2)]
    d <- d + boundary
    
    prices <- -thomas(A1, B1, C1, d)
  }
  prices <- c(0, prices, s_vector[N]-exp(-r*t_vector[M])*K)
  
  #estimate specific price
  idL <- which.min(abs(s_vector - s0))
  idU <- which.max(abs(s_vector - s0))
  
  a <- (prices[idU] - prices[idL]) / (s_vector[idU] - s_vector[idL])
  price0 <- prices[idL] + a*(s0 - s_vector[idL])
  
  endTime <- Sys.time()
  run_time <- endTime-startTime
  
  return(list("prices"=prices,"Run_time"=run_time, "price0" = price0))
}


