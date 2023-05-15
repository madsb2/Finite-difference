

Douglass <- function(param, grid){
  startTime <- Sys.time()
  
  A <- matA(param, grid)
  #Grid
  s_vector <- grid$s_vector
  t_vector <- grid$t_vector
  v_vector <- grid$v_vector

  dt <- t_vector[2]-t_vector[1]
  
  N <- length(t_vector)
  m1 <- length(s_vector)-1
  m2 <- length(v_vector)-1
  
  theta <- param$theta
  K <- param$K
  
  UU = pmax(s_vector-K,0)
  U = rep(UU,m2+1)
  
  A0 <- A$A0
  A1 <- A$A1
  A2 <- A$A2
  A  <- A$A
  
  tempb <- make_boundaries(param, grid)
  b <- tempb$b
  b0 <- tempb$b0
  b1 <- tempb$b1
  b2 <- tempb$b2
  #############
  
  
  I <- diag((m1+1)*(m2+1))

  l1 <- I-theta*dt*A1

  inv_l1 <- solve(l1)
  l2 <- I-theta*dt*A2
  inv_l2 <- solve(l2)
  
  for (n in 1:(N-1)){

    Y0 <- U + dt*(A%*%U+b)

    rhs1 <- Y0 - theta*dt*A1%*%U
    
    Y1 <- inv_l1%*%rhs1
    
    rhs2 <- Y1 - theta*dt*A2%*%U
    
    #print(rhs2[1:10])
    
    U <- inv_l2%*%rhs2
  }
  
  U <- matrix(U,nrow = m2+1,ncol = m1+1,byrow = T)
  
  endTime <- Sys.time()
  run_time <- endTime-startTime
  return(list("prices"=U,"run_time"=run_time))
}


Craig_sneyd <- function(param,grid){
  startTime <- Sys.time()
  
  A <- matA(param, grid)
  #Grid
  s_vector <- grid$s_vector
  t_vector <- grid$t_vector
  v_vector <- grid$v_vector
  
  ds <- diff(s_vector)
  dv <- diff(v_vector)
  
  dt <- t_vector[2]-t_vector[1]
  
  N <- length(t_vector)
  m1 <- length(s_vector)-1
  m2 <- length(v_vector)-1
  
  rho <- param$rho
  sigma <- param$sigma
  r <- param$r
  kappa <- param$kappa
  theta <- param$theta
  K <- param$K
  
  m <- (m2+1)*(m1+1)
  UU = pmax(s_vector-K,0)
  U = rep(UU,m2+1)
  
  A0 <- A$A0
  A1 <- A$A1
  A2 <- A$A2
  A  <- A$A
  
  tempb <- make_boundaries(param, grid)
  b <- tempb$b
  b0 <- tempb$b0
  b1 <- tempb$b1
  b2 <- tempb$b2
  
  I <- diag(m)
  l1 <- I-theta*dt*A1
  inv_l1 <- solve(l1)
  l2 <- I-theta*dt*A2
  inv_l2 <- solve(l2)
  for (n in 1:(N-1)){
    Y0 <- U + dt*A%*%U+b
    rhs1 <- Y0 - theta*dt*A1%*%U
    Y1 <- inv_l1%*%rhs1
    rhs2 <- Y1 - theta*dt*A2%*%U
    Y2 <- inv_l2%*%rhs2
    
    Y0_tilde <- Y0 + 0.5*dt*(A0%*%Y2-A0%*%U)
    
    rhs1 <- Y0_tilde -theta*dt*A1%*%U
    Y1_tilde <- inv_l1%*%rhs1
    rhs2 <- Y1_tilde - theta*dt*A2%*%U
    U <- inv_l2%*%rhs2
  }
  U <- matrix(U,nrow = m2+1,ncol = m1+1,byrow = T)
  
  endTime <- Sys.time()
  run_time <- endTime-startTime
  return(list("prices"=U,"run_time"=run_time))
}


Modified_Craig_sneyd <- function(param,grid){
  startTime <- Sys.time()
  
  A <- matA(param, grid)
  #Grid
  s_vector <- grid$s_vector
  t_vector <- grid$t_vector
  v_vector <- grid$v_vector
  
  ds <- diff(s_vector)
  dv <- diff(v_vector)
  dt <- t_vector[2]-t_vector[1]
  
  N <- length(t_vector)
  m1 <- length(s_vector)-1
  m2 <- length(v_vector)-1
  
  rho <- param$rho
  sigma <- param$sigma
  r <- param$r
  kappa <- param$kappa
  theta <- param$theta
  K <- param$K
  m <- (m2+1)*(m1+1)
  
  UU = pmax(s_vector-K,0)
  U = rep(UU,m2+1)
  
  A0 <- A$A0
  A1 <- A$A1
  A2 <- A$A2
  A  <- A$A
  
  tempb <- make_boundaries(param, grid)
  b <- tempb$b
  b0 <- tempb$b0
  b1 <- tempb$b1
  b2 <- tempb$b2
  
  I <- diag(m)
  l1 <- I-theta*dt*A1
  inv_l1 <- solve(l1)
  l2 <- I-theta*dt*A2
  inv_l2 <- solve(l2)
  for (n in 1:(N-1)){
    Y0 <- U + dt*(A%*%U+b)
    rhs1 <- Y0 - theta*dt*A1%*%U
    Y1 <- inv_l1%*%rhs1
    rhs2 <- Y1 - theta*dt*A2%*%U
    Y2 <- inv_l2%*%rhs2
    
    Y0_hat <- Y0 +theta*dt*(A0%*%Y2-A0%*%U)
    
    Y0_tilde <- Y0_hat + (0.5-theta)*dt*(A%*%Y2-A%*%U)
    
    rhs1 <- Y0_tilde - theta*dt*A1%*%U
    Y1_tilde <- inv_l1%*%rhs1
    rhs2 <- Y1_tilde - theta*dt*A2%*%U
    U <- inv_l2%*%rhs2
  }
  U <- matrix(U,nrow = m2+1,ncol = m1+1,byrow = T)
  
  endTime <- Sys.time()
  run_time <- endTime-startTime
  return(list("prices"=U,"run_time"=run_time))
}