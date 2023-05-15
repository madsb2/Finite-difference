##################################
############## A1 ################
##################################

matA1 <- function(param, grid){
  #grid
  s_vector <- grid$s_vector
  v_vector <- grid$v_vector
  ds <- diff(s_vector)
  
  #m1 and m2
  m1 <- length(s_vector)-1
  m2 <- length(v_vector)-1
  
  #Parameters
  r <- param$r
  
  #Stacked vectors
  s_stacked <- rep(s_vector, m2+1)
  v_stacked <- rep(v_vector, each = m1+1)
  
  #Vectors
  l <- d <- u <- rep(NA, (m1+1)*(m2+1))
  
  #elements 
  i_s <- rep(c(m1+2,seq(1,m1-1),m1+2), m2) #central diff
  i_last <- seq(m1+1, (m1+1)*m2, by = m1+1) # - 0.5 in last row
  rows <- 1:((m1+1)*m2) #Rows in A1
  
  l[rows] <- 0.5*s_stacked[rows]^2*v_stacked[rows] * delta(i_s,-1, ds) + r*s_stacked[rows] * beta(i_s,-1, ds)
  d[rows] <- 0.5*s_stacked[rows]^2*v_stacked[rows] * delta(i_s, 0, ds) + r*s_stacked[rows] * beta(i_s, 0, ds) - 0.5*r
  u[rows] <- 0.5*s_stacked[rows]^2*v_stacked[rows] * delta(i_s, 1, ds) + r*s_stacked[rows] * beta(i_s, 1, ds)
  
  ds_l <- 2/(ds[m1]*(ds[m1]+ds[m1]))
  ds_d <- -2/(ds[m1]*ds[m1])+ds_l
  #ds_u <- 2/(ds[m1]*(ds[m1]+ds[m1])) Added in boundary
  
  l[i_last] <- 0.5*s_stacked[i_last]^2*v_stacked[i_last]*ds_l
  d[i_last] <- 0.5*s_stacked[i_last]^2*v_stacked[i_last]*ds_d - 0.5*r
  
  #Last row is diagonal of -0.5r due to constraints
  rows <- (((m1+1)*m2+1)+1) : ((m1+1)*(m2+1))
  d[rows] <- -0.5*r

  A1 <- list(l=l, d=d, u=u)
  A1 <- lapply(A1, function(x) replace(x, is.na(x), 0))
  
  return(A1)
}

###################################
############### A2 ################
###################################

matA2 <- function(param, grid){
  #grid
  s_vector <- grid$s_vector
  v_vector <- grid$v_vector
  ds <- diff(s_vector)
  dv <- diff(v_vector)
  
  #m1 and m2
  m1 <- length(s_vector)-1
  m2 <- length(v_vector)-1
  
  #Parameters
  r <- param$r
  theta <- param$eta
  kappa <- param$kappa
  sigma <- param$sigma
  
  #Stacked vectors
  s_stacked <- rep(s_vector, m2+1)
  v_stacked <- rep(v_vector, each = m1+1)
  
  ll <- l <- d <- u <- uu <- rep(0, (m1+1)*(m2+1))
  
  #First variance in grid (Forward)
  rows <- 2:(m1+1)#1:(m1+1)
  v_i <- 1
  d[rows]  <- kappa*(theta-v_stacked[rows]) * gamma(v_i-1,0,dv) - 0.5*r
  u[rows]  <- kappa*(theta-v_stacked[rows]) * gamma(v_i-1,1,dv)
  uu[rows] <- kappa*(theta-v_stacked[rows]) * gamma(v_i-1,2,dv)
  
  #Variance > 1 (Backward)
  if(v_vector[m2] > 1){ #Checks if the 2nd largest v > 1. Remember: Last v block row is empty

    v_i <- rep(min(which(v_vector  > 1)):m2, each = m1)
    rows <- (2:(m1+1)) + (v_i-1) * (m1+1)  #Last block row is empty
    
    ll[rows] <- kappa*(theta-v_stacked[rows]) * alpha(v_i-1,-2,dv)+ 0.5*sigma^2*v_stacked[rows]*delta(v_i-2,-1,dv)
    l[rows]  <- kappa*(theta-v_stacked[rows]) * (alpha(v_i-1,-1,dv))+ 0.5*sigma^2*v_stacked[rows]*delta(v_i-2,0,dv)
    d[rows]  <- kappa*(theta-v_stacked[rows]) * (alpha(v_i-1,0,dv)) + 0.5*sigma^2*v_stacked[rows]*delta(v_i-2,1,dv) - 0.5*r
    #u[rows]  <- kappa*(theta-v_stacked[rows])*beta(v_i-1,1,dv)                         
  }
  
  #Variances in between (Central)
  if(v_i[1] != 2){
    
    v_i <- rep(2:max(which(v_vector  <= 1)), each = m1)
    rows <- (2:(m1+1)) + (v_i-1) * (m1+1)  #Last block row is empty

    l[rows] <- kappa*(theta-v_stacked[rows]) * beta(v_i-1,-1,dv) + 0.5*sigma^2*v_stacked[rows] * delta(v_i-1,-1,dv)
    d[rows] <- kappa*(theta-v_stacked[rows]) * beta(v_i-1, 0,dv) + 0.5*sigma^2*v_stacked[rows] * delta(v_i-1, 0,dv) - 0.5*r
    u[rows] <- kappa*(theta-v_stacked[rows]) * beta(v_i-1, 1,dv) + 0.5*sigma^2*v_stacked[rows] * delta(v_i-1, 1,dv)
  }
  
  #Last block row is diagonal of -0.5r due to constraints
  rows <- (((m1+1)*m2+1)+1) : ((m1+1)*(m2+1))
  d[rows] <- -0.5*r
  
  A2 <- list(ll=ll, l=l, d=d, u=u, uu=uu)
  A2 <- lapply(A2, function(x) replace(x, is.na(x), 0))
  
  return(A2)
}

###################################
############### A0 ################
###################################

matA0 <- function(param, grid){
  #grid
  s_vector <- grid$s_vector
  v_vector <- grid$v_vector
  ds <- diff(s_vector)
  dv <- diff(v_vector)
  
  #m1 and m2
  m1 <- length(s_vector)-1
  m2 <- length(v_vector)-1
  
  #Parameters
  rho <- param$rho
  sigma <- param$sigma
  
  #Stacked vectors
  s_stacked <- rep(s_vector, m2+1)
  v_stacked <- rep(v_vector, each = m1+1)
  
  ll <- ld <- lu <- dl <- dd <- du <- ul <- ud <- uu <- rep(0, (m1+1)*(m2+1))
  
  #Elements
  #i_s <- rep(c(m1+2,rep(1,m1-1),m1+2), (m2-1))
  i_s <- rep(c(m1+2,seq(1,m1-1),m1+2), (m2-1))
  i_v <- rep(2:m2, each = m1+1)
  rows <- (m1+2):((m1+1)*m2)
  
  #Lower block matrix
  ll[rows] <- rho*sigma*s_stacked[rows]*v_stacked[rows] * beta(i_v-1,-1,dv) * beta(i_s,-1,ds)
  ld[rows] <- rho*sigma*s_stacked[rows]*v_stacked[rows] * beta(i_v-1,-1,dv) * beta(i_s, 0,ds)
  lu[rows] <- rho*sigma*s_stacked[rows]*v_stacked[rows] * beta(i_v-1,-1,dv) * beta(i_s, 1,ds)
  
  #Diagonal block matrix
  dl[rows] <- rho*sigma*s_stacked[rows]*v_stacked[rows] * beta(i_v-1,0,dv) * beta(i_s,-1,ds)
  dd[rows] <- rho*sigma*s_stacked[rows]*v_stacked[rows] * beta(i_v-1,0,dv) * beta(i_s, 0,ds)
  du[rows] <- rho*sigma*s_stacked[rows]*v_stacked[rows] * beta(i_v-1,0,dv) * beta(i_s, 1,ds)
  
  #Upper block matrix
  ul[rows] <- rho*sigma*s_stacked[rows]*v_stacked[rows] * beta(i_v-1,1,dv) * beta(i_s,-1,ds)
  ud[rows] <- rho*sigma*s_stacked[rows]*v_stacked[rows] * beta(i_v-1,1,dv) * beta(i_s, 0,ds)
  uu[rows] <- rho*sigma*s_stacked[rows]*v_stacked[rows] * beta(i_v-1,1,dv) * beta(i_s, 1,ds)
  
  A0 <- list(ll=ll, ld=ld, lu=lu, dl=dl, dd=dd, du=du, ul=ul, ud=ud, uu=uu)
  A0 <- lapply(A0, function(x) replace(x, is.na(x), 0))
  return(A0)
}

###################################
############### A #################
###################################

matA <- function(param, grid){
  startTime <- Sys.time()
  
  A0 <- matA0(param, grid)
  A1 <- matA1(param, grid)
  A2 <- matA2(param, grid)
  
  A0 <- as_mat(A0, grid)
  A1 <- as_mat(A1, grid)
  A2 <- as_mat(A2, grid)
  A <- A0+A1+A2
  
  endTime <- Sys.time()
  run_time <- endTime-startTime
  return(list('A' = A,'A0'=A0,'A1'=A1,'A2'=A2, 'run_time'=run_time))
}

###################################
############### b #################
###################################

make_boundaries <- function(param, grid){
  #grid
  s_vector <- grid$s_vector
  v_vector <- grid$v_vector
  t_vector <- grid$t_vector
  ds <- diff(s_vector)
  dv <- diff(v_vector)
  
  #Parameters
  r <- param$r
  K <- param$K
  
  #m1 and m2
  m1 <- length(s_vector)-1
  m2 <- length(v_vector)-1
  
  b0 = rep(0,(m1+1)*(m2+1))
  b1 = rep(0,(m1+1)*(m2+1))
  b2 = rep(0,(m1+1)*(m2+1))
  
  # Boundary when s = S
  ds_u <- 2/(ds[m1]*(ds[m1]+ds[m1])) #Ghost point
  
  rows <- seq(m1+1, (m1+1)*m2, m1+1)
  b1[rows] <- 0.5*v_vector[1:m2]*s_vector[m1+1]^2*ds_u*ds[m1] + r*s_vector[m1+1]
  
  # Boundary when v = V
  rows <- ((m1+1)*m2+1) : ((m1+1)*(m2+1)) #Not last row, since it is in b1
  b1[rows] <- 0.5*r*s_vector
  
  b <- b0 + b1 + b2
  
  return(list("b0" =b0,"b1" =b1, "b2"=b2,"b"= b))
}

###################################
########### as matrix #############
###################################

library(Matrix)
as_mat <- function(A, grid){
  if(length(A) == 3){#A1
    m <- length(A$d)
    
    values <- c(A$l[2:m], A$d, A$u[1:(m-1)])
    rows   <- c(2:m          , 1:m , 1:(m-1))
    cols   <- c(1:(m-1)      , 1:m , 2:m)
    
    mat <- sparseMatrix(i = rows, j = cols, x = values)
    
    return(mat)
  }
  else if(length(A) == 5){#A2
    m <- length(A$d)
    m1 <- length(grid$s_vector)-1
    
    length(c(1:(m-2*(m1+1))))
    length(A$ll[(1+2*(m1+1)):m])
    values <- c(A$ll[(1+2*(m1+1)):m], A$l[(1+1*(m1+1)):m], A$d, A$u[1:(m-1*(m1+1))], A$uu[1:(m-2*(m1+1))])
    rows   <- c((1+2*(m1+1)):m       , (1+1*(m1+1)):m      , 1:m  , 1:(m-1*(m1+1))      , 1:(m-2*(m1+1)))
    cols   <- c(1:(m-2*(m1+1))       , 1:(m-1*(m1+1))      , 1:m  , (1+1*(m1+1)):m      , (1+2*(m1+1)):m)
    
    mat <- sparseMatrix(i = rows, j = cols, x = values)
    
    return(mat)
  }
  else{ #A0
    m <- length(A$dd)
    m1 <- length(grid$s_vector)-1
    
    rows <- (m1+3):(m-(m1+2))
    cols <- 1:(m-2*(m1+1)-2)+1
    
    values <- c(A$ll[rows], A$ld[rows], A$lu[rows],
                A$dl[rows], A$dd[rows], A$du[rows],
                A$ul[rows], A$ud[rows], A$uu[rows], 0)
    rows   <- c(rep(rows, 9),m)
    cols   <- c(cols+(-1+0*(m1+1)), cols+(0+0*(m1+1)), cols+(1+0*(m1+1)),
                cols+(-1+1*(m1+1)), cols+(0+1*(m1+1)), cols+(1+1*(m1+1)),
                cols+(-1+2*(m1+1)), cols+(0+2*(m1+1)), cols+(1+2*(m1+1)), m)
    
    mat <- sparseMatrix(i = rows, j = cols, x = values)
    return(mat)
  }
}


