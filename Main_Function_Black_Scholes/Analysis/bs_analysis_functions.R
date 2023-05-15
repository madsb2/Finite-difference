################################################################
############# Accuracy at different N values ###################
################################################################

an_bs_acc_all <- function(param, grid, s0 = 90, N_val, c1 = c(NA)){
  
  #Parameters
  t_vector <- grid$t_vector
  N <- length(t_vector)
  Tmin <- t_vector[1]
  Tmax <- t_vector[N]
  
  Q <- length(N_val)
  #Analytical value
  bs_val <- rep(NA, Q)
  for(i in 1:Q){
    bs_val[i] <- bs(S = s0, param$K, param$r, param$sigma, Tmax)
  }
  
  mat <- matrix(NA, Q, 3)
  thetas <- c(0, 1, 0.5) #ex, im, cn
  
  if(is.na(c1) == TRUE){ #fd
    print('fd scheme')
    for(case in 1:3){
      print(case)
      theta <- thetas[case]
      
      for(i in 1:Q){
        grid$t_vector <- seq(Tmin,Tmax,by = Tmax/(N_val[i]-1))
        fd <- bs_fd(theta, grid, param, s0)
        res <- max(1 - abs((bs_val[i] - fd$price0)/bs_val[i]),0)
        if(is.na(res) == TRUE){ res <- 0}
        mat[i,case] <- res
      }
    }
  }
  else{ #rbf
    print('rbf scheme')
    for(case in 1:3){
      print(case)
      theta <- thetas[case]
      for(i in 1:Q){
        grid$t_vector <- seq(Tmin,Tmax,by = Tmax/(N_val[i]-1))
        fd <- bs_rbf(c1, theta, grid, param, s0)
        res <- max(1 - abs((bs_val[i] - fd$price0)/bs_val[i]),0)
        if(is.na(res) == TRUE){ res <- 0}
        mat[i,case] <- res
      }
    }
    
  }
  
  mat <- matrix(c(N_val, mat), Q, 4)
  colnames(mat) <- c('N','Ex', 'Im', 'CN')
  return(mat)
}

#an_bs_acc(param, grid, 90, N_val)

################################################################
####### Accuracy of Im and CN at different s0 values ###########
################################################################
an_bs_acc_imcn <- function(param, grid, c1 = NA){
  
  #Parameters
  t_vector <- grid$t_vector
  N <- length(t_vector)
  Tmax <- t_vector[N]
  
  #S0 values to compute prices for
  s_val <- seq(0.5*param$K, 1.5*param$K, 10)
  Q <- length(s_val)
  
  #Analytical value
  bs_val <- rep(NA, Q)
  for(i in 1:Q){
    bs_val[i] <- bs(s_val[i], param$K, param$r, param$sigma, Tmax)
  }
  
  mat <- array(NA, dim = c(Q, 3, 2))
  thetas <- c(1, 0.5) #im, cn
  
  if(is.na(c1) == TRUE){ #fd
    print('fd scheme')
    for(case in 1:2){
      theta <- thetas[case]
      for(i in 1:Q){
        fd <- bs_fd(theta, grid, param, s0 = s_val[i])
        mat[i,,case] <- c(fd$price0, max(1 - abs((bs_val[i] - fd$price0)/bs_val[i]),0), fd$Run_time)
      }
    }
  }
  else{ #rbf
    print('rbf scheme')
    for(case in 1:2){
      theta <- thetas[case]
      for(i in 1:Q){
        fd <- bs_rbf(c1, theta, grid, param, s0 = s_val[i])
        mat[i,,case] <- c(fd$price0, max(1 - abs((bs_val[i] - fd$price0)/bs_val[i]),0), fd$Run_time)
      }
    }
  }
  mat <- matrix(c(s_val, bs_val,mat[,,1], mat[,,2]), Q, 8)
  colnames(mat) <- c('s0', 'bs', 'im', 'acc', 'time', 'cn', 'acc', 'time')
  return(mat)
}

#an_bs_acc_imcn(param, grid)


