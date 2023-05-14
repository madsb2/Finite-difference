#library('latex2exp')
library("plotrix") #For broken axis plot
#Figure sizes:
#550 x 500  for 1x1
#1100 x 500 for 2x1

##### Get path to functions #####
path <- getwd()
path <- strsplit(path,split="/")
if(dir.exists(paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox",sep="/"))){
  setwd(paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox","Speciale","Analysis", "bs",sep="/"))
  path1 <- paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox","Speciale","Final",sep="/")
  path2 <- paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox","Speciale","Analysis","bs",sep="/")
}else{
  setwd(paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox (Personlig)","Speciale","Analysis", "bs",sep="/"))
  path1 <- paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox (Personlig)","Speciale","Final",sep="/")
  path2 <- paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox (Personlig)","Speciale","Analysis","bs",sep="/")
}

#### Import functions from the folder 'funktioner' ####
source(paste(path1,"Funktioner_final.R",sep="/"))
source(paste(path1,"grid_final.R",sep="/"))

source(paste(path1,"bs_FD_final.R",sep="/"))

source(paste(path1,"bs_RBF_final.R",sep="/"))
source(paste(path1,"TEMP.R",sep="/")) #bs_RBF_v2 (Global?)

source(paste(path2,"bs_analysis_functions.R",sep="/"))

#################################
par(mfrow=c(1,1))

#Define parameters
param <- list(r = 0.1, sigma = 0.025, K = 100)

#Grid
Smin = 0; Smax = 8 * param$K; J <- 101
Vmin = 0; Vmax = 5; M <- 10
Tmin = 0; Tmax = 1; N <- 501
unif = 'Y'

gridparam <- list(Smin = Smin, Smax = Smax, J = J,
                  Vmin = Vmin, Vmax = Vmax,M = M,
                  Tmin = Tmin, Tmax = Tmax, N = N)
grid <- grid2(param, gridparam, unif)

########################################
############### Accuracy ###############
########################################
#Grid
Smin = 0; Smax = 8 * param$K; J <- 801 # ds = 3
Tmin = 0; Tmax = 1; N <- 501
unif = 'Y'

gridparam <- list(Smin = Smin, Smax = Smax, J = J,
                  Vmin = Vmin, Vmax = Vmax,M = M,
                  Tmin = Tmin, Tmax = Tmax, N = N)
grid <- grid2(param, gridparam, unif = unif)

print(paste('J:',J,'N:',N,'ds:',diff(grid$s_vector)[1], 'dt:', diff(grid$t_vector)[1]))

############################
## Find explicit cut-off ##
###########################
#Analytical solution
s0 <- 100
bs(S = s0, param$K, param$r, param$sigma, Tmax)

#Change grid size N
gridparam$J <- 2401
grid <- grid2(param, gridparam, unif = unif)

#Compute prices and accuracy
bs_fd(theta = 0, grid, param, s0)$price0
#1-abs(res-ex)/res

#Create plot
N_val <- seq(5, 80, 5)
mat <- an_bs_acc_all(param, grid, s0 = s0, N_val, c1 = NA) #ex, im, cn
plot(N_val, mat[,3], lty = 1, pch = 1, type ='o', ylim = c(0,1), ylab = '', xlab = '#T', main = 'Accuracy'); grid()
lines(N_val, mat[,4], lty = 2, pch = 2, type ='o')
lines(N_val, mat[,2], lty = 3, pch = 4, type ='o')
legend('right', legend = c('Im', 'CN', 'Ex'), lty = 1:3, pch = c(1,2,4))

round(100*mat,2)
#########################################
############### Im and CN ###############
#########################################
Smin = 0; Smax = 8 * param$K; J <- 2401 # ds = 3
Tmin = 0; Tmax = 1; N <- 101

gridparam <- list(Smin = Smin, Smax = Smax, J = J,
                  Vmin = Vmin, Vmax = Vmax,M = M,
                  Tmin = Tmin, Tmax = Tmax, N = N)
grid <- grid2(param, gridparam, unif = 'Y')

mat <- an_bs_acc_imcn(param, grid)
mat
#write.csv(mat, file = 'fd_accuray.csv', row.names = TRUE)
##########################

############### ds sensitivity ###############
N <- 2001
J_val <- seq(50,1000,50)
Q <- length(J_val)
tab <- array(NA, dim = c(Q, 2, 2))
for(i in 1:Q){
  J <- J_val[i]
  unif <- 'Y'
  if(unif == 'Y'){
    #Uniform
    s_vector = seq(Smin,Smax,by = Smax/(J-1))
    v_vector = seq(Vmin,Vmax,by = Vmax/(M-1))
    t_vector = seq(Tmin,Tmax,by = Tmax/(N-1))
    tau <- Tmax-t_vector
    
    grid <- list(s_vector = s_vector, t_vector = t_vector, v_vector = v_vector)
  }else{
    #Non-uniform grid
    S0 <- param$K
    V0 <- 1
    
    temp <- make_grid(param = param, S0, Smax, V0, Vmax, S_length = J, V_length = M)
    s_vector <- temp$Vec_s
    v_vector <- temp$Vec_v
    t_vector <- seq(Tmin,Tmax,by = Tmax/(N-1))
    
    grid <- list(s_vector = s_vector, t_vector = t_vector, v_vector = v_vector)
  }
  
  #Define set
  s_l <- param$K*(1-0.5) #Lower bound
  s_u <- param$K*(1+0.5) #Upper bound
  
  s_l <- find_element_index(s_vector, s_l, param$K) #Index for lower bound
  s_u <- find_element_index(s_vector, s_u, param$K) #Index for upper bound
  set <- s_l:s_u
  
  an <- bs_call(param = param, grid = grid, option = 'call')[,1]
  
  thetas <- c(1, 0.5)
  for(case in 1:2){
    theta <- thetas[case]
    fd <- bs_fd(theta, grid, param)$prices
    tab[i,, case] <- c(err_l2(an[set],  fd[set]),
                       err_inf(an[set], fd[set]))
  }
}
#Save file
mat <- matrix(c(J_val,tab[,,1], tab[,,2]), Q, 5)
colnames(mat) <- c('J', 'err_l2', 'err_inf', 'err_l2', 'err_inf')
if(unif == 'Y'){
  mat_u <- mat
  #write.csv(mat, file = 'sensitivity_ds_unif.csv', row.names = TRUE)
}else{
  mat_n <- mat
  #write.csv(mat, file = 'sensitivity_ds_not_unif.csv', row.names = TRUE)
}
#Plot
plot_sen_dx <- function(x, y1, y2, xlab, ylab, name1, name2, unif){
  if(unif == 'Y'){
    grid <- "(uniform)"
  }else{
    grid <- "(not uniform)"
  }
  plot(x, y1, xlab = xlab, ylab = ylab, type = 'o',
       ylim = c(min(c(y1,y2)), max(c(y1,y2))),
       main = paste("Sensitivity of ∆s", grid)); grid()
  points(J_val,y2, type = 'o', pch = 2, lty = 2)
  legend('topright', legend=c(name1, name2), lty=1:2, pch = rep(2,2), cex=0.8)
}
par(mfrow = c(1,2))

plot_sen_dx(mat_u[,1], mat_u[,2], mat_u[,4], 'J', expression(err[l2]), 'Im', 'CN', unif = 'Y')
plot_sen_dx(mat_n[,1], mat_n[,2], mat_n[,4], 'J', 'l2 error', 'Im', 'CN', unif = 'N')

plot_sen_dx(mat_u[,1], mat_u[,3], mat_u[,5], 'J', 'inf error', 'Im', 'CN', unif = 'Y')
plot_sen_dx(mat_n[,1], mat_n[,3], mat_n[,5], 'J', 'inf error', 'Im', 'CN', unif = 'N')
####################



