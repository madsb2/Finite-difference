library('latex2exp')
#Figure sizes:
#550 x 500  for 1x1
#1100 x 500 for 2x1

##### Get path to functions #####
path1 <- getwd()
path1 <- strsplit(path1,split="/")
if(dir.exists(paste(path1[[1]][1],path1[[1]][2],path1[[1]][3],"Dropbox",sep="/"))){
  setwd(paste(path1[[1]][1],path1[[1]][2],path1[[1]][3],"Dropbox","Speciale","Analysis", "bs",sep="/"))
  path1 <- paste(path1[[1]][1],path1[[1]][2],path1[[1]][3],"Dropbox","Speciale","Final",sep="/")
}else{
  setwd(paste(path1[[1]][1],path1[[1]][2],path1[[1]][3],"Dropbox (Personlig)","Speciale","Analysis", "bs",sep="/"))
  path1 <- paste(path1[[1]][1],path1[[1]][2],path1[[1]][3],"Dropbox (Personlig)","Speciale","Final",sep="/")
}

#Import functions from the folder 'funktioner'
source(paste(path1,"Funktioner_final.R",sep="/"))
source(paste(path1,"grid_final.R",sep="/"))

source(paste(path1,"bs_FD_final.R",sep="/"))

source(paste(path1,"bs_RBF_final.R",sep="/"))
#source(paste(path1,"TEMP.R",sep="/")) #bs_RBF_v2

#################################
par(mfrow=c(1,1))

#Define parameters
param <- list(r = 0.1, sigma = 0.025, K = 100,rho=-0.9,kappa=1.5,eta=0.04)

#Grid
Smin = 0; Smax = 8* param$K; J <- 801
Vmin = 0; Vmax = 5; M <- 250
Tmin = 0; Tmax = 1; N <- 1001


############### dt sensitivity ###############
N <- 1001
T_val <- seq(6,21,1)
#T_val <- c(10,20,40,80,160,320,640)
#T_val <- c(6,11,16,21,26,31,36)
#J <- 160*64
#J <- 611
#J <- 611*4
J <- 3000
Q <- length(T_val)
tab <- array(NA, dim = c(Q, 3, 2))
tab2 <- array(NA, dim = c(Q, 3, 2))
for(i in 1:Q){
  N <- T_val[i]
  
  #Uniform
  s_vector = seq(Smin,Smax,by = Smax/(J-1))
  v_vector = seq(Vmin,Vmax,by = Vmax/(M-1))
  t_vector = seq(Tmin,Tmax,by = Tmax/(N-1))
  tau <- Tmax-t_vector
  
  grid <- list(s_vector = s_vector, t_vector = t_vector, v_vector = v_vector)
  
  if(sum(s_vector==100)==1){
    print("Fejl Uniform")
    print(i)
  }

  #Define set
  set <- sv_set(param, grid)$s
  
  an <- bs_call(param = param, grid = grid, option = 'call')[,1]
  print(max(diff(t_vector))*param$sigma/max((diff(s_vector))^2)-0.5)
  thetas <- c(1, 0.5)
  for(case in 1:2){
    theta <- thetas[case]
    fd <- bs_fd(theta, grid, param)
    tab[i,, case] <- c(err_l2(an[set],  fd$prices[set]),
                       err_inf(an[set], fd$prices[set]),fd$Run_time)
  }
  
  #Non-uniform grid
  S0 <- param$K
  V0 <- 1
  
  temp <- make_grid(param = param, param$K, Smax, param$eta, Vmax, S_length = J, V_length = M)
  s_vector <- temp$Vec_s
  v_vector <- temp$Vec_v
  t_vector <- seq(Tmin,Tmax,by = Tmax/(N-1))
  
  grid <- list(s_vector = s_vector, t_vector = t_vector, v_vector = v_vector)
  an <- bs_call(param = param, grid = grid, option = 'call')[,1]
  
  
  if(sum(s_vector==100)==1){
    print("Fejl N-Uniform")
    print(i)
  }
  print(max(diff(t_vector))*param$sigma/max((diff(s_vector))^2)-0.5)
  #Define set
  set <- sv_set(param, grid)$s
  thetas <- c(1, 0.5)
  for(case in 1:2){
    theta <- thetas[case]
    fd <- bs_fd(theta, grid, param)
    tab2[i,, case] <- c(err_l2(an[set],  fd$prices[set]),
                        err_inf(an[set], fd$prices[set]),fd$Run_time)
  }
  if(i==Q){
    #Save file
    mat_u <- matrix(c(T_val,tab[,1:2,1], tab[,1:2,2]), Q, 5)
    colnames(mat_u) <- c('J', 'err_l2', 'err_inf', 'err_l2', 'err_inf')
    mat_n <- matrix(c(T_val,tab2[,1:2,1], tab2[,1:2,2]), Q, 5)
    colnames(mat_n) <- c('J', 'err_l2', 'err_inf', 'err_l2', 'err_inf')
  }
}


#Plot
plot_sen_dt <- function(x, y1, y2, xlab, ylab, name1, name2, unif){
  if(unif == 'Y'){
    grid <- TeX(r'(Sensitivity of uniform $\Delta t$)')
  }else{
    grid <- TeX(r'(Sensitivity of non-Uniform $\Delta t$)')
  }
  plot(x, y1, xlab = xlab, ylab = ylab, type = 'o',
       ylim = c(min(c(y1,y2)), max(c(y1,y2))),
       main = grid); 
  grid()
  points(x,y2, type = 'o', pch = 2, lty = 2)
  legend('topright', legend=c(name1, name2), lty=c(1,2), pch = c(1,2), cex=0.8)
}
par(mfrow = c(1,2))

plot_sen_dt(mat_u[,1],mat_u[,2], mat_u[,4], '#T', TeX(r'($l_2$ error)'), 'Im', 'CN', unif = 'Y')
plot_sen_dt(mat_n[,1], mat_n[,2], mat_n[,4], '#T', TeX(r'($l_2$ error)'), 'Im', 'CN', unif = 'N')

plot_sen_dt(mat_u[,1], mat_u[,3], mat_u[,5], '#T', 'Global Error', 'Im', 'CN', unif = 'Y')
plot_sen_dt(mat_n[,1], mat_n[,3], mat_n[,5], '#T', 'Global Error', 'Im', 'CN', unif = 'N')
####################
plot_con_dt <- function(x, y1, y2,case){
  x <- log(1/x)
  y1 <- log(y1) #This should be the uniform case
  y2 <- log(y2) # This should be the non-uniform case
  slope1 <- round(lm(y1~x)$coefficients[2],3)
  slope2 <- round(lm(y2~x)$coefficients[2],3)
  main1 <- paste("Slope of the uniform case: ",slope1)
  main2 <- paste("Slope of the non-uniform case: ",slope2)
  plot(x, y1, xlab = TeX(r'(Log(1/#T))'), ylab = TeX(r'($Log(L_2))'), type = 'o',
       ylim = c(min(c(y1,y2)), max(c(y1,y2))),main=case); 
  grid()
  points(x,y2, type="o",pch = 2, lty = 2)
  legend('topleft', legend=c(main1,main2 ), lty=c(1,2), pch = c(1,2), cex=0.8)
}
plot_con_dt(mat_u[,1],mat_u[,2],mat_n[,2],"Convergence of Implicit")
plot_con_dt(mat_u[,1],mat_u[,4],mat_n[,4],"Convergence of Cranc-Nicolson")

mat_dt_uniform_tabel <- matrix(c(T_val[1:10]-1,tab[1:10,1,1],rep(NA,10),tab[1:10,2,1],rep(NA,10),tab[1:10,3,1],tab[1:10,1,2],rep(NA,10),tab[1:10,2,2],rep(NA,10),tab[1:10,3,2]),byrow=F,nrow=10,ncol = 11)
mat_dt_uniform_tabel[2:10,c(3,5,8,10)] <- round(100*c(1-tab[2:10,1,1]/tab[1:9,1,1],1-tab[2:10,2,1]/tab[1:9,2,1],1-tab[2:10,1,2]/tab[1:9,1,2],1-tab[2:10,2,2]/tab[1:9,2,2]),2)
colnames(mat_dt_uniform_tabel) <- c("#s","Im: L2 ","%Improvement","Im: L INF","%Improvement","Run time","CN: L2","%Improvement","CN: L INF","%Improvement","Run time")

mat_dt_non_uniform_tabel <- matrix(c(T_val[1:10]-1,tab2[1:10,1,1],rep(NA,10),tab2[1:10,2,1],rep(NA,10),tab2[1:10,3,1],tab2[1:10,1,2],rep(NA,10),tab2[1:10,2,2],rep(NA,10),tab2[1:10,3,2]),byrow=F,nrow=10,ncol = 11)
mat_dt_non_uniform_tabel[2:10,c(3,5,8,10)] <- round(100*c(1-tab2[2:10,1,1]/tab2[1:9,1,1],1-tab2[2:10,2,1]/tab2[1:9,2,1],1-tab2[2:10,1,2]/tab2[1:9,1,2],1-tab2[2:10,2,2]/tab2[1:9,2,2]),2)
colnames(mat_dt_non_uniform_tabel) <- c("#s","Im: L2 ","%Improvement","Im: L INF","%Improvement","Run time","CN: L2","%Improvement","CN: L INF","%Improvement","Run time")



for(i in 1:Q){
  N <- T_val[i]
  
  #Uniform
  s_vector = seq(Smin,Smax,by = Smax/(J-1))
  v_vector = seq(Vmin,Vmax,by = Vmax/(M-1))
  t_vector = seq(Tmin,Tmax,by = Tmax/(N-1))
  tau <- Tmax-t_vector
  
  grid <- list(s_vector = s_vector, t_vector = t_vector, v_vector = v_vector)
  
  if(sum(s_vector==100)==1){
    print("Fejl Uniform")
    print(i)
  }
  
  #Define set
  set <- sv_set(param, grid)$s
  
  an <- bs_call(param = param, grid = grid, option = 'call')[,1]
  print(max(diff(t_vector))*param$sigma/max((diff(s_vector))^2)-0.5)
  thetas <- c(1, 0.5)
  for(case in 1:2){
    theta <- thetas[case]
    fd <- bs_rbf(sqrt(3)/2,theta, param,grid)
    tab[i,, case] <- c(err_l2(an[set],  fd$prices[set]),
                       err_inf(an[set], fd$prices[set]),fd$Run_time)
  }
  
  #Non-uniform grid
  S0 <- param$K
  V0 <- 1
  
  temp <- make_grid(param = param, param$K, Smax, param$eta, Vmax, S_length = J, V_length = M)
  s_vector <- temp$Vec_s
  v_vector <- temp$Vec_v
  t_vector <- seq(Tmin,Tmax,by = Tmax/(N-1))
  
  grid <- list(s_vector = s_vector, t_vector = t_vector, v_vector = v_vector)
  an <- bs_call(param = param, grid = grid, option = 'call')[,1]
  
  
  if(sum(s_vector==100)==1){
    print("Fejl N-Uniform")
    print(i)
  }
  print(max(diff(t_vector))*param$sigma/max((diff(s_vector))^2)-0.5)
  #Define set
  set <- sv_set(param, grid)$s
  thetas <- c(1, 0.5)
  for(case in 1:2){
    theta <- thetas[case]
    bs_rbf(sqrt(3)/2,theta, param,grid)
    tab2[i,, case] <- c(err_l2(an[set],  fd$prices[set]),
                        err_inf(an[set], fd$prices[set]),fd$Run_time)
  }
  if(i==Q){
    #Save file
    mat_u2 <- matrix(c(T_val,tab[,1:2,1], tab[,1:2,2]), Q, 5)
    colnames(mat_u2) <- c('J', 'err_l2', 'err_inf', 'err_l2', 'err_inf')
    mat_n2 <- matrix(c(T_val,tab2[,1:2,1], tab2[,1:2,2]), Q, 5)
    colnames(mat_n2) <- c('J', 'err_l2', 'err_inf', 'err_l2', 'err_inf')
  }
}

plot_sen_dt(mat_u2[,1],mat_u2[,2], mat_u2[,4], '#T', TeX(r'($l_2$ error)'), 'Im', 'CN', unif = 'Y')
plot_sen_dt(mat_n2[,1], mat_n2[,2], mat_n2[,4], '#T', TeX(r'($l_2$ error)'), 'Im', 'CN', unif = 'N')

plot_sen_dt(mat_u2[,1], mat_u2[,3], mat_u2[,5], '#T', 'Global Error', 'Im', 'CN', unif = 'Y')
plot_sen_dt(mat_n2[,1], mat_n2[,3], mat_n2[,5], '#T', 'Global Error', 'Im', 'CN', unif = 'N')

plot_con_dt(mat_u2[,1],mat_u2[,2],mat_n2[,2],"Convergence of Implicit")
plot_con_dt(mat_n2[,1],mat_u2[,4],mat_n2[,4],"Convergence of Cranc-Nicolson")



mat_dt_rbf_uniform_tabel <- matrix(c(T_val[1:7]-1,tab[1:7,1,1],rep(NA,7),tab[1:7,2,1],rep(NA,7),tab[1:7,3,1],tab[1:7,1,2],rep(NA,7),tab[1:7,2,2],rep(NA,7),tab[1:7,3,2]),byrow=F,nrow=7,ncol = 11)
mat_dt_rbf_uniform_tabel[2:7,c(3,5,8,10)] <- c(1-tab[2:7,1,1]/tab[1:6,1,1],1-tab[2:7,2,1]/tab[1:6,2,1],1-tab[2:7,1,2]/tab[1:6,1,2],1-tab[2:7,2,2]/tab[1:6,2,2])
colnames(mat_dt_rbf_uniform_tabel) <- c("#T","Im: L2 ","%Improvement","Im: L INF","%Improvement","Run time","CN: L2","%Improvement","CN: L INF","%Improvement","Run time")

mat_dt_rbf_non_uniform_tabel <- matrix(c(T_val[1:7]-1,tab2[1:7,1,1],rep(NA,7),tab2[1:7,2,1],rep(NA,7),tab2[1:7,3,1],tab2[1:7,1,2],rep(NA,7),tab2[1:7,2,2],rep(NA,7),tab2[1:7,3,2]),byrow=F,nrow=7,ncol = 11)
mat_dt_rbf_non_uniform_tabel[2:7,c(3,5,8,10)] <- c(1-tab2[2:7,1,1]/tab2[1:6,1,1],1-tab2[2:7,2,1]/tab2[1:6,2,1],1-tab2[2:7,1,2]/tab2[1:6,1,2],1-tab2[2:7,2,2]/tab2[1:6,2,2])
colnames(mat_dt_rbf_non_uniform_tabel) <- c("#T","Im: L2 ","%Improvement","Im: L INF","%Improvement","Run time","CN: L2","%Improvement","CN: L INF","%Improvement","Run time")

