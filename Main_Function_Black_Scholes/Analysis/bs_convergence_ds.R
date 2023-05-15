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
Smin = 0; Smax = 8 * param$K; J <- 801
Vmin = 0; Vmax = 5; M <- 250
Tmin = 0; Tmax = 1; N <- 100

plot_sen_dx <- function(x, y1, y2, xlab, ylab, name1, name2, unif){
  if(unif == 'Y'){
    grid <- TeX(r'(Sensitivity of uniform $\Delta s$)')
  }else{
    grid <- TeX(r'(Sensitivity of non-Uniform $\Delta s$)')
  }
  plot(x, y1, xlab = xlab, ylab = ylab, type = 'o',
       ylim = c(min(c(y1,y2)), max(c(y1,y2))),
       main = grid); 
  
  grid()
  points(J_val,y2, type = 'o', pch = 2, lty = 2)
  legend('topright', legend=c(name1, name2), lty=1:2, pch = c(1,2), cex=0.8)
}


plot_con_dx <- function(x, y1, y2,case){
  x <- log(1/x)
  
  y1 <- log(y1) #This should be the uniform case
  y2 <- log(y2) # This should be the non-uniform case
  slope1 <- round(lm(y1~x)$coefficients[2],3)
  slope2 <- round(lm(y2~x)$coefficients[2],3)
  main1 <- paste("Slope of the uniform case: ",slope1)
  main2 <- paste("Slope of the non-uniform case: ",slope2)
  plot(x, y1, xlab = TeX(r'(Log(1/#s))'), ylab = TeX(r'($Log(L_2))'), type = 'o',
       ylim = c(min(c(y1,y2)), max(c(y1,y2))),main=case); 
  grid()
  points(x,y2, type = 'o', pch = 2, lty = 2)
  legend('topleft', legend=c(main1,main2 ), lty=1:2, pch = c(1,2), cex=0.8)
}

############### ds sensitivity ###############

J_val <- seq(50,800,50)
Q <- length(J_val)
tab <- array(NA, dim = c(Q, 3, 2))
tab2 <- array(NA, dim = c(Q, 3, 2))
for(i in 1:Q){
  J <- J_val[i]

    #Uniform
    s_vector = seq(Smin,Smax,by = Smax/(J-1))
    v_vector = seq(Vmin,Vmax,by = Vmax/(M-1))
    t_vector = seq(Tmin,Tmax,by = Tmax/(N-1))
    tau <- Tmax-t_vector
    
    grid <- list(s_vector = s_vector, t_vector = t_vector, v_vector = v_vector)
    
  #Define set
  set <- sv_set(param, grid)$s
  an <- bs_call(param = param, grid = grid, option = 'call')[,1]

  
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
  

  #Define set
  set <- sv_set(param, grid)$s
  an <- bs_call(param = param, grid = grid, option = 'call')[,1]

  
  thetas <- c(1, 0.5)
  for(case in 1:2){
    theta <- thetas[case]
    fd <- bs_fd(theta, grid, param)
    tab2[i,, case] <- c(err_l2(an[set],  fd$prices[set]),
                        err_inf(an[set], fd$prices[set]),fd$Run_time)
  }
}
#Save file
mat_u <- matrix(c(J_val,tab[,1:2,1], tab[,1:2,2]), Q, 5)
colnames(mat_u) <- c('J', 'err_l2', 'err_inf', 'err_l2', 'err_inf')
mat_n <- matrix(c(J_val,tab2[,1:2,1], tab2[,1:2,2]), Q, 5)
colnames(mat_n) <- c('J', 'err_l2', 'err_inf', 'err_l2', 'err_inf')

mat_ds_uniform_tabel <- matrix(c(J_val[1:10],tab[1:10,1,1],rep(NA,10),tab[1:10,2,1],rep(NA,10),tab[1:10,3,1],tab[1:10,1,2],rep(NA,10),tab[1:10,2,2],rep(NA,10),tab[1:10,3,2]),byrow=F,nrow=10,ncol = 11)
mat_ds_uniform_tabel[2:10,c(3,5,8,10)] <- round(100*c(1-tab[2:10,1,1]/tab[1:9,1,1],1-tab[2:10,2,1]/tab[1:9,2,1],1-tab[2:10,1,2]/tab[1:9,1,2],1-tab[2:10,2,2]/tab[1:9,2,2]),2)
colnames(mat_ds_uniform_tabel) <- c("#s","Im: L2 ","Imp. (%)","Im: L INF","Imp. (%)","RT. (s)","CN: L2","Imp. (%)","CN: L INF","Imp. (%)","RT (s)")

mat_ds_non_uniform_tabel <- matrix(c(J_val[1:10],tab2[1:10,1,1],rep(NA,10),tab2[1:10,2,1],rep(NA,10),tab2[1:10,3,1],tab2[1:10,1,2],rep(NA,10),tab2[1:10,2,2],rep(NA,10),tab2[1:10,3,2]),byrow=F,nrow=10,ncol = 11)
mat_ds_non_uniform_tabel[2:10,c(3,5,8,10)] <- round(100*c(1-tab2[2:10,1,1]/tab2[1:9,1,1],1-tab2[2:10,2,1]/tab2[1:9,2,1],1-tab2[2:10,1,2]/tab2[1:9,1,2],1-tab2[2:10,2,2]/tab2[1:9,2,2]),2)
colnames(mat_ds_non_uniform_tabel) <- c("#s","Im: L2 ","Imp. (%)","Im: L INF","Imp. (%)","RT. (s)","CN: L2","Imp. (%)","CN: L INF","Imp. (%)","RT (s)")

par(mfrow = c(1,2))
plot_con_dx(mat_u[,1],mat_u[,2],mat_n[,2],"Convergence of Implicit")
plot_con_dx(mat_u[,1],mat_u[,4],mat_n[,4],"Convergence of Cranc-Nicolson")


plot_sen_dx(mat_u[,1],mat_u[,2], mat_u[,4], '#s', TeX(r'($l_2$ error)'), 'Im', 'CN', unif = 'Y')
plot_sen_dx(mat_n[,1], mat_n[,2], mat_n[,4], '#s', TeX(r'($l_2$ error)'), 'Im', 'CN', unif = 'N')

plot_sen_dx(mat_u[,1], mat_u[,3], mat_u[,5], '#s', 'Global Error', 'Im', 'CN', unif = 'Y')
plot_sen_dx(mat_n[,1], mat_n[,3], mat_n[,5], '#s', 'Global Error', 'Im', 'CN', unif = 'N')
####################



