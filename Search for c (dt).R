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
#Define parameters
param <- list(r = 0.1, sigma = 0.025, K = 100,rho=-0.9,kappa=1.5,eta=0.04)

#Grid
Smin = 0; Smax = 8* param$K; J <- 801
Vmin = 0; Vmax = 5; M <- 250
Tmin = 0; Tmax = 1; N <- 1001


#T_val <- c(6,11,16,21,26,31,36)
T_val <- seq(6,21,1)

#J <- 160*64
J <- 3000
Q <- length(T_val)
tab <- array(NA, dim = c(Q, 3, 2))
tab2 <- array(NA, dim = c(Q, 3, 2))


c_val2 <- seq(0.5,400,1/2)
#c_val <- seq(0.5,1000,1/3)
W <- length(c_val2)


CN_l2_c_uniform_dt <- matrix(0,nrow=Q,ncol=W)

CN_l2_c_n_uniform_dt <- matrix(0,nrow=Q,ncol=W)


CN_l_c_uniform_dt <- matrix(0,nrow=Q,ncol=W)

CN_l_c_n_uniform_dt <- matrix(0,nrow=Q,ncol=W)

for (j in 1:W){
  c_temp <- c_val2[j]
  print(c_val2[j])
  for(i in 1:Q){
    N <- T_val[i]
    # 
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

    temp <-bs_rbf(c_temp,0.5,param,grid)$prices[set]
    CN_l2_c_uniform_dt[i,j] <- err_l2(temp,an[set])
    CN_l_c_uniform_dt[i,j] <- err_inf(temp,an[set])
    #print(max(diff(t_vector))*param$sigma/max((diff(s_vector))^2)-0.5)
  
    
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
    #print(max(diff(t_vector))*param$sigma/max((diff(s_vector))^2)-0.5)
    #Define set
    set <- sv_set(param, grid)$s
  
    
    #Save file
    
    
      
    
    
    CN_l2_c_n_uniform_dt[i,j] <- err_l2(bs_rbf(c_val2[j],0.5,param,grid)$prices[set],an[set])
    CN_l_c_n_uniform_dt[i,j] <- err_inf(bs_rbf(c_val2[j],0.5,param,grid)$prices[set],an[set])
    
    
  }
}



#write.csv(CN_l2_c_uniform,"CN_l2_c_uniform2.csv",row.names = F)

#write.csv(CN_l_c_uniform,"CN_l_c_uniform2.csv",row.names = F)


#write.csv(CN_l2_c_n_uniform,"CN_l2_c_n_uniform2.csv",row.names = F)

#write.csv(CN_l_c_n_uniform,"CN_l_c_n_uniform2.csv",row.names = F)

L2_CN <- apply(CN_l2_c_uniform_dt,1,FUN = min)
L2_CN_n <- apply(CN_l2_c_n_uniform_dt,1, FUN=min)
L_inf_CN <- apply(CN_l_c_uniform_dt,1,FUN= min)


L_inf_CN_n <- apply(CN_l_c_n_uniform_dt,1,FUN=min)

#set_c <- c_val<2499.5
opt_c_uni_CN_L2_dt <- opt_c_uni_CN_inf_dt <- opt_c_non_CN_l2_dt <- opt_c_non_CN_L_inf_dt <- rep(0,Q)
for (i in 1:Q){
  opt_c_uni_CN_L2_dt[i] <- c_val2[L2_CN[i]==CN_l2_c_uniform_dt[i,]]
  opt_c_uni_CN_inf_dt[i] <- c_val2[L_inf_CN[i]==CN_l_c_uniform_dt[i,]]
  opt_c_non_CN_l2_dt[i] <- c_val2[L2_CN_n[i]==CN_l2_c_n_uniform_dt[i,]]
  opt_c_non_CN_L_inf_dt[i] <- c_val2[L_inf_CN_n[i]==CN_l_c_n_uniform_dt[i,]]
}

par(mfrow=c(1,2))

#Get mat_u and mat_n matrices from bs_convergence_dt

plot(T_val,opt_c_non_CN_l2_dt,type = "o",main = TeX(r"(Optimal shape parameter $L_2$)"),ylab = "Shape parameter c",xlab = "#s",ylim = c(0,400))
points(T_val,opt_c_uni_CN_L2_dt,type="o",pch=2,lty=2)
legend('bottomright', legend=c("Non-uniform", "Uniform"), lty=1:2, pch = c(1,2), cex=0.8)
grid()

plot(T_val,opt_c_non_CN_L_inf_dt,type = "o",main = TeX(r"(Optimal shape parameter $L_\infty$)"),ylab = "Shape parameter c",xlab = "#s",ylim = c(0,400))
points(T_val,opt_c_uni_CN_inf_dt,type="o",pch=2,lty=2)
legend('bottomright', legend=c("Non-uniform", "Uniform"), lty=1:2, pch = c(1,2), cex=0.8)
grid()

par(mfrow=c(1,2))

plot(T_val,(mat_u[,4]-L2_CN)/mat_u[,4]*100,type = "o",main = TeX(r"(Improvement in $L_2$)"),ylab = "Improvement (%)",xlab = "#T",ylim = c(-1,3))
points(T_val,(mat_n[,4]-L2_CN_n)/mat_n[,4]*100,type="o",pch=2,lty=2)
legend('topright', legend=c("uniform", "Non-uniform"), lty=1:2, pch = c(1,2), cex=0.8)
grid()

plot(T_val,(mat_u[,5]-L_inf_CN)/mat_u[,5]*100,type = "o",main = TeX(r"(Improvement in $L_\infty$)"),ylab = "Improvement (%)",xlab = "#s",ylim = c(-1,50))
points(T_val,(mat_n[,5]-L_inf_CN_n)/mat_n[,5]*100,type="o",pch=2,lty=2)
legend('topright', legend=c("Uniform", "non-uniform"), lty=1:2, pch = c(1,2), cex=0.8)
grid()

par(mfrow=c(1,1))
plot(log(1/T_val),log(L2_CN),type = "o",main = TeX(r"(Convergence in time ($L_2$))"),ylab = TeX(r"($\log(L_2)$)"),xlab = TeX(r'(Log(1/#s))'),ylim = c(-10,-6.5))
points(log(1/T_val),log(L2_CN_n),type="o",pch=2,lty=2)
legend('bottom', legend=c("Uniform (Slope=2.019)","Non-Uniform (Slope=2.148)"), lty=1:2, pch = c(1,2), cex=0.8)
grid()


