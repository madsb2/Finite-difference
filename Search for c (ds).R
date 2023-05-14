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



T_val <- 501



J_val <-seq(50,500,50)
Q <- length(J_val)


c_val2 <- seq(0.5,800,1/2)

W <- length(c_val2)


CN_l2_ds_uniform <- matrix(0,nrow=Q,ncol=W)

CN_l2_ds_n_uniform2 <- matrix(0,nrow=Q,ncol=W)

CN_l_ds_uniform <- matrix(0,nrow=Q,ncol=W)

CN_l_ds_n_uniform <- matrix(0,nrow=Q,ncol=W)

for (j in 1:1){
  #c_temp <- c_val[j]
  start <- Sys.time()
  print(c_val2[j]/1100)
  for(i in 1:Q){

    N <- 501
    J <- J_val[i] 
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

    
    temp <-bs_rbf(c_val2[j],0.5,param,grid)$prices[set]
    CN_l2_ds_uniform[i,j] <- err_l2(temp,an[set])
    CN_l_ds_uniform[i,j] <- err_inf(temp,an[set])
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
    
    
   
    
    
    CN_l2_ds_n_uniform[i,j] <- err_l2(bs_rbf(c_val2[j],0.5,param,grid)$prices[set],an[set])
    CN_l_ds_n_uniform[i,j] <- err_inf(bs_rbf(c_val2[j],0.5,param,grid)$prices[set],an[set])
    
    
    
  }
  slut <- Sys.time()
  print(slut-start)
}

result_CN_l2_uniform <- apply(CN_l2_ds_uniform,1,FUN=min)
result_CN_l_uniform <- apply(CN_l_ds_uniform,1,FUN=min)

result_CN_l_non_uniform <- apply(CN_l_ds_n_uniform,1,FUN=min)
result_CN_l2_non_uniform <- apply(CN_l2_ds_n_uniform,1,FUN=min)

mat_c_ds_error <- matrix(0,ncol = 10,nrow=5)
mat_c_ds_error[1,] <- result_CN_l2_uniform
mat_c_ds_error[2,] <- result_CN_l_uniform
mat_c_ds_error[3,] <- result_CN_l2_non_uniform
mat_c_ds_error[4,] <- result_CN_l_non_uniform
mat_c_ds_error[5,] <- J_val
row.names(mat_c_ds_error) <- c("L2 uni","L inf uni", "L2 N-uni","L inf N-uni","J")
#write.csv(mat_c_ds_error,"RBF_FD_BS_error.csv",col.names = F,row.names = T)

opt_par_non_CN_L <- opt_par_non_CN_L2 <- opt_par_uni_CN_L <- opt_par_uni_CN_L2  <- rep(0,Q)
for (i in 1:Q){
  opt_par_non_CN_L[i] <- c_val2[result_CN_l_non_uniform[i]==CN_l_ds_n_uniform[i,]]
  opt_par_non_CN_L2[i] <- c_val4[result_CN_l2_non_uniform[i]==CN_l2_ds_n_uniform[i,]]
  opt_par_uni_CN_L[i] <- c_val2[result_CN_l_uniform[i]==CN_l_ds_uniform[i,]]
  opt_par_uni_CN_L2[i] <- c_val2[result_CN_l2_uniform[i]==CN_l2_ds_uniform[i,]]
}



par(mfrow=c(1,2))

plot(J_val,opt_par_uni_CN_L,type = "o",main = "Uniform grid",ylab = "Shape parameter c",xlab = "#s",ylim = c(0,80))
points(J_val,opt_par_uni_CN_L2,type="o",pch=2,lty= 2)
legend('bottomright', legend=c(TeX(r"($L_\infty$)" ),TeX(r"($L_2$)")), lty=1:2, pch = c(1,2), cex=0.9)
grid()

plot(J_val,cc,type = "o",main ="Non-uniform grid",ylab = "Shape parameter c",xlab = "#s",ylim = c(0,1000))
points(J_val,opt_par_non_CN_L2,type="o",pch=2,lty=2)
legend('bottomright', legend=c(TeX(r"($L_\infty$)" ),TeX(r"($L_2$)")), lty=1:2, pch = c(1,2), cex=0.9)
grid()

#Get mat_u and mat_n from bs_convergence_ds

mat_c_ds <- matrix(0,ncol=10,nrow=5)
mat_c_ds[1,] <- (mat_u[1:10,4]-result_CN_l2_uniform)/mat_u[1:10,4]*100
mat_c_ds[2,] <- (mat_u[1:10,5]-result_CN_l_uniform)/mat_u[1:10,5]*100
mat_c_ds[3,] <- (mat_n[1:10,4]-result_CN_l2_non_uniform)/mat_n[1:10,4]*100
mat_c_ds[4,] <- (mat_n[1:10,5]-CN_l_ds_n_uniform2)/mat_n[1:10,5]*100
mat_c_ds[5,] <- J_val
row.names(mat_c_ds) <- c("L2 uni","L inf uni", "L2 N-uni","L inf N-uni","J")
#write.csv(mat_c_ds,"comparison_FD_RBF.csv",col.names = F,row.names = T)

mat_par_ds <- matrix(0, ncol=16, nrow=5)
mat_par_ds[1,] <- opt_par_uni_CN_L2
mat_par_ds[2,] <- opt_par_uni_CN_L
mat_par_ds[3,] <- opt_par_non_CN_L2
mat_par_ds[4,] <- opt_par_non_CN_L
mat_par_ds[5,] <- J_val
row.names(mat_par_ds) <- c("L2 uni","L inf uni", "L2 N-uni","L inf N-uni","J")
#write.csv(mat_par_ds,"optimal_c_ds.csv",col.names = F,row.names = T)



par(mfrow=c(1,2))

plot(J_val,mat_c_ds[1,],type = "o",main = TeX(r"(Improvement in $L_2$)"),ylab = "Improvement (%)",xlab = "#s",ylim = c(-1,3))
points(J_val,mat_c_ds[3,],type="o",pch=2,lty=2)
legend('topright', legend=c("uniform", "Non-uniform"), lty=1:2, pch = c(1,2), cex=0.8)
grid()

plot(J_val,mat_c_ds[2,],type = "o",main = TeX(r"(Improvement in $L_\infty$)"),ylab = "Improvement (%)",xlab = "#s",ylim = c(-1,50))
points(J_val,mat_c_ds[4,],type="o",pch=2,lty=2)
legend('topright', legend=c("Uniform", "non-uniform"), lty=1:2, pch = c(1,2), cex=0.8)
grid()

par(mfrow=c(1,1))
plot(log(1/J_val),log(result_CN_l2_uniform),type = "o",main = TeX(r"(Convergence in $L_2$)"),ylab = TeX(r"($\log(L_2)$)"),xlab = TeX(r'(Log(1/#s))'),ylim=c(-10,-1))
points(log(1/J_val),log(result_CN_l2_non_uniform2),type="o",pch=2,lty=2)
legend('bottom', legend=c("Uniform (Slope=2.019)","non-uniform (slope=2.44)"), lty=1:2, pch = c(1,2), cex=0.8)
grid()



