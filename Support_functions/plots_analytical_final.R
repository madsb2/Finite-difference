#Get path to functions
path <- getwd()
path <- strsplit(path,split="/")
if(dir.exists(paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox",sep="/"))){
  path <- paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox","Speciale","Final",sep="/")
}else{
  path <- paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox (Personlig)","Speciale","Final",sep="/")
}

#Import functions from the folder 'funktioner'
source(paste(path,"Funktioner_final.R",sep="/"))
source(paste(path,"Heston_analytical_final.R",sep="/"))

call_option <- function(s,K){
  return(max(s-K,0))
}

############################
#Define parameters
param <- list(r = 0.025, sigma = 0.3, K = 10,theta = 0.5,
              rho=-0.9,kappa=1.5,eta=0.04)

#Grid
Smin = 0; Smax = 30; J <- 50
Vmin = 0; Vmax = 2; M <- 50
Tmin = 0; Tmax = 0.5; N <- 50

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
  V0 <- param$eta
  
  temp <- make_grid(param = param, S0, Smax, V0, Vmax, S_length = J, V_length = M)
  s_vector <- temp$Vec_s
  v_vector <- temp$Vec_v
  t_vector = seq(Tmin,Tmax,by = Tmax/(N-1))
  
  grid <- list(s_vector = s_vector, t_vector = t_vector, v_vector = v_vector)
}
##############
#Option pay-off
payoff <- rep(NA,J)
for(j in 1:J){
  payoff[j] <- call_option(s_vector[j],param$K)
}


#Black Scholes call
vol <- c(0.1, 0.5, 1.2); col_vec <- c("red", "purple", "blue")
param$sigma <- vol[1]
bs_prices <- bs_call(param, grid)
plot(s_vector, bs_prices[,1], col = col_vec[1], lwd=1.5, type = 'l', xlab = 's', ylab = '', main = 'Black Scholes prices'); grid()
for(i in 2:3){
  param$sigma <- vol[i]
  bs_prices <- bs_call(param, grid)
  lines(s_vector, bs_prices[,1], col = col_vec[i], lwd=1.5)
}
lines(s_vector, payoff, lty = 'dashed')

legend("topleft", inset = 0.05, legend=c(paste("σ",vol[1]), paste("σ",vol[2]), paste("σ",vol[3]), "Pay-off"),
       col=c("red", "purple", "blue", "black"), lty=c(1,1,1,3), lwd = c(rep(1.5,3),1), cex=0.8)

#Heston call
heston_prices <- t(Heston_call(param,grid))

source(paste(path,"Funktioner_final.R",sep="/"))
FDplot1(s_vector, v_vector, heston_prices,xlabel = 's', ylabel = 'v', title = 'Heston priecs')
#close3d()

######################
#Black Scholes call
dim(heston_prices)
bs_prices <- matrix(NA, J, M)
for(m in 1:M){
  param$sigma <- v_vector[m]
  bs_prices[,m] <- bs_call(param, grid)[,1]
}

FDplot2(s_vector,v_vector,heston_prices,bs_prices,xlabel = 's', ylabel = 'v',title1 = 'Heston prices', title2 = 'Black-Scholes prices')

