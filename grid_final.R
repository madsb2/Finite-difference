map_s <- function(xi,K,c1){
  return(K+c1*sinh(xi))
}

map_v <- function(xi,d){
  return(d*sinh(xi))
}

#make_grid(50,25,800,100,100,100/5,5,0.04,5/500)

####################
#S0: Current value of S
#V0: Current value of V
make_grid <- function(param, S0, Smax, V0, Vmax, S_length, V_length){
  
  m1 <- S_length-1
  m2 <- V_length-1
  S <- Smax
  V <- Vmax
  K <- param$K
  c1 <- K/5
  d <- V/500
  
  dx <- 1.0/m1 * (asinh((S-K)/c1)-asinh(-K/c1))
  Uniform_s <- rep(0,m1+1)
  Uniform_v <- rep(0,m2+1)
  Vec_s <- rep(0,m1+1)
  Vec_v <- rep(0,m2+1)
  ds <- rep(0,m1)
  dv <- rep(0,m2)
  
  for (i in 1:(m1+1)){
    Uniform_s[i] <- asinh(-K/c1)+(i-1)*dx 
  }
  for (i in 1:(m1+1)){
    Vec_s[i] <- map_s(Uniform_s[i],K,c1)
  }
  #Vec_s <- head(sort(c(Vec_s,S0)),-1)
  #Vec_s <- sort(Vec_s)
  for (i in 1:(m1)){
    ds[i] <- Vec_s[i+1]-Vec_s[i] 
  }
  
  d_eta <- 1/m2*asinh(V/d)
  for (i in 1:(m2+1)){
    Uniform_v[i] <- (i-1)*d_eta
  }
  for (i in 1:(m2+1)){
    Vec_v[i] <- map_v(Uniform_v[i],d)
  }
  #Vec_v <- head(sort(c(Vec_v,V0)),-1)
  #Vec_v <- sort(Vec_v)
  for (i in 1:m2){
    dv[i] <- Vec_v[i+1]-Vec_v[i]
  }
  
  X <- matrix(rep(Vec_s,m2+1),nrow = m2+1,byrow = T)
  Y <- matrix(rep(Vec_v,m1+1),ncol=m1+1,byrow = F)
  
  Vec_s[1] <- 0 #negative numbers sometimes show up
  return(list("Vec_s"=Vec_s,"Vec_v"=Vec_v,"ds"=ds,"dv"=dv,"X"=X,"Y"=Y))
}

# test <- make_grid(50,25,800,100,100,100/5,5,0.04,5/500)
# library("pracma")
# qq <- meshgrid(test$Vec_s,test$Vec_v)
# par(mfrow=c(1,2))
# plot(qq$X,qq$Y,xlab = "S",ylab = "v",ylim =c(0,5),xlim=c(0,800),main="Non-uniform")
# plot(qq$X,rep(0,length(qq$Y)),xlab = "S",ylab="",yaxt="n",xlim=c(0,800),main="Non-uniform",type = "n")
# for(i in 1:length(test$Vec_s)){
#   abline(v=test$Vec_s[i],col="red")
# }
# for(i in 1:length(test$Vec_v)){
#   abline(h=test$Vec_v[i],col="red")
# }
# xMin = 0; xMax = 800; dx = xMax/50#0.125
# x_vector = seq(xMin,xMax,dx)
# v_vector <- seq(0,5,5/25)
# plot(qq$X,qq$Y,type="n",xlab = "S",ylab = "v",ylim =c(0,5),xlim=c(0,800),main="uniform")
# for(i in 1:length(x_vector)){
#   abline(v=x_vector[i],col="red")
# }
# for(i in 1:length(v_vector)){
#   abline(h=v_vector,col="red")
# }
# plot(qq$X,rep(0,length(qq$Y)),xlab = "S",ylab="",yaxt="n",xlim=c(0,800),main="Non-uniform",type = "n")
# for(i in 1:length(test$Vec_s)){
#   abline(v=test$Vec_s[i],col="red")
# }
# plot(qq$X,rep(0,length(qq$Y)),type="n",xlab = "S",ylab="",yaxt="n",xlim=c(0,800),main="uniform")
# for(i in 1:length(x_vector)){
#   abline(v=x_vector[i],col="red")
# }


#################################################################
####################### grid function ###########################
#################################################################

grid2 <- function(param, gridparam, unif = 'Y'){
  #Grid parameters
  Smin <- gridparam$Smin
  Smax <- gridparam$Smax
  J <- gridparam$J
  Vmin <- gridparam$Vmin
  Vmax <- gridparam$Vmax
  M <- gridparam$M
  Tmin <- gridparam$Tmin
  Tmax <- gridparam$Tmax
  N <- gridparam$N
  
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
    if(exists("eta", param)){V0 <- param$eta}else{V0 <- 1}
    
    temp <- make_grid(param = param, S0, Smax, V0, Vmax, S_length = J, V_length = M)
    s_vector <- temp$Vec_s
    v_vector <- temp$Vec_v
    t_vector <- seq(Tmin,Tmax,by = Tmax/(N-1))
    
    grid <- list(s_vector = s_vector, t_vector = t_vector, v_vector = v_vector)
  }
  return(grid)
}