##### Get path to functions #####
path <- getwd()
path <- strsplit(path,split="/")
if(dir.exists(paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox",sep="/"))){
  path0 <- paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox","Speciale","Analysis", "heston",sep="/")
  path1 <- paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox","Speciale","Final",sep="/")
  path2 <- paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox","Speciale","Visual",sep="/")
}else{
  path0 <- paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox (Personlig)","Speciale","Analysis", "heston",sep="/")
  path1 <- paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox (Personlig)","Speciale","Final",sep="/")
  path2 <- paste(path[[1]][1],path[[1]][2],path[[1]][3],"Dropbox (Personlig)","Speciale","Visual",sep="/")
}

#### Import functions ####
source(paste(path1,"grid_final.R",sep="/"))
source(paste(path1,"plot_final.R",sep="/"))
source(paste(path1,"Funktioner_final.R",sep="/"))
source(paste(path1,"Heston_analytical_final.R",sep="/"))

source(paste(path1,"coefficienter_final.R",sep="/"))
source(paste(path1,"matrix_efficient0.R",sep="/"))
source(paste(path1,"ADI_schemes_global.R",sep="/"))
source(paste(path0,"heston_analysis_functions.R",sep="/"))

##########################
tab_impr <- function(vec){
  return(-c(0, diff(vec) / head(vec, -1))*100)
}

##############################

#Define parameters
param1 <- list(kappa = 1.5, eta = 0.04, sigma = 0.3,
               rho = -0.9, r = 0.025, K = 100, Tmax = 1, theta = 0.5)
param2 <- list(kappa = 3, eta = 0.12, sigma = 0.04,
               rho = 0.6, r = 0.01, K = 100, Tmax = 1, theta = 0.5)
param3 <- list(kappa = 0.6067, eta = 0.0707, sigma = 0.2928,
               rho = -0.7571, r = 0.03, K = 100, Tmax = 3, theta = 0.5)
param4 <- list(kappa = 2.5, eta = 0.06, sigma = 0.5,
               rho = -0.1, r = 0.0507, K = 100, Tmax = 0.25, theta = 0.5)
cases <- list(param1 = param1, param2 = param2, param3 = param3, param4 = param4)

#Grid
Smin = 0; Smax = 800; J <- NA
Vmin = 0; Vmax = 5; M <- NA
Tmin = 0; Tmax = 1; N <- 101

gridparam <- list(Smin = Smin, Smax = Smax, J = J,
                  Vmin = Vmin, Vmax = Vmax,M = M,
                  Tmin = Tmin, Tmax = Tmax, N = N)

#########################################
############ Price accuracy #############
#########################################
wd <- paste(path0, "s800v5", sep="/"); setwd(wd)
ADI <- c('DO', 'CS', 'MCS')
unif <- 'Y'
gridparam$N <- 41
gridparam$J <- 51
gridparam$M <- 51

#Compute prices
run <- 'N'
if(run == 'Y'){
  for(cas in 1:4){
    for(adi in 1:3){
      an_heston_prices(cases, gridparam, m1=gridparam$J-1, m2=gridparam$M-1, unif, ADI[adi], cas, wd)
    }
  }
}

mat <- array(NA, dim = c(4, 2, 4)) #Method x price-error x cases
colnames(mat) <- c('Price', 'Accuracy'); row.names(mat) <- c('HES',ADI)

s0 <- 90
v0 <- 0.3
for(cas in 1:4){#Cases
  gridparam$Tmax <- cases[[cas]]$Tmax
  grid <- grid2(cases[[cas]], gridparam, unif)
  mat[1,1,cas] <- an_heston_price0_sol(cases[[cas]], gridparam, s0, v0)
  for(adi in 1:3){#ADI
    fd <- readRDS(paste(ADI[adi],'_case',cas,'_J',gridparam$J,'_M',gridparam$M,'_N',gridparam$N,'_unif',unif, sep = ""))$prices
    mat[adi+1,1,cas] <- an_heston_price0_fd(grid, fd, s0, v0) #Price
    mat[adi+1,2,cas] <- 1-abs(mat[1,1,cas] - mat[adi+1,1,cas])/mat[1,1,cas]
  }
}
round(mat,4)

setwd(path0)
#write.csv(mat,paste('acc_mat_unif',unif,'.csv',sep = ''))

##################################
########### M-V 3D error #########
##################################
wd <- paste(path0, "s800v5", sep="/")
case <- 1
ADI <- c('DO', 'CS', 'MCS')
unif <- 'N'

#Grid
Smin = 0; Smax = 800; J <- NA
Vmin = 0; Vmax = 5; M <- NA
Tmin = 0; Tmax = 1; N <- 101

gridparam <- list(Smin = Smin, Smax = Smax, J = J,
                  Vmin = Vmin, Vmax = Vmax,M = M,
                  Tmin = Tmin, Tmax = Tmax, N = N)

m1 <- seq(10, 100, 10) #J-1
m2 <- seq(10, 100, 10) #M-1
#m2 <- c(5,seq(10, 100, 10))#seq(10, 100, 10)

m1_vec <- rep(m1, length(m2))
m2_vec <- rep(m2, each = length(m1))
#matrix(c(m1_vec,m2_vec), ncol=2)
for(adi in 1:3){
  an_heston_prices(cases, gridparam, m1_vec, m2_vec, unif, ADI[adi], case, wd)
}

z <- an_heston_errors(cases, gridparam, m1_vec, m2_vec, unif, ADI[1], case, wd)

z_reg  <- array(z$region, dim = c(length(m1),length(m2),2)) #m1 x m2 x err
#z_glob <- array(z$global, dim = c(length(m1),length(m2),2))#1: l2, 2:inf

View(z_reg[,,2])

heston_err_plot(m1, m2, z_reg[,,1], title = bquote(L[2] ~ error ~ (.(ADI[1]))), xlabel = '#s', ylabel = '#v')
heston_err_plot(m1, m2, z_reg[,,2], title = bquote(L['\u221E'] ~ error ~ (.(ADI[1]))), xlabel = '#s', ylabel = '#v')

savepic <- 'N'
if(savepic == 'Y'){
  for(adi in 1:3){
    z <- an_heston_errors(cases, gridparam, m1_vec, m2_vec, unif, ADI[adi], case, wd)
    z_reg  <- array(z$region, dim = c(length(m1),length(m2),2))
    
    w <- paste(path0, "plots", sep="/"); setwd(w)
    heston_err_plot(m1, m2, z_reg[,,2], title = bquote(L['\u221E'] ~ error ~ (.(ADI[adi]))), xlabel = '#s', ylabel = '#v')
    rgl.snapshot(paste('Case',case,'Linf error',ADI[adi],'.png', sep = '')); close3d()
    
    heston_err_plot(m1, m2, z_reg[,,1], title = bquote(L[2] ~ error ~ (.(ADI[adi]))), xlabel = '#s', ylabel = '#v')
    rgl.snapshot(paste('Case',case,'L2 error',ADI[adi],'.png', sep = '')); close3d()
  }
}

#########################################
####### Analysis of a single case #######
#########################################
wd <- paste(path0, "s800v5", sep="/")
setwd(wd)
case <- 1
ADI <- 'DO'
unif <- 'N'
gridparam$N <- 101
gridparam$J <- 101
gridparam$M <- 51
grid <- grid2(cases[[case]], gridparam, unif)

temp_hes <- readRDS(paste('HES_case',case,'_J',gridparam$J,'_M',gridparam$M,'_N',gridparam$N,'_unif',unif, sep = ""))
temp_adi <- readRDS(paste(ADI,'_case',case,'_J',gridparam$J,'_M',gridparam$M,'_N',gridparam$N,'_unif',unif, sep = ""))

an_heston_price0_sol(cases[[case]], gridparam, s0 = 90, v0 = 0.3)
an_heston_price0_fd(grid, temp_hes$prices, s0 = 90, v0 = 0.3)


set <- sv_set(cases[[case]],grid)

#Error in our set
heston_err_plot(grid$v_vector[set$v], grid$s_vector[set$s], abs(temp_hes$prices[set$v,set$s] - temp_adi$prices[set$v,set$s]),
                xlabel = 'v', ylabel = 's', title = 'Absolute error')

#Error on whole pricing surface
heston_err_plot(grid$v_vector, grid$s_vector, abs(temp_hes$prices - temp_adi$prices),
                xlabel = 'v', ylabel = 's', title = 'Absolute error')

#rgl.snapshot('abs error.png'); close3d()
close3d()
getwd()

#########################################
############## convergence ##############
#########################################
wd <- paste(path0, "s800v5", sep="/")
unif <- 'N'
case <- 1
gridparam$N <- 501#501

m2 <- seq(10,100,10)
m1 <- 2*m2

#Compute prices
run <- 'N'
if(run == 'Y'){
  for(cas in 1:4){
    for(adi in 1:3){
      an_heston_prices(cases, gridparam, m1, m2, unif, ADI[adi], cas, wd)
    }
  }
}


#Create table for 2x2 plot
ADI <- c('DO', 'CS', 'MCS')
mat <- array(NA, dim = c(length(m2), 6, 4)) #
colnames(mat) <- c('l2', 'inf', 'l2', 'inf', 'l2', 'inf')
row.names(mat) <- m2
for(cas in 1:4){
  for(adi in 1:3){
    for(i in 1:length(m2)){
      z <- an_heston_errors(cases, gridparam, m1=m2[i]*2, m2=m2[i], unif, ADI[adi], cas, wd)
      mat[i, 1:2+2*(adi-1), cas] <- z$region[,1:2]
    }
  }
}

#Plot 2x2
error <- 2 #1 is l_2, 2 is inf

par(mfrow = c(2,2)); par(mar = c(4.8, 4.1, 3, 1))
ylabs <- c(expression(log(L[2])), expression(log(L['\u221E'])))
slope <- rep(NA, 3)

for(case in 1:4){
  plot(log(1/m2), log(mat[,error, case]), ylab = ylabs[error],xlab = 'log(1/#v)',
       main = paste('Case',case), type = 'o', lty = 1, pch = 1,
       ylim = c(min(log(mat[,c(error, error+2, error+4), case])),
                max(log(mat[,c(error, error+2, error+4), case])) ));grid()
  
  fit <- lm(log(mat[,error, case]) ~ log(1/m2))
  slope[1] <- round(fit$coefficients[2],4)
  
  for(i in 2:3){ #CS and MCS
    lines(log(1/m2), log(mat[,(i-1)*2+error, case]), type = 'o', lty = i, pch = i)
    fit <- lm(log(mat[,(i-1)*2+error, case]) ~ log(1/m2))
    slope[i] <- round(fit$coefficients[2],4)
  }
  legend('topleft',legend = c(paste(ADI[1],slope[1]),
                              paste(ADI[2],slope[2]),
                              paste(ADI[3],slope[3])),
         lty = 1:3, pch = c(1,2,4), cex = 0.8)
}

#########################################
######### Computational cost ############
#########################################
wd <- paste(path0, "s800v5", sep="/")
ADI = 'DO'
unif <- 'N'

#Chose grid size
gridparam$N <- 201
m2 <- seq(10,100,10)
m1 <- 2*m2

#Create table
mat <- array(NA, dim = c(length(m2), 15, 4)) #m2 x ADI schemes x cases
colnames(mat) <- rep(c('L_2', 'Imp.', '$inf', 'Imp.', 'RT.'),3)
row.names(mat) <- m2

schemes <- c('DO', 'CS', 'MCS')
for(cas in 1:4){
  for(adi in 1:3){
    for(i in 1:length(m2)){
      z <- an_heston_errors(cases, gridparam,m1=2*m2[i],m2=m2[i],unif,schemes[adi],cas,wd)
      mat[i, 1+5*(adi-1), cas] <- z$region[,1] #l2
      mat[i, 3+5*(adi-1), cas] <- z$region[,2] #inf
      mat[i, 5+5*(adi-1), cas] <- z$region[,3] #time - beware min and sec together
      
    }
    mat[, 2+5*(adi-1), cas] <- tab_impr(mat[, 1+5*(adi-1), cas]) #Impr l2
    mat[, 4+5*(adi-1), cas] <- tab_impr(mat[, 3+5*(adi-1), cas]) #Impr inf
  }
}
savemat <- 'N'
if(savemat == 'Y'){
  setwd(path0)
  for(cas in 1:4){
    write.csv(mat[,,cas],paste('comp_cost_ds_',cas,'.csv',sep = ''))
  }
}

mat

#Plot
par(mfrow=c(1,2));#par(mar = c(4.8, 4.1, 3, 8), xpd = TRUE)
ylabs2 <- c(expression(L[2]), expression(L['\u221E']))

case <- 4
m2 <- seq(30,100,10)
m1 <- 2*m2

for(error in 1:2){ #Inf, then l2
  ylimit <- c(min(mat[m2/10,c(1+(error-1)*2, 6+(error-1)*2, 11+(error-1)*2),case]),
              max(mat[m2/10,c(1+(error-1)*2, 6+(error-1)*2, 11+(error-1)*2),case]))
  plot( m2, mat[m2/10,1+(error-1)*2,case],  type = 'o', lty = 1, pch = 1, xlab = '#v', ylab = ylabs2[error], main = paste('Case',case), ylim = ylimit);grid() #DO
  lines(m2, mat[m2/10,6+(error-1)*2,case],  type = 'o', lty = 2, pch = 2, xlab = '#v', ylab = ylabs2[error]) #CS
  lines(m2, mat[m2/10,11+(error-1)*2,case], type = 'o', lty = 3, pch = 4, xlab = '#v', ylab = ylabs2[error]) #MCS
  
  legend('topright',legend = c('DO','CS', 'MCS'),lty = seq(1,3,1), pch = c(1,2,4), cex = 0.8)
}


#########################################
########## Sensitivity of dt ############
#########################################

wd <- paste(path0, "s800v5", sep="/")
case <- 1
unif <- 'N'

#Grid
Smin = 0; Smax = 800; J <- 101
Vmin = 0; Vmax = 5; M <- 51
Tmin = 0; Tmax = 1

gridparam <- list(Smin = Smin, Smax = Smax, J = J,
                  Vmin = Vmin, Vmax = Vmax,M = M,
                  Tmin = Tmin, Tmax = Tmax, N = N)

#Compute prices
N_val <- seq(5, 15, 2)
gridparam$J <- 201
gridparam$M <- 101

ADI <- c('DO', 'CS', 'MCS')
run <- 'N'
if(run == 'Y'){
  for(cas in 1:4){
    for(i in 1:3){ #ADI
      an_heston_dt(cases, gridparam, N_val, unif, ADI[i], cas, wd)
    }
  }
}

#Create table
mat <- array(NA, dim=c(length(N_val), 15, 4)) # N x adi x case
colnames(mat) <- c(rep(c('L2', 'Imp.', 'inf', 'Imp.', 'RT'),3)); row.names(mat) <- N_val
for(cas in 1:4){#case
  for(adi in 1:3){#adi
      z <- an_heston_errors_dt(cases, gridparam, N_val, unif, ADI[adi], cas, wd)
      mat[, 1+(adi-1)*5,cas] <- z$region[,1] #l2
      mat[, 3+(adi-1)*5,cas] <- z$region[,2] #inf
      mat[, 5+(adi-1)*5,cas] <- z$region[,3] #RT
      
      mat[, 2+(adi-1)*5,cas] <- tab_impr(mat[, 1+(adi-1)*5,cas]) #Imp. l2
      mat[, 4+(adi-1)*5,cas] <- tab_impr(mat[, 3+(adi-1)*5,cas]) #Imp. inf
  }
}
mat
savemat <- 'N'
if(savemat == 'Y'){
  setwd(path0)
  for(cas in 1:4){
    write.csv(mat[,,cas],paste('comp_cost_dt_',cas,'.csv',sep = ''))
  }
}

#Plot
par(mfrow = c(2,2)); par(mar = c(4.8, 4.1, 3, 1))
ylabs <- c(expression(log(L[2])), expression(log(L['\u221E'])))
error <- 1
slope <- rep(NA, 3)
for(cas in 1:4){ #Inf, then l2
  ylimit <- c(min(log(mat[,c(1+(error-1)*2, 6+(error-1)*2, 11+(error-1)*2),cas])),
              max(log(mat[,c(1+(error-1)*2, 6+(error-1)*2, 11+(error-1)*2),cas]))+3.5)

  #Log
  plot( log(1/N_val), log(mat[,1+(error-1)*2,cas]),  type = 'o', lty = 1, pch = 1, xlab = 'log(1/#T)', ylab = ylabs[error], main = paste('Case',cas), ylim = ylimit);grid() #DO
  lines(log(1/N_val), log(mat[,6+(error-1)*2,cas]),  type = 'o', lty = 2, pch = 2) #CS
  lines(log(1/N_val), log(mat[,11+(error-1)*2,cas]), type = 'o', lty = 3, pch = 4) #MCS
  
  fit <- lm(log(mat[,1+(error-1)*2,cas]) ~ log(1/N_val))
  slope[1] <- round(fit$coefficients[2],4)
  fit <- lm(log(mat[,6+(error-1)*2,cas]) ~ log(1/N_val))
  slope[2] <- round(fit$coefficients[2],4)
  fit <- lm(log(mat[,11+(error-1)*2,cas]) ~ log(1/N_val))
  slope[3] <- round(fit$coefficients[2],4)
  
  legend('topleft',legend = c(paste(ADI[1],slope[1]),
                              paste(ADI[2],slope[2]),
                              paste(ADI[3],slope[3])),
         lty = 1:3, pch = c(1,2,4), cex = 0.8)
}


for(cas in 1:4){ #Inf, then l2
  ylimit <- c(min(mat[,c(1+(error-1)*2, 6+(error-1)*2, 11+(error-1)*2),cas]),
              max(mat[,c(1+(error-1)*2, 6+(error-1)*2, 11+(error-1)*2),cas]))
  
  #Log
  plot( 1/N_val, mat[,1+(error-1)*2,cas],  type = 'o', lty = 1, pch = 1, xlab = 'log(1/#T)', ylab = ylabs[error], main = paste('Case',cas), ylim = ylimit);grid() #DO
  lines(1/N_val, mat[,6+(error-1)*2,cas],  type = 'o', lty = 2, pch = 2) #CS
  lines(1/N_val, mat[,11+(error-1)*2,cas], type = 'o', lty = 3, pch = 4) #MCS
  
  legend('topleft',legend = ADI,
         lty = 1:3, pch = c(1,2,4), cex = 0.8)
}



#################
N_val <- seq(5, 50, 5)

#Create table
mat <- array(NA, dim=c(length(N_val), 15, 4)) # N x adi x case
colnames(mat) <- c(rep(c('L2', 'Imp.', 'inf', 'Imp.', 'RT'),3)); row.names(mat) <- N_val
for(cas in 1:4){#case
  for(adi in 1:3){#adi
    z <- an_heston_errors_dt(cases, gridparam, N_val, unif, ADI[adi], cas, wd)
    mat[, 1+(adi-1)*5,cas] <- z$region[,1] #l2
    mat[, 3+(adi-1)*5,cas] <- z$region[,2] #inf
    mat[, 5+(adi-1)*5,cas] <- z$region[,3] #RT
    
    mat[, 2+(adi-1)*5,cas] <- tab_impr(mat[, 1+(adi-1)*5,cas]) #Imp. l2
    mat[, 4+(adi-1)*5,cas] <- tab_impr(mat[, 3+(adi-1)*5,cas]) #Imp. inf
  }
}

#Plot
case <- 1

par(mfrow=c(1,2));#par(mar = c(4.8, 4.1, 3, 8), xpd = TRUE)
ylabs2 <- c(expression(L[2]), expression(L['\u221E']))

for(error in 1:2){ #Inf, then l2
  ylimit <- c(min(mat[,c(1+(error-1)*2, 6+(error-1)*2, 11+(error-1)*2),case]),
              max(mat[,c(1+(error-1)*2, 6+(error-1)*2, 11+(error-1)*2),case]))
  plot( N_val, mat[,1+(error-1)*2,case],  type = 'o', lty = 1, pch = 1, xlab = '#T', ylab = ylabs2[error], ylim = ylimit, main = paste('Case',case));grid() #DO
  lines(N_val, mat[,6+(error-1)*2,case],  type = 'o', lty = 2, pch = 2, xlab = '#T', ylab = ylabs2[error]) #CS
  lines(N_val, mat[,11+(error-1)*2,case], type = 'o', lty = 3, pch = 4, xlab = '#T', ylab = ylabs2[error]) #MCS

  legend('topright',legend = ADI, lty = seq(1,3,1), pch = c(1,2,4), cex = 0.8)
}




#########################################
########## Accuracy for M fixed #########
#########################################
wd <- paste(path0, "s800v5-2", sep="/")

unif <- 'N'
case <- 4

#Coord set
m1 <- seq(20, 300, 20)
m2 <- rep(50,length(m1))

adi_set <- c('DO', 'CS', 'MCS')
for(i in 1:3){
  ADI <- adi_set[i]
  an_heston_prices(cases, gridparam, m1, m2, unif, ADI, case, wd)
}

z_DO <- an_heston_errors(cases, gridparam, m1, m2, unif, 'DO', case, wd)$region
z_CS <- an_heston_errors(cases, gridparam, m1, m2, unif, 'CS', case, wd)$region
z_MCS <- an_heston_errors(cases, gridparam, m1, m2, unif, 'MCS', case, wd)$region

dim(z_DO)

imp <- matrix(NA, length(m1), 2)
imp[2:length(m1),1] <- (z_DO[1:(length(m1)-1),1] - z_DO[2:length(m1),1]) / z_DO[1:(length(m1)-1),1]
imp[2:length(m1),2] <- (z_DO[1:(length(m1)-1),2] - z_DO[2:length(m1),2]) / z_DO[1:(length(m1)-1),2]

mat <- matrix(cbind(m1, z_DO[,1], imp[,1], z_DO[,2], imp[,2], z_DO[,3]), ncol = 6)
colnames(mat) <- c('m1', 'l2', 'impr', 'inf', 'impr', 'run time'); mat


#Compute prices
for(cas in 1:4){
  #an_heston_prices(cases, gridparam, m1, m2, unif = 'N', ADI, case=cas, wd)
}


#########################################
######### All schemes and cases #########
#########################################
wd <- paste(path0, "s800v5-2", sep="/")
gridparam$N <- 101
unif <- 'N'

#Define grid points
m2_val <- seq(10,100,10)
gridparam
#Compute prices (4 hours on m2=seq(10,100,10))
schemes <- c('DO', 'CS', 'MCS')
getwd()
for(scm in 1:3){
  ADI <- schemes[scm]
  for(cas in 1:2){
    an_heston_prices(cases, gridparam, 2*m2_val, m2_val, unif = unif, ADI, case=cas, wd)
  }
}

