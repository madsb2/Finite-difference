################################################
################## Plot 1 ######################
################################################
library('rgl')

FDplot1 <- function(x,y,z,xlabel = 'x', ylabel = 'y',title = 'Finite difference'){
  open3d()
  persp3d(x,y,z,col='blue',xlab = '',ylab = '', zlab = '',axes = FALSE)
  axes3d(edge = c('x+','y','z++'),nticks = 6,cex = 0.7)
  title3d(main=title,ylab = ylabel, zlab = '', line = 4)
  ##############
  mtext3d(text = xlabel, edge = "x+-", line = 10, at = x[length(x)]) 
  
  #Adding gridlines
  n <- min(15,length(x), length(y))
  xgrid <- x[seq(1, length(x), length.out = n)]
  ygrid <- y[seq(1, length(y), length.out = n)]
  zgrid <- matrix(NA,n,n)
  for(xg in 1:n){
    i <- which(x == xgrid[xg],arr.ind = TRUE)
    for(yg in 1:n){
      j <- which(y == ygrid[yg],arr.ind = TRUE)
      zgrid[xg,yg] <- z[i,j]
    }
  }
  persp3d(xgrid,ygrid,zgrid,front='lines',back='lines',add=TRUE)
  
  um <- matrix(c(-.7,-.7,0,0,
                 .3,-.3,.9,0,
                 -.6,.6,.3,0,
                 0,0,0,1),4,4,byrow = TRUE)
  
  view3d(userMatrix = um)
  #rgl.snapshot("Heston prices.png")
}

################################################
################## Plot 2 ######################
################################################

FDplot2 <- function(x,y,z1,z2,xlabel = 'x', ylabel = 'y',title1 = 'Finite difference', title2 = 'Analytical'){
  open3d()
  mfrow3d(1,2)
  persp3d(x,y,z1,col="deepskyblue1",xlab = '',ylab = '', zlab = '',axes = FALSE)
  axes3d(edge = c('x+','y','z++'),nticks = 6,cex = 0.7)
  title3d(main=title1,xlab = '',ylab = ylabel, zlab = '', line = 4)
  #mtext3d(xlabel, "x+-", line = 10, at = x[round(length(x)/2,0)]*2)
  text3d(x = x[round(length(x)/2,0)], y = max(y)+y[round(length(y)/2.8,0)], z = -1, texts = xlabel)
  #Adding gridlines
  n <- min(15,length(x), length(y))
  xgrid <- x[seq(1, length(x), length.out = n)]
  ygrid <- y[seq(1, length(y), length.out = n)]
  zgrid <- matrix(NA,n,n)
  for(xg in 1:n){
    i <- which(x == xgrid[xg],arr.ind = TRUE)
    for(yg in 1:n){
      j <- which(y == ygrid[yg],arr.ind = TRUE)
      zgrid[xg,yg] <- z1[i,j]
    }
  }
  persp3d(xgrid,ygrid,zgrid,front='lines',back='lines',add=TRUE)
  
  um <- matrix(c(-.7,-.7,0,0,
                 .3,-.3,.9,0,
                 -.6,.6,.3,0,
                 0,0,0,1),4,4,byrow = TRUE)
  
  view3d(userMatrix = um)
  
  #Analytical
  persp3d(x,y,z2,col="deepskyblue1",xlab = '',ylab = ylabel, zlab = '',axes = FALSE)
  axes3d(edge = c('x+','y','z++'),nticks = 6,cex = 0.7)
  title3d(main=title2,xlab = '',ylab = ylabel, zlab = '', line = 4)
  #mtext3d(xlabel, "x+-", line = 10, at = x[round(length(x)/2,0)]*2)
  text3d(x = x[round(length(x)/2,0)], y = max(y)+y[round(length(y)/2.8,0)], z = -1, texts = xlabel)
  
  for(xg in 1:n){
    i <- which(x == xgrid[xg],arr.ind = TRUE)
    for(yg in 1:n){
      j <- which(y == ygrid[yg],arr.ind = TRUE)
      zgrid[xg,yg] <- z2[i,j]
    }
  }
  persp3d(xgrid,ygrid,zgrid,front='lines',back='lines',add=TRUE)
  view3d(userMatrix = um)
  rgl.snapshot("Heston vs BS prices.png")
}

################################################
############# Heston error plot ################
################################################

heston_err_plot <- function(x,y,z,xlabel = 's', ylabel = 'v',title = 'Errors', z_upper = NA){
  if(is.na(z_upper)==TRUE){z_upper <- max(z)}
  open3d()
  persp3d(x,y,z,col='deepskyblue1',xlab = '',ylab = '', zlab = '',axes = FALSE, zlim = c(0,z_upper))
  axes3d(edge = c('x+','y+','z+'),nticks = 6,cex = 0.7)
  title3d(main=title,ylab = ylabel, zlab = '', line = 4)
  
  mtext3d(text = xlabel, edge = "x+-", line = 10, at = x[round(length(x)/8,0)])
  #Adding gridlines
  n <- min(15,length(x), length(y))
  xgrid <- x[seq(1, length(x), length.out = n)]
  ygrid <- y[seq(1, length(y), length.out = n)]
  zgrid <- matrix(NA,n,n)
  for(xg in 1:n){
    i <- which(x == xgrid[xg],arr.ind = TRUE)
    for(yg in 1:n){
      j <- which(y == ygrid[yg],arr.ind = TRUE)
      zgrid[xg,yg] <- z[i,j]
    }
  }
  persp3d(xgrid,ygrid,zgrid,front='lines',back='lines',add=TRUE, zlim = c(0,z_upper))
  
  um <- matrix(c(-.7,.7,-.1,0,
                 -.25,-.25,.9,0,
                 .6,.6,.4,0,
                 0,0,0,1),4,4,byrow = TRUE)
  
  view3d(userMatrix = um) #par3d()$userMatrix
}

#############################################
heston_err_plot2 <- function(x,y,z1, z2,xlabel = '#s', ylabel = '#v',title1 = 'Errors', title2 = 'Errors'){
  
  open3d()
  rgl.viewpoint
  mfrow3d(1,2)
  par3d('windowRect')
  ?par3d
  
  persp3d(x,y,z1,col='deepskyblue1',xlab = '',ylab = '', zlab = '',axes = FALSE)
  axes3d(edge = c('x+','y+','z+'),nticks = 6,cex = 0.7)
  title3d(main=title1,ylab = ylabel, zlab = '', line = 4)
  
  mtext3d(text = xlabel, edge = "x+-", line = 10, at = x[round(length(x)/8,0)])
  #Adding gridlines
  n <- min(15,length(x), length(y))
  xgrid <- x[seq(1, length(x), length.out = n)]
  ygrid <- y[seq(1, length(y), length.out = n)]
  zgrid <- matrix(NA,n,n)
  for(xg in 1:n){
    i <- which(x == xgrid[xg],arr.ind = TRUE)
    for(yg in 1:n){
      j <- which(y == ygrid[yg],arr.ind = TRUE)
      zgrid[xg,yg] <- z1[i,j]
    }
  }
  persp3d(xgrid,ygrid,zgrid,front='lines',back='lines',add=TRUE)
  
  um <- matrix(c(-.7,.7,-.1,0,
                 -.25,-.25,.9,0,
                 .6,.6,.4,0,
                 0,0,0,1),4,4,byrow = TRUE)
  
  view3d(userMatrix = um) #par3d()$userMatrix
  
  next3d()
  persp3d(x,y,z2,col='deepskyblue1',xlab = '',ylab = '', zlab = '',axes = FALSE)
  axes3d(edge = c('x+','y+','z+'),nticks = 6,cex = 0.7)
  title3d(main=title2,ylab = ylabel, zlab = '', line = 4)
  
  mtext3d(text = xlabel, edge = "x+-", line = 10, at = x[round(length(x)/8,0)])
  #Adding gridlines

  for(xg in 1:n){
    i <- which(x == xgrid[xg],arr.ind = TRUE)
    for(yg in 1:n){
      j <- which(y == ygrid[yg],arr.ind = TRUE)
      zgrid[xg,yg] <- z2[i,j]
    }
  }
  persp3d(xgrid,ygrid,zgrid,front='lines',back='lines',add=TRUE)
  
  um <- matrix(c(-.7,.7,-.1,0,
                 -.25,-.25,.9,0,
                 .6,.6,.4,0,
                 0,0,0,1),4,4,byrow = TRUE)
  
  view3d(userMatrix = um) #par3d()$userMatrix
}