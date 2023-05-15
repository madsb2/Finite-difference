Heston.Fourier <- function(spot,timetoexp,strike,r,divyield,V,theta,kappa,epsilon,rho,greek=1){
  
  X<-log(spot/strike)+(r-divyield)*timetoexp
  kappahat<-kappa-0.5*rho*epsilon
  xiDummy<-kappahat^2+0.25*epsilon^2
  
  integrand<-function(k){
    xi<-sqrt(k^2*epsilon^2*(1-rho^2)+2i*k*epsilon*rho*kappahat+xiDummy)
    Psi.P<--(1i*k*rho*epsilon+kappahat)+xi
    Psi.M<-(1i*k*rho*epsilon+kappahat)+xi
    alpha<--kappa*theta*(Psi.P*timetoexp + 2*log((Psi.M + Psi.P*exp(-xi*timetoexp))/(2*xi)))/epsilon^2
    beta<--(1-exp(-xi*timetoexp))/(Psi.M + Psi.P*exp(-xi*timetoexp))
    numerator<-exp((-1i*k+0.5)*X+alpha+(k^2+0.25)*beta*V)
    
    if(greek==1) dummy<-Re(numerator/(k^2+0.25))
    if(greek==2) dummy<-Re((0.5-1i*k)*numerator/(spot*(k^2+0.25)))
    if(greek==3) dummy<--Re(numerator/spot^2)
    if(greek==4) dummy<-Re(numerator*beta)
    
    #integrand<-dummy
    return(dummy)
  } 
  
  dummy<-integrate(integrand,lower=-100,upper=100,stop.on.error = FALSE)$value
  
  if (greek==1) dummy<-exp(-divyield*timetoexp)*spot-strike*exp(-r*timetoexp)*dummy/(2*pi)
  
  if(greek==2) dummy<-exp(-divyield*timetoexp)-strike*exp(-r*timetoexp)*dummy/(2*pi)
  
  if(greek==3) dummy<--strike*exp(-r*timetoexp)*dummy/(2*pi)
  
  if(greek==4) dummy<--strike*exp(-r*timetoexp)*dummy/(2*pi)
  
  #Heston.Fourier<-dummy
  return(dummy)
}

Heston_call <- function(param, grid){
  s_vector <- grid$s_vector; J <- length(s_vector)
  v_vector <- grid$v_vector; M <- length(v_vector)
  t_vector <- grid$t_vector; N <- length(t_vector)
  
  strike <- param$K
  r <- param$r
  theta <- param$eta
  kappa <- param$kappa
  rho <- param$rho
  sigma <- param$sigma
  
  prices <- matrix(NA,J,M)
  
  for(j in 2:J){
    for(m in 1:M){
      prices[m,j] <- Heston.Fourier(spot = s_vector[j],timetoexp = t_vector[N],strike = strike,r = r,
                divyield = 0,V = v_vector[m],theta = theta, kappa = kappa, epsilon = sigma, rho = rho)
    }
  }
  prices[,1] <- 0
  return(prices)
}

#param <- list(r = 0.025, sigma = 0.3, K = 10,theta = 0.8,
#              rho=-0.9,kappa=1.5,eta=0.4)
#S <- 15
#V <- 0.5

Andreasen.Fourier <- function(spot,timetoexp,strike,Z,lambda,beta,epsilon){
  
  X<-log(spot/strike)
  
  integrand<-function(k){
    neweps<-lambda*epsilon
    xi<-sqrt(k^2*neweps^2+beta^2+0.25*neweps^2)
    Psi.P<--beta+xi
    Psi.M<-beta+xi
    A<--beta*(Psi.P*timetoexp + 2*log((Psi.M + Psi.P*exp(-xi*timetoexp))/(2*xi)))/(epsilon^2)
    B<-(1-exp(-xi*timetoexp))/(Psi.M + Psi.P*exp(-xi*timetoexp))
    integrand<-Re(exp( (-1i*k+0.5)*X+A-(k^2+0.25)*B*lambda^2*Z)/(k^2+0.25))
  } 
  
  dummy<-integrate(integrand,lower=-Inf,upper=Inf)$value
  
  #Andreasen.Fourier<-spot-strike*dummy/(2*pi)
  return(spot-strike*dummy/(2*pi))
  
}

