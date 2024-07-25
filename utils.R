library(mvtnorm)
library(sn)
library(pracma)
library(fdrtool)
library(rstan)
library(gdata)
library(bayesplot)
library(nimble)
library(extraDistr)
library(matrixStats)
library(MCMCglmm)
library(FMsmsnReg)
library(ggplot2)
library(moments)
library(cubature)
library(matrixcalc)

#####################################
#### Observed Mahalanobis distance
#####################################
# y = the response data of dimension (nxd)
# mu = the mean vector
# Sigma = the covariance parameter of dimension (dxd)
# delta = the vector of skewness
# d = dimension of the random vector
MD.SMMCSN = function(y, mu, Sigma, delta, d){
  if(ncol(y) != d){
    stop("The number of columns of y must be d")
  }
  if(d < 1){
    stop("d must be greater than 0")
  }
  m = 0
  b = sqrt(2/pi)
  muz = b*delta
  Sigmaz = diag(d)-muz%*%t(muz)
  smz = t(chol(Sigmaz))
  sm = t(chol(Sigma))
  xi = mu-sm%*%solve(smz)%*%muz
  om = sm%*%solve(Sigmaz)%*%t(sm)
  om1 = solve(om)
  for(i in 1:n){
    m[i] = t(y[i,]-xi)%*%om1%*%(y[i,]-xi)
  }
  return(m)
}

########################################
#### Healy plot for Mahalanobis distance
########################################
# m = Observed Mahalanobis distance
# d = dimension of the random vector
# dist = distribution (MCSN, MCST, MCSS or MCSCN)
# nu = shape parameter (NULL if dist = MCSN)
healy.SMMCSN = function(m, d, dist, nu = NULL){
  if(all(c("MCST","MCSS","MCSCN", "MCSN") != dist)){
    stop("The distribution is not defined.")
  }
  if(d < 1){
    stop("d must be greater than 0.")
  }
  if(!is.null(nu) & dist == "MCSN"){
    stop("nu must be NULL if the distribution is the MCSN.")
  }
  fn = ecdf(m)
  ms = sort(m)
  aux1 = fn(ms)
  if(dist == "MCSN"){
    aux2 = pchisq(ms,d)
  }
  if(dist == "MCST"){
    aux2 = pf(ms/d,d,nu)
  }
  if(dist == "MCSS"){
    aux2 = pchisq(ms,d) - 2^(nu)*gamma(d/2+nu)*pchisq(ms,d+2*nu)/(ms^(nu)*gamma(d/2))
  }
  if(dist == "MCSCN"){
    nu1 = nu[1]
    nu2 = nu[2]
    aux2 = nu1*pchisq(nu2*ms,d) + (1-nu1)*pchisq(ms,d)
  }
  ptg = ggplot() + geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x = aux1, y = aux2)) +
    xlab("Empirical distribution") + 
    ylab("Theoretical distribution") +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "grey"),
          text=element_text(size=15,family="serif"))
  print(ptg)
}

#####################################
#### QQ plot for Mahalanobis distance
#####################################
# res = Observed Mahalanobis distance
# d = dimension of the random vector
# dist = distribution (MCSN, MCST, MCSS or MCSCN)
# nu = shape parameter (NULL if dist = MCSN)
# confi = coefficient of confidence of the envelopes
# repl = number of replicates used to construct the envelopes
# splotlim = superior limit of the ggplot
qq.SMMCSN = function(res, d, dist, nu = NULL, confi = 0.95, 
                     repl = 100, splotlim = NULL){
  if(all(c("MCST","MCSS","MCSCN", "MCSN") != dist)){
    stop("The distribution is not defined")
  }
  if(d < 1){
    stop("d must be greater than 0")
  }
  if(!is.null(nu) & dist == "MCSN"){
    stop("nu must be NULL if the distribution is the MCSN.")
  }
  n = length(res)
  rres = matrix(0, nrow = n, ncol = repl)
  sres = sort(res)
  k = 0
  if(dist == "MCSN"){
    while(k<repl){
      k = k+1
      aux20 = rchisq(n,d)
      rres[,k] = aux20
    }
  }
  if(dist == "MCST"){
    while(k<repl){
      k = k+1
      aux20 = d*rf(n, d, nu)
      rres[,k] = aux20
    }
  }
  if(dist == "MCSS"){
    Fd = function(x, d, nu){
      aux = pchisq(x, d) - (2^nu*gamma(d/2+nu)/(x^nu*gamma(d/2)))*pchisq(x, d+2*nu)
      #if(x<0.000001){aux = 0}
      return(aux)
    }
    while(k<repl){
      k = k+1
      aux20 = 0
      for(i in 1:n){
        u = runif(1)
        aux20[i] = optim(1, function(x) (Fd(x, d, nu)-u)^2)$par
      }
      rres[,k] = aux20
    }
  }
  if(dist == "MCSCN"){
    nu1 = nu[1]
    nu2 = nu[2]
    Fd = function(x, d, nu){
      aux = nu1*pchisq(x*nu2, d) + (1-nu1)*pchisq(x, d)
      return(aux)
    }
    while(k<repl){
      k = k+1
      aux20 = 0
      for(i in 1:n){
        u = runif(1)
        aux20[i] = optim(1, function(x) (Fd(x, d, nu)-u)^2)$par
      }
      rres[,k] = aux20
    }
  }
  
  srres = matrix(0, nrow = n, ncol = repl)
  
  for(k in 1:repl) {srres[,k]=sort(rres[,k])}
  descq = matrix(0, nrow = n, ncol = 3)
  for(k in 1:n){
    descq[k,1] = quantile(srres[k,],probs = (1-confi)/2)
    descq[k,2] = median(srres[k,])
    descq[k,3] = quantile(srres[k,],probs = 1-(1-confi)/2)
  }
  if(is.null(splotlim)){splotlim = max(sres)+10}
  if(dist == "MCSN"){
    aux50 = qchisq(ppoints(n),d)
    ptg = ggplot() + geom_point(aes(x = aux50, y = sres)) +
      geom_line(aes(x = aux50, y = descq[,1])) +
      geom_line(aes(x = aux50, y = descq[,3])) +
      geom_line(aes(x = aux50, y = descq[,2]), linetype = "dashed") +
      labs(y = "Observed Mahalanobis distance", 
           x = "Chi-square quantiles")+
      ylim(0,splotlim) +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15)
      )
  }
  if(dist == "MCST"){
    aux50 = d*qf(ppoints(n), d, nu)
    ptg = ggplot() + geom_point(aes(x = aux50, y = sres)) +
      geom_line(aes(x = aux50, y = descq[,1])) +
      geom_line(aes(x = aux50, y = descq[,3])) +
      geom_line(aes(x = aux50, y = descq[,2]), linetype = "dashed") +
      labs(y = "Observed Mahalanobis distance", 
           x = "Theoretical Mahalanobis distance quantiles")+
      ylim(0,splotlim) +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15)
      )
  }
  if(dist == "MCSS"){
    auxpp = ppoints(n)
    aux50 = 0
    for(i in 1:length(auxpp)){
      aux50[i] = optim(1, function(x) (Fd(x, d, nu)-auxpp[i])^2)$par
    }
    ptg = ggplot() + geom_point(aes(x = aux50, y = sres)) +
      geom_line(aes(x = aux50, y = descq[,1])) +
      geom_line(aes(x = aux50, y = descq[,3])) +
      geom_line(aes(x = aux50, y = descq[,2]), linetype = "dashed") +
      labs(y = "Observed Mahalanobis distance", 
           x = "Theoretical Mahalanobis distance quantiles")+
      ylim(0,splotlim) +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15)
      )
  }
  if(dist == "MCSCN"){
    auxpp = ppoints(n)
    aux50 = 0
    for(i in 1:length(auxpp)){
      aux50[i] = optim(1, function(x) (Fd(x, d, nu)-auxpp[i])^2)$par
    }
    ptg = ggplot() + geom_point(aes(x = aux50, y = sres)) +
      geom_line(aes(x = aux50, y = descq[,1])) +
      geom_line(aes(x = aux50, y = descq[,3])) +
      geom_line(aes(x = aux50, y = descq[,2]), linetype = "dashed") +
      labs(y = "Observed Mahalanobis distance", 
           x = "Theoretical Mahalanobis distance quantiles")+
      ylim(0,splotlim) +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15)
      )
  }
  print(ptg)
}
# y = the response data of dimension (nxd)
# mu = the mean vector
# Sigma = the covariance parameter of dimension (dxd)
# delta = the vector of skewness
dMCSN = function(y, mu, Sigma, delta){
  d = length(y)
  b = sqrt(2/pi)
  lamb = delta/sqrt(1-t(delta)%*%delta)
  sm = t(chol(Sigma))
  muz = b*delta
  Sigmaz = diag(d)-muz%*%t(muz)
  smz = t(chol(Sigmaz))
  xi = mu-sm%*%solve(smz)%*%muz
  omb = sm%*%solve(smz)
  om = omb%*%t(omb)
  dens = 2*dmvnorm(as.vector(y-xi),rep(0,d),om)*pnorm(t(lamb)%*%solve(t(chol(om)))%*%(y-xi))
  return(dens)
}
# n = number of simulated values
# mu = the mean vector
# Sigma = the covariance parameter of dimension (dxd)
# delta = the vector of skewness
rMCSN = function(n, mu, Sigma, delta){
  d = ncol(Sigma)
  b = sqrt(2/pi)
  lamb = delta/sqrt(1-sum(delta^2))
  sm = t(chol(Sigma))
  muz = b*delta
  Sigmaz = diag(d)-muz%*%t(muz)
  smz = t(chol(Sigmaz))
  xi = mu-sm%*%solve(smz)%*%muz
  omb = sm%*%solve(smz)
  simul = matrix(0,n,d)
  for(i in 1:n){
    simul[i,] = xi+omb%*%(delta*abs(rnorm(1))+t(chol(diag(d)-delta%*%t(delta)))%*%t(rmvnorm(1, mean = rep(0,d))))
  }
  return(simul)
}
# y = the response data of dimension (nxd)
# mu = the mean vector
# Sigma = the covariance parameter of dimension (dxd)
# delta = the vector of skewness
# nu = the shape parameter
dMCST = function(y,mu,Sigma,delta,nu){
  integrand = function(x,mu,Sigma,delta,nu){
    as.numeric(dMCSN(y,mu,(1/x)*Sigma,delta))*dgamma(x,nu/2,nu/2)
  }
  aux = cubintegrate(integrand,0,Inf,mu=mu,Sigma=Sigma,delta=delta,nu=nu,
                     method = "hcubature")$integral
  return(aux)
}
# y = the response data of dimension (nxd)
# mu = the mean vector
# Sigma = the covariance parameter of dimension (dxd)
# delta = the vector of skewness
# nu = the shape parameter
dMCSS = function(y,mu,Sigma,delta,nu){
  integrand = function(x,mu,Sigma,delta,nu){
    as.numeric(dMCSN(y,mu,(1/x)*Sigma,delta))*dbeta(x,nu,1)
  }
  aux = cubintegrate(integrand,0,1,mu=mu,Sigma=Sigma,delta=delta,nu=nu,
                     method = "hcubature")$integral
  return(aux)
}
# y = the response data of dimension (nxd)
# mu = the mean vector
# Sigma = the covariance parameter of dimension (dxd)
# delta = the vector of skewness
# nu = the vector of shape parameters
dMCSCN = function(y,mu,Sigma,delta,nu){
  aux = nu[1]*dMCSN(y,mu,(1/nu[2])*Sigma,delta)+(1-nu[1])*dMCSN(y,mu,Sigma,delta)
  return(aux)
}