code = nimbleCode({
  for(i in 1:t){
    mu[i] ~ dnorm(0, tau = 1/1000)
  }
  delta[1] ~ dunif(-1,1)
  delta[2] ~ dunif(-sqrt(1-delta[1]^2),sqrt(1-delta[1]^2))
  ssa[1:t,1:t] <- diag(t)
  Sigma[1:t,1:t] ~ dinvwish(ssa[1:t,1:t],t+1)
  Sz[1:t,1:t] <- diag(t)-(2/pi)*(delta[1:t]%*%t(delta[1:t]))
  Szsq[1:t,1:t] <- inverse(t(chol(Sz[1:t,1:t])))
  Ssq[1:t,1:t] <- t(chol(Sigma[1:t,1:t]))
  Delta[1:t] <- Ssq[1:t,1:t]%*%Szsq[1:t,1:t]%*%delta[1:t]
  tau[1:t,1:t] <- Ssq[1:t,1:t]%*%Szsq[1:t,1:t]%*%(diag(t)-delta[1:t]%*%t(delta[1:t]))%*%t(Szsq[1:t,1:t])%*%t(Ssq[1:t,1:t])
  nu ~ T(dexp(theta), 2, Inf)
  theta ~ dunif(0.05, 0.49)
  for(i in 1:n){
    h[i] ~ T(dnorm(0, 1), 0, Inf)
    u[i] ~ dgamma(nu/2,nu/2)
    muy[1:t,i] <- mu[1:t] + (1/sqrt(u[i]))*Delta[1:t]*(h[i]-sqrt(2/pi))
    precy[1:t,1:t,i] <- u[i]*inverse(tau[1:t,1:t])
    y[i,] ~ dmnorm(muy[1:t,i], prec = precy[1:t,1:t,i])
  }
})