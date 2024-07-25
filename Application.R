#rm(list = ls(all.names = TRUE))
source("utils.R") # Auxiliary functions

# Data set
library(MVT)
data("WindSpeed")
dad = WindSpeed
y = cbind(dad[,1], dad[,3])
ggplot() +
  geom_point(aes(x = y[,1], y = y[,2])) +
  theme_bw() +
  xlab("vs") + ylab("kw") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))
ggplot() +
  geom_histogram(aes(x = y[,1], y=..density..), 
                 colour="black", fill="white") +
  xlab(expression(y[1])) + ylab("Density") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))
ggplot() +
  geom_histogram(aes(x = y[,2], y=..density..), 
                 colour="black", fill="white") +
  xlab(expression(y[1])) + ylab("Density") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))
n = nrow(y)
t = ncol(y)

#########################################
# Multivariate centered skew-normal model
#########################################

source("MCSN.R") # code used in nimbleMCMC()
constants = list(n = n, t = t, pi = pi)
data <- list(y = y)
inits <- function() list(mu = colMeans(y), Sigma = cov(y),
                         delta = rep(0, t),
                         h = abs(rnorm(n)))
samples = nimbleMCMC(code = code, data = data, 
                     inits = inits, constants = constants,
                     monitors = c("mu","Sigma", "delta"),
                     thin = 50, niter = 120000, 
                     nburnin = 80000, nchains = 1)
library(coda)
theta = colMeans(samples) # posterior means
tsd = colSds(as.matrix(samples)) # posterior standard errors

# Posterior statistics, standard errors and HPD
mes = theta[7:8]; round(mes,2)
round(tsd[7:8],2)
round(HPDinterval(as.mcmc(samples[,7:8])),2)
ses = matrix(theta[1:4], ncol = 2, byrow = F); round(ses,2)
round(tsd[1:4],2)
round(HPDinterval(as.mcmc(samples[,1:4])),2)
des = theta[5:6]; round(des,2)
round(tsd[5:6],2)
round(HPDinterval(as.mcmc(samples[,5:6])),2)

# Mahalanobis distance
res = MD.SMMCSN(y, mes, ses, des, t)
healy.SMMCSN(res,t,"MCSN") # Healy plot
qq.SMMCSN(res, t, "MCSN") # QQ plot

# Bivariate density
xx = seq(-43,52, 2)
yy = seq(-33,57,2)
length(xx)
z = matrix(0, length(xx), length(yy))
for(i in 1:length(xx)){
  print(i)
  for(j in 1:length(yy)){
    z[i,j] = dMCSN(c(xx[i],yy[j]), mu = mes, ses, delta = des)
  }
}

x1 = expand.grid(xx,yy)
ggplot() +
  geom_point(aes(x = y[,1], y = y[,2]), alpha = 0.5) +
  geom_contour(aes(x = x1[,1], y = x1[,2], z = vec(z)),
               colour = "black") +
  xlab("vs") + ylab("kw") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))

# Calculate the marginal parameters
b = sqrt(2/pi)
sm = t(chol(ses))
muz = b*des
Sigmaz = diag(t)-muz%*%t(muz)
om = sm%*%solve(Sigmaz)%*%sm
omsqrt = t(chol(om))
d1s = (1/omsqrt[1,1])*c(1,0)%*%omsqrt%*%des; d1s
d2s = (1/omsqrt[2,2])*c(0,1)%*%omsqrt%*%des; d2s

# Mariginal densities
z = rep(0,length(xx))
for(j in 1:length(xx)){
  print(j)
  z[j] = dMCSN(xx[j], mes[1], ses[1,1], d1s)
}
df = data.frame(y = y[,1])
ggplot() +
  geom_histogram(aes(x = y[,1], y=..density..), colour="black", fill="white") +
  geom_line(aes(x = xx,y = z)) +
  xlab("vs") + ylab("Density") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))


z = rep(0,length(yy))
for(j in 1:length(yy)){
  print(j)
  z[j] = dMCSN(yy[j], mes[2], ses[2,2], d2s)
}
ggplot() +
  geom_histogram(aes(x = y[,2], y=..density..), colour="black", fill="white") +
  geom_line(aes(x = yy,y = z)) +
  xlab("kw") + ylab("Density") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))

####################################
# Multivariate centered skew-t model
####################################
source("MCST.R") # code used in nimbleMCMC()
constants = list(n = n, t = t, pi = pi)
data <- list(y = y)
inits <- function() list(mu = colMeans(y), Sigma = cov(y),
                         delta = rep(0, t),
                         h = abs(rnorm(n)), nu = runif(1,3,10),
                         theta = runif(1,0.05,0.49),
                         u = rgamma(n, 4/2, 4/2))
samples <- nimbleMCMC(code = code, data = data, inits = inits,
                      constants = constants,
                      monitors = c("mu","Sigma", "delta", "nu"),
                      thin = 50, niter = 120000, 
                      nburnin = 80000,
                      nchains = 1)
theta = colMeans(samples) # posterior means
theta2 = posterior.mode(samples) # posterior modes
tsd = colSds(as.matrix(samples)) # posterior standard errors

# Posterior statistics, standard errors and HPD
mes = theta[7:8]; round(mes,2)
round(tsd[7:8],2)
round(HPDinterval(as.mcmc(samples[,7:8])),2)
ses = matrix(theta[1:4], ncol = 2, byrow = F); round(ses,2)
round(tsd[1:4],2)
round(HPDinterval(as.mcmc(samples[,1:4])),2)
des = theta2[5:6]; round(des,2)
round(tsd[5:6],2)
round(HPDinterval(as.mcmc(samples[,5:6])),2)
nues = theta2[9]; round(nues,2)
round(tsd[9],2)
round(HPDinterval(as.mcmc(samples[,9])),2)

# Mahalanobis distance
res = MD.SMMCSN(y, mes, ses, des, t)
healy.SMMCSN(res,t,"MCST", nu = nues) # Healy plot
qq.SMMCSN(res, t, "MCST", nu = nues) # QQ plot

# Bivariate density
xx = seq(-43,52, 2)
yy = seq(-33,57,2)
z = matrix(0, length(xx), length(yy))
for(i in 1:length(xx)){
  print(i)
  for(j in 1:length(yy)){
    z[i,j] = dMCST(c(xx[i],yy[j]), mu = mes, ses, delta = des, 
                   nu = nues)
  }
}

x1 = expand.grid(xx,yy)
ggplot() +
  geom_point(aes(x = y[,1], y = y[,2]), alpha = 0.5) +
  geom_contour(aes(x = x1[,1], y = x1[,2], z = vec(z)),
               colour = "black") +
  xlab("vs") + ylab("kw") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))

# Calculate the marginal parameters
b = sqrt(2/pi)
sm = t(chol(ses))
muz = b*des
Sigmaz = diag(t)-muz%*%t(muz)
om = sm%*%solve(Sigmaz)%*%sm
omsqrt = t(chol(om))
d1s = (1/omsqrt[1,1])*c(1,0)%*%omsqrt%*%des; d1s
d2s = (1/omsqrt[2,2])*c(0,1)%*%omsqrt%*%des; d2s

# Mariginal densities
z = rep(0,length(xx))
for(j in 1:length(xx)){
  print(j)
  z[j] = dMCST(xx[j], mes[1], ses[1,1], d1s, nu = nues)
}
df = data.frame(y = y[,1])
ggplot() +
  geom_histogram(aes(x = y[,1], y=..density..), colour="black", 
                 fill="white") +
  geom_line(aes(x = xx,y = z)) +
  xlab("vs") + ylab("Density") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))



z = rep(0,length(yy))
for(j in 1:length(yy)){
  print(j)
  z[j] = dMCST(yy[j], mes[2], ses[2,2], d2s, nu = nues)
}
ggplot() +
  geom_histogram(aes(x = y[,2], y=..density..), 
                 colour="black", fill="white") +
  geom_line(aes(x = yy,y = z)) +
  xlab("kw") + ylab("Density") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))

##########################################
### Multivariate centered skew slash model
##########################################

source("MCSS.R") # code used in nimbleMCMC()
constants = list(n = n, t = t, pi = pi)
data <- list(y = y)
inits <- function() list(mu = colMeans(y), Sigma = cov(y),
                         delta = rep(0, t),
                         h = abs(rnorm(n)), nu = runif(1,2,10),
                         theta = runif(1,0.05,0.99),
                         u = rbeta(n, 4, 1))
samples <- nimbleMCMC(code = code, data = data, inits = inits,
                      constants = constants,
                      monitors = c("mu","Sigma", "delta", "nu"),
                      thin = 50, niter = 120000, 
                      nburnin = 80000,
                      nchains = 1)
theta = colMeans(samples) # posterior means
theta2 = posterior.mode(samples) # posterior modes
tsd = colSds(as.matrix(samples)) # posterior standard errors

# Posterior statistics, standard errors and HPD
mes = theta[7:8]; round(mes,2)
round(tsd[7:8],2)
round(HPDinterval(as.mcmc(samples[,7:8])),2)
ses = matrix(theta[1:4], ncol = 2, byrow = F); round(ses,2)
round(tsd[1:4],2)
round(HPDinterval(as.mcmc(samples[,1:4])),2)
des = theta2[5:6]; round(des,2)
round(tsd[5:6],2)
round(HPDinterval(as.mcmc(samples[,5:6])),2)
nues = theta2[9]; round(nues,2)
round(tsd[9],2)
round(HPDinterval(as.mcmc(samples[,9])),2)

# Mahalanobis distance
res = MD.SMMCSN(y, mes, ses, des, t)
healy.SMMCSN(res, t, "MCSS", nu = nues) # Healy plot
qq.SMMCSN(res, t, "MCSS", nu = nues) # QQ plot

# Bivariate density
xx = seq(-43,52, 2)
yy = seq(-33,57,2)
z = matrix(0, length(xx), length(yy))
for(i in 1:length(xx)){
  print(i)
  for(j in 1:length(yy)){
    z[i,j] = dMCSS(c(xx[i],yy[j]), mu = mes, ses, delta = des, 
                   nu = nues)
  }
}

x1 = expand.grid(xx,yy)
ggplot() +
  geom_point(aes(x = y[,1], y = y[,2]), alpha = 0.5) +
  geom_contour(aes(x = x1[,1], y = x1[,2], z = vec(z)),
               colour = "black") +
  xlab("vs") + ylab("kw") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))

# Calculate the marginal parameters
b = sqrt(2/pi)
sm = t(chol(ses))
muz = b*des
Sigmaz = diag(t)-muz%*%t(muz)
om = sm%*%solve(Sigmaz)%*%sm
omsqrt = t(chol(om))
d1s = (1/omsqrt[1,1])*c(1,0)%*%omsqrt%*%des
d2s = (1/omsqrt[2,2])*c(0,1)%*%omsqrt%*%des

# Mariginal densities
z = rep(0,length(xx))
for(j in 1:length(xx)){
  print(j)
  z[j] = dMCSS(xx[j], mes[1], ses[1,1], d1s, nu = nues)
}
df = data.frame(y = y[,1])
ggplot() +
  geom_histogram(aes(x = y[,1], y=..density..), colour="black", fill="white") +
  geom_line(aes(x = xx,y = z)) +
  xlab("vs") + ylab("Density") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))


z = rep(0,length(yy))
for(j in 1:length(yy)){
  print(j)
  z[j] = dMCSS(yy[j], mes[2], ses[2,2], d2s, nu = nues)
}
ggplot() +
  geom_histogram(aes(x = y[,2], y=..density..), colour="black", fill="white") +
  geom_line(aes(x = yy,y = z)) +
  xlab('kw') + ylab("Density") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))

######################################################
# Multivariate centered skew contaminated normal model
######################################################

source("MCSCN.R") # code used in nimbleMCMC()
constants = list(n = n, t = t, pi = pi)
data <- list(y = y)
inits <- function() list(mu = colMeans(y), Sigma = cov(y),
                         delta = rep(0, t),
                         h = abs(rnorm(n)), nu1 = runif(1),
                         nu2 = runif(1),
                         v = rbinom(n,1,0.5))
samples <- nimbleMCMC(code = code, data = data, inits = inits,
                      constants = constants,
                      monitors = c("mu","Sigma", "delta", "nu1",
                                   "nu2"),
                      thin = 50, niter = 120000, 
                      nburnin = 80000,
                      nchains = 1)
theta = colMeans(samples) # posterior means
theta2 = posterior.mode(samples) # posterior modes

tsd = colSds(as.matrix(samples)) # posterior standard errors

# Posterior statistics, standard errors and HPD
mes = theta[7:8]; round(mes,2)
round(tsd[7:8],2)
round(HPDinterval(as.mcmc(samples[,7:8])),2)
ses = matrix(theta[1:4], ncol = 2, byrow = F); round(ses,2)
round(tsd[1:4],2)
round(HPDinterval(as.mcmc(samples[,1:4])),2)
des = theta2[5:6]; round(des,2)
round(tsd[5:6],2)
round(HPDinterval(as.mcmc(samples[,5:6])),2)
nues1 = theta[9]; round(nues1,2)
round(tsd[9],2)
round(HPDinterval(as.mcmc(samples[,9])),2)
nues2 = theta[10]; round(nues2,2)
round(tsd[10],2)
round(HPDinterval(as.mcmc(samples[,10])),2)

# Mahalanobis distance
res = MD.SMMCSN(y, mes, ses, des, t)
healy.SMMCSN(res, t, "MCSCN", nu = c(nues1,nues2)) # Healy plot
qq.SMMCSN(res, t, "MCSCN", nu = c(nues1,nues2)) # QQ plot


# Bivariate density
xx = seq(-43,52, 2)
yy = seq(-33,57,2)
z = matrix(0, length(xx), length(yy))
for(i in 1:length(xx)){
  print(i)
  for(j in 1:length(yy)){
    z[i,j] = dMCSCN(c(xx[i],yy[j]), mu = mes, ses, delta = des, 
                   nu = c(nues1,nues2))
  }
}

x1 = expand.grid(xx,yy)
ggplot() +
  geom_point(aes(x = y[,1], y = y[,2]), alpha = 0.5) +
  geom_contour(aes(x = x1[,1], y = x1[,2], z = vec(z)),
               colour = "black", bins = 20) +
  xlab("vs") + ylab("kw") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))

# Calculate the marginal parameters
b = sqrt(2/pi)
sm = t(chol(ses))
muz = b*des
Sigmaz = diag(t)-muz%*%t(muz)
om = sm%*%solve(Sigmaz)%*%sm
omsqrt = t(chol(om))
d1s = (1/omsqrt[1,1])*c(1,0)%*%omsqrt%*%des; d1s
d2s = (1/omsqrt[2,2])*c(0,1)%*%omsqrt%*%des; d2s

# Mariginal densities
z = rep(0,length(xx))
for(j in 1:length(xx)){
  print(j)
  z[j] = dMCSCN(xx[j], mes[1], ses[1,1], d1s, nu = c(nues1,nues2))
}
df = data.frame(y = y[,1])
ggplot() +
  geom_histogram(aes(x = y[,1], y=..density..), colour="black", fill="white") +
  geom_line(aes(x = xx,y = z)) +
  xlab("vs") + ylab("Density") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))



z = rep(0,length(yy))
for(j in 1:length(yy)){
  print(j)
  z[j] = dMCSCN(yy[j], mes[2], ses[2,2], d2s, nu = c(nues1,nues2))
}
ggplot() +
  geom_histogram(aes(x = y[,2], y=..density..), 
                 colour="black", fill="white") +
  geom_line(aes(x = yy,y = z)) +
  xlab("kw") + ylab("Density") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15,family="serif"))
