Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1) # only necessary for Linux without the nodejs library / headers
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

require(rstan)
install.packages("mnormt", repos = "https://cloud.r-project.org/")
require(mnormt)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


smd.weighted <- function(ps,x,z){
  x.t <- subset(x,z==1)
  x.c <- subset(x,z==0)
  om <- z/ps + (1-z)/(1-ps)
  om.t <- subset(om,z==1)
  om.c <- subset(om,z==0)
  mean.t <- sum(om.t*x.t)/sum(om.t)
  mean.c <- sum(om.c*x.c)/sum(om.c)
  var.t <- sum(om.t)/(sum(om.t)^2 - sum(om.t^2)) * sum(om.t*(x.t - mean.t)^2)
  var.c <- sum(om.c)/(sum(om.c)^2 - sum(om.c^2)) * sum(om.c*(x.c - mean.c)^2)
  return(abs(mean.t - mean.c)/sqrt((var.t + var.c)/2))
}

m <- 50
nrep <- 4
expit <- function(x){return(1/(1+exp(-x)))}
gamma <- c(1,1,1)
delta <- c(1,1,1,1)
index <- rep(1:m,each=nrep)
betaZ <- 0
RE_X1 <- rnorm(m,0,.4)
RE_X2 <- rnorm(m,0,.4)
X1 <- rnorm(m*nrep,0,.1) + RE_X1[index]
X2 <- rnorm(m*nrep,0,.1) + RE_X2[index]
beta.ests <-  var.beta.ests <- replicate(3,data.frame())
nsim <- 1000
prop.score.1 <- prop.score.2 <- replicate(nsim,data.frame())
SMD.X1 <- SMD.X2 <- replicate(3,data.frame())
seeds <- c()
data <- replicate(3,data.frame())
max.Rhat.Zre <- max.Rhat.Z <- matrix(NA,3,nsim)

for(k in 2){
  data[[k]] <- replicate(1000,data.frame())
  SMD.X1[[k]] <- SMD.X2[[k]] <- matrix(NA,nsim,3)
  object1 <- stan_model("/project/6003552/widloro/git/exposure_bernoulli_model_re.stan")
  load("/project/6003552/widloro/git/TBsim_binbin_case2Exp.RData")
  indexes <- which(max.Rhat.Zre[1,] > 1.06)
  set.seed(123456)
  for(w in indexes){
    a <- rnorm(m,0,1)
    if(k == 1) {b <- rnorm(m,0,1)}
    if(k == 2) {ab <- rmnorm(m,rep(0,2),matrix(c(1,.5,.5,1),2,2,byrow=TRUE));
    a <- ab[,1];b=ab[,2]}
    if(k == 3) {b <- a}
    Z <- rbinom(m*nrep,1,expit(cbind(1,X1,X2)%*%gamma + a[index])) 
    muY <- Z * betaZ + cbind(1,X1,X2,X1*X2)%*%delta + b[index] + rnorm(m*nrep)
    Y <-  rbinom(m*nrep,1,expit(muY)) 
    data[[k]][[w]] <- list(N = m*nrep,M=m,I=index, Z = Z, Y = Y, 
                           X=cbind(1,X1,X2),q=dim(cbind(1,X1,X2))[2])
    ############################ ajuste ######################################
    fitExp.2 <- sampling(object1, data = data[[k]][[w]],chains = 2,
                         iter = 8000)
    chain.al <- extract(fitExp.2,'alpha')
    chain.re <- extract(fitExp.2,'indRE')
    prop.score.2[[w]] <- expit(c(cbind(1,X1,X2)%*%apply(chain.al$alpha,2,mean) + apply(chain.re$indRE,2,mean)[index]))
    SMD.X1[[k]][w,3] <- smd.weighted(ps=prop.score.2[[w]],x=X1,z=data[[k]][[w]]$Z)
    SMD.X2[[k]][w,3] <- smd.weighted(ps=prop.score.2[[w]],x=X2,z=data[[k]][[w]]$Z)
    max.Rhat.Zre[k,w] <- max(stan_rhat(fitExp.2)$`data`)
    if(w %in% seq(50,nsim,len=20)){print(w);print(timestamp());save.image("/project/6003552/widloro/git/TBsim_binbin_case2Exp.RData")}
  }
  
  object2 <- stan_model("/project/6003552/widloro/git/exposure_bernoulli_model.stan")
  for(w in 1:nsim){
    fitExp.1 <- sampling(object2, data = data[[k]][[w]],  chains = 2,
                     iter = 4000)
    chain.al <- extract(fitExp.1,'alpha')
    prop.score.1[[w]] <- expit(c(cbind(1,X1,X2)%*%apply(chain.al$alpha,2,mean))) 
    SMD.X1[[k]][w,1] <- smd.weighted(ps=rep(1,m*nrep),x=X1,z=data[[k]][[w]]$Z)
    SMD.X1[[k]][w,2] <- smd.weighted(ps=prop.score.1[[w]],x=X1,z=data[[k]][[w]]$Z)
    SMD.X2[[k]][w,1] <- smd.weighted(ps=rep(1,m*nrep),x=X2,z=data[[k]][[w]]$Z)
    SMD.X2[[k]][w,2] <- smd.weighted(ps=prop.score.1[[w]],x=X2,z=data[[k]][[w]]$Z)
    max.Rhat.Z[k,w] <- max(stan_rhat(fitExp.1)$`data`)
    if(w %in% seq(50,nsim,len=20)){print(w);print(timestamp());save.image("/project/6003552/widloro/git/TBsim_binbin_case2Exp.RData")}
  }
}
