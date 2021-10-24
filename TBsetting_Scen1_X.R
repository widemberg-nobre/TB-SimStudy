install.packages("rstan", repos = "https://cloud.r-project.org/")
require(rstan)
require(mnormt)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

object2 <- stan_model("/project/6003552/widloro/git/exposure_bernoulli_model.stan")

for(k in 1:3){
  if(k==1){load("/project/6003552/widloro/git/TBsim_binbin_case1Exp_Scen1_de_X.RData")}
  if(k==2){load("/project/6003552/widloro/git/TBsim_binbin_case2Exp_Scen1_de_X.RData")}
  if(k==3){load("/project/6003552/widloro/git/TBsim_binbin_case3Exp_Scen1_de_X.RData")}
  for(w in 1:nsim){
    fitExp.1 <- sampling(object2, data = data[[k]][[w]],  chains = 2,iter = 4000)
    chain.al <- extract(fitExp.1,'alpha')
    prop.score.1[[w]] <- expit(c((data[[k]][[w]]$X)%*%apply(chain.al$alpha,2,mean))) 
    SMD.X1[[k]][w,1] <- smd.weighted(ps=rep(1,m*nrep),x=(data[[k]][[w]]$X)[,2],z=data[[k]][[w]]$Z)
    SMD.X1[[k]][w,2] <- smd.weighted(ps=prop.score.1[[w]],x=(data[[k]][[w]]$X)[,2],z=data[[k]][[w]]$Z)
    SMD.X2[[k]][w,1] <- smd.weighted(ps=rep(1,m*nrep),x=(data[[k]][[w]]$X)[,3],z=data[[k]][[w]]$Z)
    SMD.X2[[k]][w,2] <- smd.weighted(ps=prop.score.1[[w]],x=(data[[k]][[w]]$X)[,3],z=data[[k]][[w]]$Z)
    max.Rhat.Z[k,w] <- max(stan_rhat(fitExp.1)$`data`)
    if(w %in% seq(50,nsim,len=20)){print(w);print(timestamp());save.image("/project/6003552/widloro/git/TBsim_binbin_Exp_Scen1_de_X.RData")}
  }
}
