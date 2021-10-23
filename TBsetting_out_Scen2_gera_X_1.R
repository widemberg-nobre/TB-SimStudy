Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1) # only necessary for Linux without the nodejs library / headers
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

require(rstan)
install.packages("mnormt", repos = "https://cloud.r-project.org/")
require(mnormt)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

beta.ests <-  var.beta.ests  <- OddsRatio.ests <- replicate(3,data.frame())

object4 <- stan_model("/project/6003552/widloro/git/outcome_bernoulli_model.stan")

for(k in 1:3){
  if(k==1){("/project/6003552/widloro/git/TBsim_binbin_case1Exp_Scen2_de_X.RData")}
  if(k==2){("/project/6003552/widloro/git/TBsim_binbin_case2Exp_Scen2_de_X.RData")}
  if(k==3){("/project/6003552/widloro/git/TBsim_binbin_case3Exp_Scen2_de_X.RData")}
object4 <- stan_model("/project/6003552/widloro/git/outcome_bernoulli_model.stan")
  beta.ests[[k]] <- OddsRatio.ests[[k]] <- var.beta.ests[[k]] <-  matrix(NA,nsim,4)
  for(w in 1:nsim){
    data.3 <- list(N = m*nrep,M=m,I=index, 
                   Z = data[[k]][[w]]$Z, Y = data[[k]][[w]]$Y, 
                   X=cbind(1,prop.score.1[[w]]),
                   q=dim(cbind(1,prop.score.1[[w]]))[2])
    data.4 <- list(N = m*nrep,M=m,I=index, 
                   Z = data[[k]][[w]]$Z, Y = data[[k]][[w]]$Y, 
                   X=cbind(1,prop.score.2[[w]]),
                   q=dim(cbind(1,prop.score.2[[w]]))[2])
    fitOut1 <- sampling(object4, data = data.3,  chains = 2,
                        iter = 6000)
    fitOut2 <- sampling(object4, data = data.4,  chains = 2,
                        iter = 6000)
    chain.ate1 <- extract(fitOut1,'ate')
    chain.ate2 <- extract(fitOut2,'ate')
    
    beta.ests[[k]][w,1] <- mean(chain.ate1$ate)
    beta.ests[[k]][w,2] <- mean(chain.ate2$ate)
    var.beta.ests[[k]][w,1] <- var(chain.ate1$ate)
    var.beta.ests[[k]][w,2] <- var(chain.ate2$ate)
    chain.beta1 <- extract(fitOut1,'beta')
    chain.beta2 <- extract(fitOut2,'beta')
    OddsRatio.ests[[k]][w,1] <- mean(exp(chain.beta1$beta))
    OddsRatio.ests[[k]][w,2] <- mean(exp(chain.beta2$beta))
    if(w %in% seq(50,nsim,len=20)){print(w);print(timestamp());save.image("/project/6003552/widloro/git/TBsim_binbin_Scen2_de_X_MD1_e_MD2.RData")}
  }
}
save.image("/project/6003552/widloro/git/TBsim_binbin_Scen2_de_X_MD1_e_MD2.RData")
