# Poisson Deviance Loss

Poisson.Deviance <- function(pred, obs){
  2*(sum(pred)-sum(obs)+sum(log((obs/pred)^(obs))))/length(pred)
}


Gamma.Deviance <-function(pred, obs, wei){
  2*sum(wei*((obs/pred)-1-log(obs/pred)))/length(pred)
}

gen.trun(par = c(0),
         family = PO,
         type = "left") ## generate the truncated poisson distribution functions


minus_ZTP_logLik <- function(pred, obs){
  -sum(dPOtr(obs, mu=pred,log = TRUE)) / length(pred)
}


z1<-function(lambda){
  qnorm((1-dpois(0,lambda)-dpois(1,lambda))/(1-dpois(0,lambda)))
}

z2<-function(phi,mu,sev,num){
  alpha<-1/phi
  beta<-1/(phi*mu)
  qnorm(pgamma(q = sev, shape = num*alpha, rate = num*beta))
}

z2_ln<-function(phi,ln_mu,ln_sev){
  qnorm(pnorm(q = ln_sev, mean=ln_mu, sd=sqrt(phi)))
}

z2_lga<-function(phi,lga_mu,ln_sev){
  alpha<-1/phi
  beta<-1/(phi*lga_mu)
  qnorm(pgamma(q = ln_sev, shape = alpha, rate = beta))
}

rho_f<-function(ind,z1,z2,rho){
  logp_1 = pnorm(q = -(z1-rho*z2)/sqrt(1-rho^2), log.p = TRUE)
  logp_2 = pnorm(q = (z1-rho*z2)/sqrt(1-rho^2), log.p = TRUE)
  -sum(ind*logp_1+(1-ind)*logp_2)
}

cdf_T = function(t, lambda){
  if(t >1){
    cdf_value = (1-dpois(0,lambda)-lambda*exp(-lambda*t))/(1-dpois(0,lambda))
  }
  else{
    cdf_value = (1-dpois(0,lambda*t)-dpois(1,lambda*t))/(1-dpois(0,lambda))
  }
  return(cdf_value)
}

cdf_T(0.9,0.1)

inverse_cdf_T <- function(q, lambda){
  uniroot(function(x) q-cdf_T(x,lambda),lower=0,upper=10^9)$root
}
inverse_cdf_T(0.04,0.1)
  
pred_S <- function(lambda, mu, phi, rho_hat, TT,dis){
  S <- rep(0, TT) ## aggregated claims amount
  for(i in 1:TT){
    T_1 = rexp(1, lambda)
    if(T_1 > 1){
      S[i] = 0
    }
    else{
      z = mvrnorm(n=1, c(0,0), matrix(c(1, rho_hat, rho_hat, 1), 2, 2))# Generate  $(z_{i,1},z_{i,2})$ from the bivariate normal distribution with correlation $\hat{\rho}$
      t = inverse_cdf_T(pnorm(z[1]), lambda = lambda)# Recover $T_{n+1}$ by inverting its distribution function
      N = 1
      while(t <= 1){
        N = N + 1
        t = t + rexp(1, lambda)
      }
      # if (t<=1) {N=2}
      # if (t>1) {N=1}
      if (dis=="gamma"){
         alpha = 1/phi # shape parameter of gamma distribution
         beta = 1/(phi*mu) # rate parameter of gamma distribution
         S[i] = N*qgamma(pnorm(z[2]), shape = alpha*N, rate = beta*N)
        }
      if (dis=="lnorm"){
          S[i] = N*exp(qnorm(pnorm(z[2]), mean=mu, sd=sqrt(phi) ) )
      }
      if (dis=="lgamma"){
        alpha = 1/phi # shape parameter of gamma distribution
        beta = 1/(phi*mu) # rate parameter of gamma distribution
        S[i] = N*exp(qgamma(pnorm(z[2]), shape=alpha, rate=beta ) )
      }
      }
    }
  return(S)
}

D1 <- function(u, v, rho){
  #rho = (exp(2*g_rho)-1)/(exp(2*g_rho)+1)
  d_1 = pnorm(q = (qnorm(v)-rho*qnorm(u))/(sqrt(1-rho^2)),
              mean = 0,
              sd = 1)
  return(d_1)
}

ifm_rho_kramer <- function(ClaimNb, lambda, Sev, mu, 
                           phi, ln){
  
  alpha = 1/phi # shape parameter of gamma distribution
  beta = alpha/mu # rate parameter of gamma distribution
  
  # calculate u, v, w from Kramer(2013)
  if (ln==F) {
    u = pgamma(q = Sev, shape = ClaimNb*alpha, rate = ClaimNb*beta)
  }
  if (ln==T){
    u = pnorm(q = Sev, mean = mu, sd = sqrt(phi) )
  }
  v = pPOtr(q = ClaimNb, mu=lambda)
  w = pPOtr(q = ClaimNb - 1, mu=lambda)

  #minus joint log-likelihood
  object_fun <- function(rho){
    -sum(log(D1(u, v, rho)-D1(u, w, rho)))
  }
  
  ifm_joint <- optim(par=0,
                     fn = object_fun,
                     lower=-1,
                     upper=1,
                     method = 'Brent')
  # rho_hat = (exp(2*ifm_joint$par)-1)/(exp(2*ifm_joint$par)+1)
  rho_hat = ifm_joint$par
  return(rho_hat)
}

pred_kramer_S <- function(lambda, mu, phi, rho, TT,ln){
  S <- rep(0,TT)
  alpha = 1/phi # shape parameter of gamma distribution
  beta = 1/(phi*mu) # rate parameter of gamma distribution
  for(i in 1:TT){
    if (ln==F){
      sev <- rgamma(1,shape = alpha, rate = beta)
      u = pgamma(q = sev, shape = alpha, rate = beta)
    }
    if (ln==T){
      sev <- exp(rnorm(1, mean = mu, sd = sqrt(phi)) )
      u = pnorm(q = log(sev), mean = mu, sd = sqrt(phi) )
    }
    Pr_N<-rep(NA,4)
    Pr_N[0+1] <- ppois(0, lambda)
    for (k in 1:2){
      v = pPOtr(q = k, mu=lambda)
      w = pPOtr(q = k-1, mu=lambda)
      Pr_N[k+1] <- (D1(u,v,rho)-D1(u,w,rho))*(1-Pr_N[0+1])
    }
    Pr_N[4]<-1-sum(Pr_N[1:3])
    R=which(as.vector(rmultinom(1,1,prob=Pr_N))==1)-1
    S[i]=R*sev
    if (is.na(S[i])==1) stop("NA generated") 
    }
 return(S)
}

pred_condition_S <- function(lambda, mu, phi, coef_N, TT, dis){
  S <- rep(0,TT)
  for(i in 1:TT){
    N <- rpois(1, lambda)
    if (N == 0){
      S[i] = 0
    } 
    else{
      if (dis=="gamma"){
        mu1 = mu * exp(N*coef_N)
        alpha = 1/phi # shape parameter of gamma distribution
        beta = 1/(phi*mu1) # rate parameter of gamma distribution
        S[i] <- N*rgamma(1, shape = N*alpha, rate = N*beta)
      }
      if (dis=="lnorm"){
        mu1 = mu + N*coef_N
        S[i] <- N*exp(rnorm(1, mean = mu1, sd = sqrt(phi)))
      }
      if (dis=="lgamma"){
        mu1 = mu + N*coef_N
        alpha = 1/phi
        beta = 1/(phi*mu1)
        S[i] <- N*exp(rgamma(1, shape = alpha, rate=beta))
      }
    }
  }
  return(S)
}

agg_per<-function(simulated, actual){
  actual<-data.frame(ClaimAmount=actual)
  actual$group_ind<-rep(1:100, 100)[1:nrow(actual)]
  group_sum<-aggregate(actual$ClaimAmount, by=list(actual$group_ind), sum)
  simulated<-data.frame(simulated)
  simulated$group_ind<-rep(1:100, 100)[1:nrow(actual)]
  group_simulated<-matrix(NA,nrow=nrow(simulated),ncol=TT)
  for (i in 1:100){
    group_simulated[i,]<-apply(simulated[simulated$group_ind==i,1:500],2,sum)
  }
  group_percent<-rep(NA,100)
  for (i in 1:100){
    ecdf_temp<-ecdf(group_simulated[i,])
    group_percent[i]<-ecdf_temp(group_sum$x[i])
  }
  return(group_percent)
}

