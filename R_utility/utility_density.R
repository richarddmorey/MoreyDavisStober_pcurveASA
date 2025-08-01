###############
# For Chi/ChiSq distribution
##############

dtchisq <- function(x2, df=1, mu=0, tr=qchisq(.05, df, lower.tail = FALSE)){
  area = pchisq(tr, df, mu^2, lower.tail = FALSE)
  dchisq(x2, df, mu^2) / area * (x2 >= tr)
}

rtchisq <- function(n, df=1, mu=0, tr=qchisq(.05, df, lower.tail = FALSE))
  qchisq(runif(n,pchisq(tr,df,mu^2),1),df,mu^2)

qtchisq = function(alpha,tuning_power = 1/3, alpha_bound = .05, df=1){
  crit = qchisq(alpha_bound,df,lower.tail = FALSE)
  lambda0 = find_ncp_chi2_uniroot(tuning_power = tuning_power,alpha = alpha_bound,df=df)
  # Recompute tuning_power for numerical reasons
  tp0 = pchisq(crit,df,lambda0,lower.tail = FALSE)
  qchisq(1 - tp0*(1-alpha),df,lambda0) 
}

dtfoldednorm <- function(z, mu=0, sigma=1, tr=qnorm(.975)){
  area = pnorm(-abs(tr),mu,sigma) + pnorm(abs(tr), mu, sigma, lower.tail = FALSE)
  (dnorm(z, mu, sigma) + dnorm(z, -mu, sigma))/area * (z>=tr)    
}

rtfoldednorm <- function(n, mu=0, sigma=1, tr=qnorm(.975)){
  lower_area = pnorm(-abs(tr),mu,sigma)
  upper_area = pnorm(abs(tr), mu, sigma, lower.tail = FALSE)
  u = runif(n,0,lower_area + upper_area)
  u[u>lower_area] = u[u>lower_area] + (1 - lower_area - upper_area)
  abs(qnorm(u, mu, sigma))
}

# Used for sum(abs(Z)) test
get_crit_dtfoldednorm = function(k, alpha=.05, alpha_bound = .05){
  d = dsumf1(\(x) dtfoldednorm(x, 0), k, a = qnorm(alpha_bound, lower.tail = FALSE))
  mean(d$x[which(d$cs>(1-alpha), arr.ind = TRUE)[1] + c(-1,0)])
}

# Used for sum(abs(Z)) test
# This version allows computation of test 2 critical values
get_crit_dtfoldednorm2 = function(k, alpha=.05, alpha_bound = .05, mu = 0, lower.tail = mu > 0){
  d = dsumf1(\(x) dtfoldednorm(x, mu), k, a = qnorm(alpha_bound, lower.tail = FALSE))
  if(lower.tail) alpha = 1-alpha
  mean(d$x[which(d$cs>(1-alpha), arr.ind = TRUE)[1] + c(-1,0)])
}

# Used for sum(X2) test
get_crit_dtchisq = memoise::memoise(function(k, alpha=.05, alpha_bound = .05, df = 1){
  d = dsumf1(\(x) dtchisq(x, df, 0), k, a = qchisq(1-alpha_bound,df))
  mean(d$x[which(d$cs>(1-alpha), arr.ind = TRUE)[1] + c(-1,0)])
})

# Used for sum(X2) test
# This version allows computation of test 2 critical values
get_crit_dtchisq2 = memoise::memoise(function(k, alpha=.05, alpha_bound = .05, mu = 0, lower.tail = mu > 0, df=1){
  d = dsumf1(\(x) dtchisq(x, df, mu), k, a = qchisq(1-alpha_bound,df))
  if(lower.tail) alpha = 1-alpha
  mean(d$x[which(d$cs>(1-alpha), arr.ind = TRUE)[1] + c(-1,0)])
})

# Density function for test statistic in 2014 Test 1 (EV)
d2014_test1 = Vectorize(function(z2, mu, k, alpha = .05, df  = 1){
  crit_z2 = qchisq(alpha,df,lower.tail = FALSE)
  if(z2<0) return(0)
  pow = pchisq(crit_z2,df,mu^2,lower.tail = FALSE)
  p = exp(-z2/2)
  z2_0 = qchisq(p*alpha,df,lower.tail = FALSE)
  d=dchisq(z2_0,df,mu^2) / pow / dchisq(z2_0,df) * alpha * exp(-z2/2) / 2
  if(is.na(d)) d = 0
  return(d)
},'z2')

# Density function for test statistic in 2015 Test 1 (EV*)
d2015_test1 = Vectorize(function(z, mu, k, alpha = .05, df = 1){
  crit_z2 = qchisq(alpha,df,lower.tail = FALSE)
  pow = pchisq(crit_z2,df,mu^2,lower.tail = FALSE)
  p = pnorm(z*sqrt(k))
  z2 = qchisq(p*alpha,df,lower.tail = FALSE)
  d=dchisq(z2,df,mu^2) / pow / dchisq(z2,df) * alpha * dnorm(z*sqrt(k)) * sqrt(k)
  if(is.na(d)) d = 0
  return(d)
},'z')

# Density function for test statistic in 2014 Test 2 (LEV)
d2014_test2 = Vectorize(function(z2, mu, k, tuning_power = 1/3, alpha = .05, df = 1){
  tr = qchisq(alpha, df, lower.tail = FALSE)
  lambda = mu^2
  lambda0 = find_ncp_chi2_uniroot(tuning_power = tuning_power, alpha = alpha, df = df)
  c0 = pchisq(tr,df,lambda0) # Recompute the power for numerical reasons
  if(z2 <= 0) return(0)
  x = qchisq((1-c0)*(exp(-z2/2) + c0/(1-c0)),df,lambda0) 
  f = .5 * dchisq(x,df,lambda) * (1-c0) * exp(-z2/2) / dchisq(x,df,lambda0) / pchisq(tr,df,lambda,lower.tail = FALSE)
  return(ifelse(is.nan(f),0,f))
},"z2")

# Density function for test statistic in 2015 Test 2 (EV*)
d2015_test2 = Vectorize(function(z2, mu, k, tuning_power = 1/3, alpha = .05, df = 1){
  tr = qchisq(alpha, df, lower.tail = FALSE)
  lambda = mu^2
  lambda0 = find_ncp_chi2_uniroot(tuning_power = tuning_power, alpha = alpha, df = df)
  c0 = pchisq(tr,df,lambda0) # Recompute the power for numerical reasons
  x = qchisq((1-c0)*(pnorm(z2*sqrt(k)) + c0/(1-c0)),df,lambda0) 
  f = sqrt(k) * dchisq(x,df,lambda) * (1-c0) * dnorm(z2*sqrt(k)) / dchisq(x,df,lambda0) / pchisq(tr,df,lambda,lower.tail = FALSE)
  return(ifelse(is.nan(f),0,f))
},"z2")

###############
# For F distribution
##############

rtf <- function(n, df1=1, df2, ncp=0, tr=qf(.05, df1, df2, lower.tail = FALSE))
  qf(runif(n,pf(tr,df1,df2,ncp),1),df1,df2,ncp)

dtf <- function(x, df1=1, df2, ncp=0, tr=qf(.05, df1, df2, lower.tail = FALSE)){
  area = pf(tr, df1, df2, ncp, lower.tail = FALSE)
  df(x, df1, df2, ncp) / area * (x >= tr)
}

qtf = function(alpha,tuning_power = 1/3, alpha_bound = .05, df1=1, df2){
  crit = qf(alpha_bound,df1,df2,lower.tail = FALSE)
  lambda0 = find_ncp_F_uniroot(tuning_power = tuning_power,alpha = alpha_bound,df1=df1,df2=df2)
  # Recompute tuning_power for numerical reasons
  tp0 = pf(crit,df1,df2,lambda0,lower.tail = FALSE)
  qf(1 - tp0*(1-alpha),df1,df2,lambda0) 
}

dtfoldedt <- function(x, df, ncp=0, tr=qt(.975, df)){
  area = pt(-abs(tr),df,ncp) + pt(abs(tr), df, ncp, lower.tail = FALSE)
  (dt(x, df, ncp) + dt(-x, df, ncp))/area * (x>=tr)    
}

rtfoldedt <- function(n, df, ncp=0, tr=qt(.975, df)){
  lower_area = pt(-abs(tr),df,ncp)
  upper_area = pt(abs(tr), df, ncp, lower.tail = FALSE)
  u = runif(n,0,lower_area + upper_area)
  u[u>lower_area] = u[u>lower_area] + (1 - lower_area - upper_area)
  abs(qt(u, df, ncp))
}

# Used for sum(F) test
get_crit_dtf = function(k, alpha=.05, alpha_bound = .05, df1 = 1, df2){
  d = dsumf1(\(x) dtf(x, df1, df2, 0), k, a = qf(1-alpha_bound,df1, df2), sc=30, nx = 2^13)
  mean(d$x[which(d$cs>(1-alpha), arr.ind = TRUE)[1] + c(-1,0)])
}

# Used for sum(F) test
# This version allows computation of test 2 critical values
get_crit_dtf2 = function(k, alpha=.05, alpha_bound = .05, ncp = 0, lower.tail = ncp > 0, df1=1, df2){
  d = dsumf1(\(x) dtf(x, df1, df2, ncp), k, a = qf(1-alpha_bound,df1,df2), sc=30, nx = 2^13)
  if(lower.tail) alpha = 1-alpha
  mean(d$x[which(d$cs>(1-alpha), arr.ind = TRUE)[1] + c(-1,0)])
}


## Three functions for transformation sum(F/(1+F))

dtransf1 = function(x, df1, df2, ncp, alpha_bound = .05, log=FALSE){
  f_crit = qf(alpha_bound,df1,df2,lower.tail = FALSE)
  x_crit = f_crit / (1 + f_crit)
  area = pf(f_crit,df1,df2,ncp,lower.tail = FALSE)
  i = x>x_crit & x<1
  fval = x[i]/(1-x[i])
  d = rep(NA,length(x))
  d[i] = df(fval, df1, df2, ncp) / (1-x[i])^2  / area
  d[!i] = 0
  if(isTRUE(log)){
    return(log(d))
  }else{
    return(d)
  }
}

# Used for sum(F) test
get_crit_dtransf1 = function(k, alpha=.05, alpha_bound = .05, df1 = 1, df2){
  tr = qf(1-alpha_bound,df1, df2)
  tr = tr / (1 + tr)
  d = dsumf1(\(x) dtransf1(x, df1, df2, ncp=0), k, a = tr, b = 1, moments_interval = c(tr,1))
  mean(d$x[which(d$cs>(1-alpha), arr.ind = TRUE)[1] + c(-1,0)])
}

# Used for sum(F) test
# This version allows computation of test 2 critical values
get_crit_dtransf1_LEV = function(k, alpha=.05, alpha_bound = .05, ncp = 0, lower.tail = ncp > 0, df1=1, df2){
  tr = qf(1-alpha_bound,df1, df2)
  tr = tr / (1 + tr)
  d = dsumf1(\(x) dtransf1(x, df1, df2, ncp), k, a = tr, b = 1, moments_interval = c(tr,1))
  if(lower.tail) alpha = 1-alpha
  mean(d$x[which(d$cs>(1-alpha), arr.ind = TRUE)[1] + c(-1,0)])
}

# Density function for test statistic in 2014 Test 1 (EV)
d2014_test1_F = Vectorize(function(x2, ncp, k, alpha = .05, df1  = 1, df2){
  crit_F = qf(alpha,df1,df2,lower.tail = FALSE)
  if(x2<0) return(0)
  pow = pf(crit_F,df1,df2,ncp,lower.tail = FALSE, log=TRUE)
  p = -x2/2
  x2_0 = qf(p + log(alpha),df1,df2,lower.tail = FALSE, log.p = TRUE)
  logd=df(x2_0,df1,df2,ncp,log=TRUE) - pow - df(x2_0,df1,df2,log=TRUE) + log(alpha) -x2/2 - log(2)
  ifelse(is.na(logd),0,exp(logd))
},'x2')

# Density function for test statistic in 2014 Test 2 (LEV)
d2014_test2_F = Vectorize(function(x2, ncp, k, tuning_power = 1/3, alpha = .05, df1 = 1, df2){
  crit_F = qf(alpha, df1, df2, lower.tail = FALSE)
  ncp0 = find_ncp_F_uniroot(tuning_power = tuning_power, alpha = alpha, df1 = df1, df2 = df2)
  c0 = pf(crit_F,df1,df2,ncp0) # Recompute the power for numerical reasons
  if(x2 <= 0) return(0)
  x = qf((1-c0)*(exp(-x2/2) + c0/(1-c0)),df1,df2,ncp0) 
  f = .5 * df(x,df1,df2,ncp) * (1-c0) * exp(-x2/2) / df(x,df1,df2,ncp0) / pf(crit_F,df1,df2,ncp,lower.tail = FALSE)
  return(ifelse(is.nan(f),0,f))
},"x2")

# Density function for test statistic in 2015 Test 1 (EV*)
d2015_test1_F = Vectorize(function(z, ncp, k, alpha = .05, df1 = 1, df2){
  crit_F = qf(alpha, df1, df2, lower.tail = FALSE)
  pow = pf(crit_F,df1,df2,ncp,lower.tail = FALSE, log.p = TRUE)
  p = pnorm(z*sqrt(k))
  z2 = qf(p*alpha,df1,df2,lower.tail = FALSE)
  logd = df(z2,df1,df2,ncp,log = TRUE) - pow - df(z2,df1,df2, log = TRUE) + 
    log(alpha) + dnorm(z*sqrt(k), log = TRUE) + .5*log(k)
  #d=df(z2,df1,df2,ncp) / pow / df(z2,df1,df2) * alpha * dnorm(z*sqrt(k)) * sqrt(k)
  ifelse(is.na(logd),0,exp(logd))
},'z')

# Density function for test statistic in 2015 Test 2 (LEV*)
d2015_test2_F = Vectorize(function(z, ncp, k, tuning_power = 1/3, alpha = .05, df1 = 1, df2){
  crit_F = qf(alpha,df1,df2,lower.tail = FALSE)
  ncp0 = find_ncp_F_uniroot(tuning_power = tuning_power, alpha = alpha, df1 = df1, df2 = df2)
  c0 = pf(crit_F,df1,df2,ncp0) # Recompute the power for numerical reasons
  x = qf((1-c0)*(pnorm(z*sqrt(k)) + c0/(1-c0)),df1,df2,ncp0) 
  f = sqrt(k) * df(x,df1,df2,ncp) * (1-c0) * dnorm(z*sqrt(k)) / df(x,df1,df2,ncp0) / pf(crit_F,df1,df2,ncp,lower.tail = FALSE)
  return(ifelse(is.nan(f),0,f))
},"z")
