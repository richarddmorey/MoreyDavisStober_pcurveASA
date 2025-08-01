#' Estimate the power of the 2014 and 2015 p curve Test 1 ("evidential value")
#' using Fast Fourier Transforms and (optionally and slowly) Monte Carlo
#'
#' @param mu1 Effect size for first group of studies
#' @param k1 Number of studies with effect size mu1
#' @param mu2 Effect size for second group of studies
#' @param k2 Number of studies with effect size mu2
#' @param test the test for which the power is estimated: "2015" (2015 test 1), 
#'        or "z2sum" (sum of the squared z scores)
#' @param alpha alpha for test 
#' @param alpha_bound alpha for bound for p hacking/pub bias assumption
#' @param mc_iter Set number of iterations with which to check the result with Monte Carlo. Default is no checking
#' @param only_pow (logical) Return only power as a single number?
#' @param return_samples (logical) Return Monte Carlo samples (only if only_pow is FALSE)?
#' @param ... parameters to pass to fft functions
#'
#' @return
#' @export
#'
#' @examples
power_test1 <- function(mu1, k1, mu2=NULL, k2=0, df = 1, test = "2015", alpha = .05, alpha_bound = .05, mc_iter = 0, only_pow = mc_iter < 1, return_samples = FALSE, ...){
  test = match.arg(test,choices = c("2015","z2sum","2014"))
  k1 = ceiling(max(k1,0))
  k2 = ceiling(max(k2,0))
  if(k1 + k2 < 1)
    stop("One of k1 or k2 must be > 0.")
  if(is.null(mu2) & k2 != 0)
    stop("If k2 is not 0, you must specify mu2")
  mu1 = abs(mu1)
  mu2 = abs(mu2)
  if(alpha<=0 | alpha>=1)
    stop("Invalid alpha.")
  if(alpha_bound<=0 | alpha_bound>=1)
    stop("Invalid alpha_bound.")
  k = sort(c(k1,k2))
  mu = c(mu1, mu2)[order(c(k1,k2))]
  if(test == "2015"){
    left_tr = NULL
    crit = qnorm(alpha)
    dens_fun = \(x, mu) d2015_test1(x, mu, k = sum(k), alpha = alpha_bound, df = df)
  }else if(test == "z2sum"){
    left_tr = qchisq(alpha_bound,df,lower.tail = FALSE)
    crit = get_crit_dtchisq(k = sum(k), alpha = alpha, alpha_bound = alpha_bound, df = df)
    dens_fun = \(x, mu) dtchisq(x, df=df, mu = mu, tr = qchisq(1-alpha_bound,df=df))
  }else if(test == '2014'){
    left_tr = 0
    crit = qchisq(alpha,2*sum(k), lower.tail = FALSE)
    dens_fun = \(x, mu) d2014_test1(x, mu, alpha = alpha_bound, df=df)
  }else{
    stop("Invalid test: must be one of '2014', '2015', or 'z2sum'")
  }
  if(k[1]==0){
    dsum <- dsumf1(\(x)  dens_fun(x, mu[2]), k1 = k[2], left_tr = left_tr, ...)
  }else{
    dsum = dsumf2(
      f1 = \(x) dens_fun(x, mu[1]), 
      f2 = \(x) dens_fun(x, mu[2]),
      k1 = k[1],
      k2 = k[2],
      left_tr = left_tr,
      ...
    )
  }
  pow = mean(dsum$cs[which(dsum$x>crit, arr.ind = TRUE)[1] + c(-1,0)])
  if(test != "2015")
    pow = 1-pow
  if(isTRUE(only_pow)) return(pow)
  mc_pow = NULL
  mc_samples = NULL
  if(mc_iter>0){
    if(df == 1){
      z0 = replicate(mc_iter, c(rtfoldednorm(k[1],mu[1]),rtfoldednorm(k[2],mu[2]))^2)
    }else{
      z0 = replicate(mc_iter, c(rtchisq(k[1],mu[1],df=df),rtchisq(k[2],mu[2],df=df)))
    }
    if(test == "2015"){
      mc_samples = apply(z0, 2, \(z) pcurve_X2(z)[,1,2])[2,]
      mc_pow = mean(mc_samples<alpha)
    }else if(test == "z2sum"){
      mc_samples=colSums(z0)
      mc_pow = mean(mc_samples>crit)
    }else if(test == '2014'){
      mc_samples = apply(z0, 2, \(z) pcurve_X2(z)[,1,1])[2,]
      mc_pow = mean(mc_samples<alpha)
    }
  }
  list(
    test = test,
    pow = pow,
    mu = mu,
    k = k,
    mc_iter = mc_iter,
    mc_pow = mc_pow,
    mc_se = sqrt(mc_pow * (1-mc_pow) / mc_iter),
    mc_samples = if(isTRUE(return_samples)){mc_samples}else{ NULL} 
  )
}

power_test2 <- function(mu1, k1, mu2=NULL, k2=0, df=1, test = "2015", alpha = .05, tuning_power = 1/3, alpha_bound = .05, mc_iter = 0, only_pow = mc_iter < 1, return_samples = FALSE, ...){
  test = match.arg(test,choices = c("2015","z2sum","2014"))
  k1 = ceiling(max(k1,0))
  k2 = ceiling(max(k2,0))
  if(k1 + k2 < 1)
    stop("One of k1 or k2 must be > 0.")
  if(is.null(mu2) & k2 != 0)
    stop("If k2 is not 0, you must specify mu2")
  tr = qchisq(alpha_bound, df, lower.tail = FALSE)
  lambda_tuning = find_ncp_chi2_uniroot(tuning_power = tuning_power, alpha = alpha_bound, df = df)
  mu1 = abs(mu1)
  mu2 = abs(mu2)
  if(alpha<=0 | alpha>=1)
    stop("Invalid alpha.")
  if(alpha_bound<=0 | alpha_bound>=1)
    stop("Invalid alpha_bound.")
  k = sort(c(k1,k2))
  mu = c(mu1, mu2)[order(c(k1,k2))]
  if(test == "2015"){
    left_tr = NULL
    crit = qnorm(alpha)
    dens_fun = \(x, mu) d2015_test2(x, mu, k = sum(k), alpha = alpha_bound,df=df)
  }else if(test == "z2sum"){
    left_tr = tr
    crit = get_crit_dtchisq2(k = sum(k), alpha = alpha, alpha_bound = alpha_bound, mu = sqrt(lambda_tuning),df=df)
    dens_fun = \(x, mu) dtchisq(x, df=df, mu = mu, tr = qchisq(1-alpha_bound, df))
  }else if(test == '2014'){
    left_tr = 0
    crit = qchisq(alpha,2*sum(k), lower.tail = FALSE)
    dens_fun = \(x, mu) d2014_test2(x, mu, alpha = alpha_bound,df=df)
  }else{
    stop("Invalid test: must be one of '2014', '2015', or 'z2sum'")
  }
  if(k[1]==0){
    dsum <- dsumf1(\(x)  dens_fun(x, mu[2]), k1 = k[2], left_tr = left_tr, ...)
  }else{
    dsum = dsumf2(
      f1 = \(x) dens_fun(x, mu[1]), 
      f2 = \(x) dens_fun(x, mu[2]),
      k1 = k[1],
      k2 = k[2],
      left_tr = left_tr,
      ...
    )
  }
  pow = mean(dsum$cs[which(dsum$x>crit, arr.ind = TRUE)[1] + c(-1,0)])
  if(test == "2014")
    pow = 1-pow
  if(isTRUE(only_pow)) return(pow)
  mc_pow = NULL
  mc_samples = NULL
  if(mc_iter>0){
    if(df == 1){
      z0 = replicate(mc_iter, c(rtfoldednorm(k[1],mu[1]),rtfoldednorm(k[2],mu[2]))^2)
    }else{
      z0 = replicate(mc_iter, c(rtchisq(k[1],mu[1],df=df),rtchisq(k[2],mu[2],df=df)))
    }
    if(test == "2015"){
      mc_samples = apply(z0, 2, \(z) pcurve_X2(z)[,1,2])[2,]
      mc_pow = mean(mc_samples<alpha)
    }else if(test == "z2sum"){
      mc_samples=colSums(z0)
      mc_pow = mean(mc_samples<crit)
    }else if(test == '2014'){
      mc_samples = apply(z0, 2, \(z) pcurve_X2(z)[,1,1])[2,]
      mc_pow = mean(mc_samples<alpha)
    }
  }
  list(
    test = test,
    pow = pow,
    mu = mu,
    k = k,
    mc_iter = mc_iter,
    mc_pow = mc_pow,
    mc_se = sqrt(mc_pow * (1-mc_pow) / mc_iter),
    mc_samples = if(isTRUE(return_samples)){mc_samples}else{ NULL} 
  )
}

## Estimate power of test using Tippett's method on reversed p values ("familywise")
power_fw_correction = function(mu1, k1, mu2=NULL, k2=0, fw_alpha = 0.05, tuning_power = 1/3, alpha_bound = 0.05, df = 1){
  fw_crit = qtchisq(1-(1-fw_alpha)^(1/(k1+k2)), tuning_power = tuning_power, alpha_bound = alpha_bound, df = df)
  crit = qchisq(alpha_bound,df,lower.tail = FALSE)
  c0 = pchisq(crit,df,mu1^2)
  pow1 = (pchisq(fw_crit,df,mu1^2) -  c0)/(1-c0)
  if(k2>=1){
    c0 = pchisq(crit,df,mu2^2)
    pow2 = (pchisq(fw_crit,df,mu2^2) -  c0)/(1-c0)
    pow = 1-((1-pow1)^k1 * (1-pow2)^k2)
  }else{
    pow = 1-(1-pow1)^k1
  }
  return(pow)
}


## Estimate power of test using hypercube rejection region around origin
power_min_test_correction = function(mu1, k1, mu2=NULL, k2=0, fw_alpha = 0.05, tuning_power = 1/3, alpha_bound = 0.05, df = 1){
  
  crit = qchisq(alpha_bound,df,lower.tail = FALSE)
  c0 = pchisq(crit,df,mu1^2)
  fw_crit = qtchisq(fw_alpha^(1/(k1+k2)), tuning_power = tuning_power,alpha_bound = alpha_bound, df = df)
  
  pow1 = (pchisq(fw_crit,df,mu1^2) -  c0)/(1-c0)
  if(k2>=1){
    c0 = pchisq(crit,df,mu2^2)
    pow2 = (pchisq(fw_crit,df,mu2^2) -  c0)/(1-c0)
    pow = pow1^k1 * pow2^k2
  }else{
    pow = pow1^k1
  }
  return(pow)
}
