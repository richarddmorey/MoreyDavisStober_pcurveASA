#' Estimate the power of the 2014 and 2015 p curve Test 1 ("evidential value")
#' (with F statistics)
#' using Fast Fourier Transforms and (optionally and slowly) Monte Carlo
#'
#' In Marden (1982), the summed test statistic is defined as the ratio of the sums
#' of squares, NOT dividing the numerator and denominator by their degrees of freedom.
#' This will yield the same test if the degrees of freedom are the same across all 
#' K studies, but NOT if the degrees of freedom vary. This function will need modification
#' if you want to use degrees of freedom that differ across studies.
#'
#' @param ncp1 Effect size for first group of studies
#' @param k1 Number of studies with effect size ncp1
#' @param ncp2 Effect size for second group of studies
#' @param k2 Number of studies with effect size ncp2
#' @param test the test for which the power is estimated: "2015" (2015 test 1), 
#'        "2014" (2014 test 1), or "Fsum" 
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
power_test1_F <- function(ncp1, k1, ncp2=NULL, k2=0, df1 = 1, df2, test = "2015", alpha = .05, alpha_bound = .05, mc_iter = 0, only_pow = mc_iter < 1, return_samples = FALSE, ...){
  test = match.arg(test,choices = c("2015","Fsum","2014"))
  k1 = ceiling(max(k1,0))
  k2 = ceiling(max(k2,0))
  if(k1 + k2 < 1)
    stop("One of k1 or k2 must be > 0.")
  if(is.null(ncp2) & k2 != 0)
    stop("If k2 is not 0, you must specify ncp2")
  ncp1 = abs(ncp1)
  ncp2 = abs(ncp2)
  if(alpha<=0 | alpha>=1)
    stop("Invalid alpha.")
  if(alpha_bound<=0 | alpha_bound>=1)
    stop("Invalid alpha_bound.")
  k = sort(c(k1,k2))
  ncp = c(ncp1, ncp2)[order(c(k1,k2))]
  if(test == "2015"){
    left_tr = NULL
    crit = qnorm(alpha)
    dens_fun = \(x, ncp) d2015_test1_F(x, ncp = ncp, k = sum(k), alpha = alpha_bound, df1 = df1, df2 = df2)
  }else if(test == "Fsum"){
    left_tr = qf(alpha_bound,df1,df2,lower.tail = FALSE)
    crit = get_crit_dtf(k = sum(k), alpha = alpha, alpha_bound = alpha_bound, df1 = df1, df2 = df2)
    dens_fun = \(x, ncp) dtf(x, df1=df1,df2=df2, ncp = ncp, tr = qf(1-alpha_bound,df1=df1,df2=df2))
  }else if(test == '2014'){
    left_tr = 0
    crit = qchisq(alpha,2*sum(k), lower.tail = FALSE)
    dens_fun = \(x, ncp) d2014_test1_F(x, ncp=ncp, alpha = alpha_bound, df1 = df1, df2 = df2)
  }else{
    stop("Invalid test: must be one of '2014', '2015', or 'Fsum'")
  }
  if(k[1]==0){
    dsum <- dsumf1(\(x)  dens_fun(x, ncp[2]), k1 = k[2], left_tr = left_tr, nx=2^13, sc=30, ...)
  }else{
    dsum = dsumf2(
      f1 = \(x) dens_fun(x, ncp[1]), 
      f2 = \(x) dens_fun(x, ncp[2]),
      k1 = k[1],
      k2 = k[2],
      left_tr = left_tr,
      nx=2^13, sc=30,
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
    if(df1 == 1){
      z0 = replicate(mc_iter, c(rtfoldedt(k[1],ncp[1],df=df2),rtfoldedt(k[2],ncp[2],df=df2))^2)
    }else{
      z0 = replicate(mc_iter, c(rtf(k[1],ncp=ncp[1],df1=df1,df2=df2),rtf(k[2],ncp=ncp[2],df1=df1,df2=df2)))
    }
    if(test == "2015"){
      mc_samples = apply(z0, 2, \(z) pcurve_F(z)[,1,2])[2,]
      mc_pow = mean(mc_samples<alpha)
    }else if(test == "Fsum"){
      mc_samples=colSums(z0)
      mc_pow = mean(mc_samples>crit)
    }else if(test == '2014'){
      mc_samples = apply(z0, 2, \(z) pcurve_F(z)[,1,1])[2,]
      mc_pow = mean(mc_samples<alpha)
    }
  }
  list(
    test = test,
    pow = pow,
    ncp = ncp,
    k = k,
    mc_iter = mc_iter,
    mc_pow = mc_pow,
    mc_se = sqrt(mc_pow * (1-mc_pow) / mc_iter),
    mc_samples = if(isTRUE(return_samples)){mc_samples}else{ NULL} 
  )
}

power_test2_F <- function(ncp1, k1, ncp2=NULL, k2=0, df1=1, df2=df2, test = "2015", alpha = .05, tuning_power = 1/3, alpha_bound = .05, mc_iter = 0, only_pow = mc_iter < 1, return_samples = FALSE, ...){
  test = match.arg(test,choices = c("2015","Fsum","2014"))
  k1 = ceiling(max(k1,0))
  k2 = ceiling(max(k2,0))
  if(k1 + k2 < 1)
    stop("One of k1 or k2 must be > 0.")
  if(is.null(ncp2) & k2 != 0)
    stop("If k2 is not 0, you must specify ncp2")
  tr = qf(alpha_bound, df1, df2, lower.tail = FALSE)
  ncp_tuning = find_ncp_F_uniroot(tuning_power = tuning_power, alpha = alpha_bound, df1 = df1, df2 = df2)
  ncp1 = abs(ncp1)
  ncp2 = abs(ncp2)
  if(alpha<=0 | alpha>=1)
    stop("Invalid alpha.")
  if(alpha_bound<=0 | alpha_bound>=1)
    stop("Invalid alpha_bound.")
  k = sort(c(k1,k2))
  ncp = c(ncp1, ncp2)[order(c(k1,k2))]
  if(test == "2015"){
    left_tr = NULL
    crit = qnorm(alpha)
    dens_fun = \(x, ncp) d2015_test2_F(x, ncp=ncp, k = sum(k), alpha = alpha_bound, df1=df1, df2=df2)
  }else if(test == "Fsum"){
    left_tr = tr
    crit = get_crit_dtf2(k = sum(k), alpha = alpha, alpha_bound = alpha_bound, ncp = ncp_tuning, df1=df1, df2=df2)
    dens_fun = \(x, ncp) dtf(x, df1=df1, df2=df2, ncp = ncp, tr = qf(1-alpha_bound, df1, df2))
  }else if(test == '2014'){
    left_tr = 0
    crit = qchisq(alpha,2*sum(k), lower.tail = FALSE)
    dens_fun = \(x, ncp) d2014_test2_F(x, ncp, alpha = alpha_bound, df1=df1, df2=df2)
  }else{
    stop("Invalid test: must be one of '2014', '2015', or 'Fsum'")
  }
  if(k[1]==0){
    dsum <- dsumf1(\(x)  dens_fun(x, ncp[2]), k1 = k[2], left_tr = left_tr, nx=2^13, sc=30, ...)
  }else{
    dsum = dsumf2(
      f1 = \(x) dens_fun(x, ncp[1]), 
      f2 = \(x) dens_fun(x, ncp[2]),
      k1 = k[1],
      k2 = k[2],
      left_tr = left_tr,
      nx=2^13, sc=30,
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
    if(df1 == 1){
      z0 = replicate(mc_iter, c(rtfoldedt(k[1],ncp=sqrt(ncp[1]),df=df2),rtfoldedt(k[2],ncp=sqrt(ncp[2]),df=df2))^2)
    }else{
      z0 = replicate(mc_iter, c(rtf(k[1],ncp=ncp[1],df1=df1,df2=df2),rtchisq(k[2],ncp=ncp[2],df1=df1,df2=df2)))
    }
    if(test == "2015"){
      mc_samples = apply(z0, 2, \(z) pcurve_F(z)[,1,2])[2,]
      mc_pow = mean(mc_samples<alpha)
    }else if(test == "Fsum"){
      mc_samples=colSums(z0)
      mc_pow = mean(mc_samples<crit)
    }else if(test == '2014'){
      mc_samples = apply(z0, 2, \(z) pcurve_F(z)[,1,1])[2,]
      mc_pow = mean(mc_samples<alpha)
    }
  }
  list(
    test = test,
    pow = pow,
    ncp = ncp,
    k = k,
    mc_iter = mc_iter,
    mc_pow = mc_pow,
    mc_se = sqrt(mc_pow * (1-mc_pow) / mc_iter),
    mc_samples = if(isTRUE(return_samples)){mc_samples}else{ NULL} 
  )
}
