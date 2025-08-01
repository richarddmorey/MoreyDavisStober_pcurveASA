#' Find the noncentrality parameter such that the chi-squared RV has 
#' some probability of being above the alpha critical value.
#'
#' This function uses optimize(), which allows the penalizing of solutions
#' that would lead to numerical problems.
#'
#' @param tuning_power probability desired for being above the lower bound
#' @param alpha alpha used to determining lower bound (under the null)
#' @param df1 Numerator degrees of freedom for the F RVs
#' @param df2 Denominator degrees of freedom for the F RVs
#' @param lb lower bound for the ncp search 
#' @param ub upper bound for the ncp search 
#' @param ... Extra arguments for optimize()
find_ncp_F_opt = Vectorize(function(tuning_power = 1/3, alpha = .05, df1 = 1, df2, lb = 0, ub, ...){
  if(tuning_power <  alpha) stop('tuning_power must be >= alpha.')
  if(tuning_power == alpha) return(0)
  crit.F = qf( alpha, df1, df2, lower.tail = FALSE )
  f = function( ncp0 ){
    a = pf( crit.F, df1, df2, ncp0, lower.tail = FALSE) 
    if(a>tuning_power) return(.Machine$double.xmax)
    (a-tuning_power)^2
  }
  ncp0 = optimize(f,c(lb,ub),tol = .Machine$double.eps^0.8,...)
  return(ncp0$minimum)
},"df1","df2")

#' Find the noncentrality parameter such that the F RV has 
#' some probability of being above the alpha critical value.
#'
#' This function will subsequently call find_ncp_F_optim (which uses 
#' optimize instead of uniroot) if the
#' search produces a value that would lead to a probability greater than
#' tuning_power (which would cause numerical problems). However, I've 
#' generally found uniroot to give better solutions so we do a first pass
#' with uniroot.
#'
#' @param tuning_power probability desired for being above the lower bound
#' @param alpha alpha used to determining lower bound (under the null)
#' @param df1 Numerator degrees of freedom for the F RVs
#' @param df2 Denominator degrees of freedom for the F RVs
#' @param ub upper-bound for the ncp search 
#' @param eps value to subtract off the solution for the lower bound used for the 
#'            subsequent find_ncp_F_optim call
#' @param ... Extra arguments for uniroot
find_ncp_F_uniroot = Vectorize(function(tuning_power = 1/3, alpha = .05, df1 = 1, df2, ub = 100, eps = .01, ...){
  if(tuning_power <  alpha) stop('tuning_power must be >= alpha.')
  if(tuning_power == alpha) return(0) # Obviously the NCP for alpha is 0
  crit.F = qf( alpha, df1, df2, lower.tail = FALSE )
  f = function( ncp0 ){
    a = pf( crit.F, df1, df2, ncp0, lower.tail = FALSE) 
    a-tuning_power
  }
  r = uniroot(f,c(0,ub),...)
  lambda0 = r$root
  if(r$f.root>0)
    lambda0 = find_ncp_F_opt(tuning_power = tuning_power, alpha = alpha, df1 = df1, df2 = df2, lb = lambda0 - eps, ub = lambda0)
  return(lambda0)
},"df1","df2")

#' Compute the 2014 and 2015 p curve tests for t values.
#'
#' @param ts vector of z values
#' @param df degrees of freedom for the t values
#' @param tuning_power Power to use for test 2
#' @param alpha_bound alpha for bound for p hacking/pub bias assumption
#' @param log.p (logical) return the log of the p values?
#' @param recompute.power (logical) Recompute the power after computing the noncentrality 
#'        parameter? This can increase numberical stability.
#'
#' @return
#' @export
#'
#' @examples
pcurve_t <- function(ts, df, tuning_power = 1/3, alpha_bound = .05, log.p = FALSE, recompute.power = FALSE){
  pcurve_F(
    ts^2,
    df1 = 1,
    df2 = df,
    tuning_power = tuning_power,
    alpha_bound = alpha_bound,
    log.p = log.p,
    recompute.power = recompute.power)
}


#' Compute the 2014 and 2015 p curve tests for F values.
#'
#' @param Fs vector of F values
#' @param df1 Numerator degrees of freedom for the F values
#' @param df2 Denominator degrees of freedom for the F values
#' @param tuning_power Power to use for test 2
#' @param alpha_bound alpha for bound for p hacking/pub bias assumption
#' @param log.p (logical) return the log of the p values?
#' @param recompute.power (logical) Recompute the power after computing the noncentrality 
#'        parameter? This can increase numberical stability.
#'
#' @return
#' @export
#'
#' @examples
pcurve_F <- function(Fs, df1=1, df2, tuning_power = 1/3, alpha_bound = .05, log.p = FALSE, recompute.power = FALSE){
  log.p = isTRUE(log.p)
  crit.F = qf(alpha_bound,df1,df2,lower.tail = FALSE)
  lambda0 = find_ncp_F_uniroot(tuning_power, alpha_bound, df1 = df1, df2 = df2)
  if(isTRUE(recompute.power)){
    tuning_power2 = pf(crit.F, df1, df2, lambda0, lower.tail = FALSE)
  }else{
    tuning_power2 = tuning_power
  }
  log_p = pf(Fs, df1, df2, log.p = TRUE, lower.tail = FALSE)
  if(any(log_p >= log(alpha_bound)))
    stop('Some p values were nonsignificant. Please remove them.')
  n = length(Fs)
  # Log p values for test 1
  log_pp = log_p - log(alpha_bound)
  # Log p values for test 2
  log_qq = -log(tuning_power2) +
    log(1-tuning_power2) + 
    log( expm1(pf(Fs, df1, df2, lambda0, log.p = TRUE) - log(1-tuning_power2))) 
  
  # Test 1, 2014
  test1_2014_teststat = -2 * sum(log_pp)
  test1_2014_pvalue = pchisq(test1_2014_teststat, 2*n, log.p = log.p, lower.tail = FALSE)
  
  # Test 2, 2014
  test2_2014_teststat = -2 * sum(log_qq)
  test2_2014_pvalue = pchisq(test2_2014_teststat, 2*n, log.p = log.p, lower.tail = FALSE)
  
  # Test 1, 2015
  test1_2015_teststat = sum(qnorm(log_pp, log.p = TRUE)/sqrt(n))
  test1_2015_pvalue = pnorm(test1_2015_teststat, log.p = log.p)
  
  # Test 2, 2015
  test2_2015_teststat = sum(qnorm(log_qq, log.p = TRUE)/sqrt(n))
  test2_2015_pvalue = pnorm(test2_2015_teststat, log.p = log.p)
  
  if(isTRUE(recompute.power))
    message('Recomputed power was ', as.character(tuning_power2))
  
  x = c(
    test1_2014_teststat = test1_2014_teststat,
    test1_2014_pvalue = test1_2014_pvalue,
    test2_2014_teststat = test2_2014_teststat,
    test2_2014_pvalue = test2_2014_pvalue, 
    test1_2015_teststat = test1_2015_teststat,
    test1_2015_pvalue = test1_2015_pvalue, 
    test2_2015_teststat = test2_2015_teststat,
    test2_2015_pvalue = test2_2015_pvalue
  )
  
  z = array(x,dim = c(2,2,2))
  dimnames(z) = list(statistic = c("test","p"), test = c(1,2), pub = c("2014","2015"))
  return(z)
}

