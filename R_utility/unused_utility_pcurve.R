## Functions erf, erfc, erfi, and erfci are for code clarity

#' The Gauss error function
#'
#' @param x A real value
#'
#' @return Returns a probability
#' @export
erf = function(x){
  2 * pnorm(x * sqrt(2)) - 1  
}


#' The complementary Gauss error function
#'
#' @param x A real value
#'
#' @return Returns a probability
#' @export
erfc = function(x){
  1-erf(x)
}

#' The inverse of the Gauss error function
#'
#' @param p A probability 
#'
#' @return Returns a real value
#' @export
erfi = function(p){
  qnorm((p + 1)/2)/sqrt(2)
}

#' The inverse of the complementary Gauss error function
#'
#' @param p 
#'
#' @return
#' @export
erfci = function(p){
  erfi(1-p)
}

power_2015_test1_integrate <- Vectorize(function(mu, k, alpha = .05, df = 1){
  crit = qnorm(alpha)
  ex = integrate(\(x) x*d2015_test1(x,mu = mu, k = k, df = df),lower = -Inf, upper = Inf)[[1]]
  ex2 = integrate(\(x) x^2*d2015_test1(x,mu = mu, k = k, df = df),lower = -Inf, upper = Inf)[[1]]
  vx = ex2 - ex^2
  pnorm(crit, k*ex,sqrt(k)*sqrt(vx))
},c("mu","k"))

power_2015_test1_integrate_general <- function(mu, alpha = .05, df = 1){
  crit = qnorm(alpha)
  n = length(mu)
  ex = sapply(mu, \(mu) integrate(\(x) x*d2015_test1(x,mu = mu, k = k, df = df),lower = -Inf, upper = Inf)[[1]])
  ex2 = sapply(mu, \(mu) integrate(\(x) x^2*d2015_test1(x,mu = mu, k = k, df = df),lower = -Inf, upper = Inf)[[1]])
  vx = ex2 - ex^2
  pnorm(crit, sum(ex), sqrt(sum(vx)))
}

# Adapted from Simonsohn's code
getncp.c =function(df, power, alpha)   {      
  xc=qchisq(p=alpha, df=df, lower.tail = FALSE) 
  error = function(ncp_est, power, x, df)      pchisq(x, df = df, ncp = ncp_est) - (1-power)   
  return(uniroot(error, c(0, 1000), x = xc, df = df, power=power)$root)   }
# End Simonsohn's code

#' Perform P curve analyses of Morey (in progress) for TWO Z values
#'
#' @param Z Z values for the P curve analysis
#' @param tuning.power The probability that the chi^2 random variable takes on a value greater than the alpha-level critical value of the central chi^2 distribution
#' @param alpha The probability that a central chi^2 value takes on a value greater than the critical value
#' @param upper.bound.ncp The upper bound for the search for the noncentrality parameter
#' @param cat.results Print the statistic to the console for easy pasting into the online app
#' @param method Method for computing p values ("Fisher", "likelihood", or "both")
#' @param MC.iterations For "likelihood", how many Monte Carlo iterations to perform. More iterations gives better precision. 
#' @param ... arguments to pass to related functions
#' 
#' @return Returns a 3D array of statistics and p values.
#' @export
#'
#' @examples
#' Z = c(1.96, 3.6)
#' pcurve_RDM_Z2(Z, 1/3, .05, method = "both")
pcurve_RDM_Z2 = function(Z, tuning.power, alpha, upper.bound.ncp = 30, cat.results = FALSE, method = "Fisher", MC.iterations = 10000, ...){
  
  # defines Z, Z2, k, and p.values
  list2env( pvalues_Z(Z, alpha), envir = environment(), ... )
  
  stopifnot(method %in% c("Fisher", "likelihood", "both"))
  stopifnot(length(p.values)==2)
  
  result = array(NA,c(2,2,2))
  
  if(method %in% c("Fisher", "both")){
    stat = sum(log(p.values))
    test1.F = Fisher_combine_p( p.values / alpha)
    test1.F[['stat']] = stat
    p = 1 - sum2logp_cdf_chi2(stat, tuning.power, alpha, ...)
    test2.F = c(stat = stat, p = p)
    result[,,1] = c(test1.F, test2.F)
  }
  if(method %in% c("likelihood", "both")){
    result[,,2] = pcurve_RDM_likelihood_Z(Z, tuning.power, alpha, upper.bound.ncp, MC.iterations, ...)
  }
  
  dimnames(result) = list(c("Statistic","p"),
                          c("Test 1", "Test 2"),
                          c("Fisher","(log) Likelihood"))
  
  # cat() values to the console so that they
  # can be easily copy/pasted into the online app
  if(cat.results)
    cat(paste("Z =",Z),"\n", sep = "\n")
  
  return(result)
}


loglike_moments_integrate_chi2 = function(upper = 10, moment = 1, tuning.power, alpha, lambda0 = NULL, upper.bound.ncp = 30, ...){
  
  if( is.null(lambda0) )
    lambda0 = find_ncp_chi2(tuning.power, alpha, upper.bound.ncp)
  
  lower = likelihood_ratio_p_chi2(alpha, tuning.power, alpha, log = TRUE, upper.bound.ncp = upper.bound.ncp)
  
  integrate(function(x, alpha, ...){
    x^moment*log_likelihood_ratio_density(x, tuning.power = tuning.power, alpha = alpha, lambda0 = lambda0, upper.bound.ncp = upper.bound.ncp, ...)
  }, lower = lower, upper = upper, alpha = alpha, ...)[[1]]
}


log_likelihood_ratio_inverse = Vectorize(function(llr, tuning.power, alpha, lambda0 = NULL, upper.bound.ncp = 30, optim.lower = exp(-30) ){
  
  if( is.null(lambda0) )
    lambda0 = find_ncp_chi2(tuning.power, alpha, upper.bound.ncp)
  
  fn = function(p){
    res = (p_density_chi2(p, tuning.power, alpha, lambda0 = lambda0, log = TRUE, upper.bound.ncp = upper.bound.ncp) + 
             log(alpha) + 
             -llr)^2
  }
  gr = function(p){
    2 * (p_density_chi2(p, tuning.power, alpha, lambda0 = lambda0, log = TRUE, upper.bound.ncp = upper.bound.ncp) + 
           log(alpha) + 
           -llr) * 
      d_likelihood_ratio_p_chi2(p, tuning.power, alpha, upper.bound.ncp)  
  }
  
  optim(par = 0.025, f = fn , gr = gr, method = "L-BFGS-B", lower = optim.lower, upper = alpha)$par
  #uniroot(solve.func, interval = c(uniroot.lower,alpha))$root
},"llr")

sum2llr_density_chi2 = Vectorize(function(sumllr, tuning.power, alpha, lambda0 = NULL, null = FALSE, log = FALSE, upper.bound.ncp = 30, ...){
  
  lower1 = likelihood_ratio_p_chi2(alpha, tuning.power, alpha, log = TRUE, upper.bound.ncp = upper.bound.ncp)
  
  if(sumllr < 2 * lower1) return(ifelse(log, -Inf, 0))
  
  if(null & !is.null(tuning.power))
    warning("null was true, but tuning.power was not set to NULL. Null assumed.")
  
  if(is.null(lambda0))
    lambda0 = find_ncp_chi2(tuning.power, alpha, upper.bound.ncp)
  
  # numerical convolution
  fv_convo = function(y, z, ...){
    x = z - y
    exp(
      log_likelihood_ratio_density(x, tuning.power, alpha, lambda0, log = TRUE, null = null, ...) +
        log_likelihood_ratio_density(y, tuning.power, alpha, lambda0, log = TRUE, null = null, ...)
    )
  }
  
  res = log(integrate(fv_convo, lower = lower1,  
                      upper = sumllr - lower1, z = sumllr, ...)[[1]])
  
  if(log){
    return(res)
  }else{
    return(exp(res))
  }
}, "sumllr")

sum2llr_cdf_chi2 = Vectorize(function(sumllr, tuning.power, alpha, lambda0 = NULL, null = FALSE, upper.integration.limit = 10, ...){
  
  lower1 = likelihood_ratio_p_chi2(alpha, tuning.power, alpha, log = TRUE, upper.bound.ncp = upper.bound.ncp)
  
  if(sumllr > upper.integration.limit) return(1)
  if(sumllr < 2*lower1) return(0)
  
  integrate(sum2llr_density_chi2, 2*lower1, sumllr, tuning.power = tuning.power, alpha = alpha, lambda0 = lambda0, null = null, ...)[[1]]  
}, "sumllr")

sum2llr_quantile_chi2 = Vectorize(function(p, tuning.power, alpha, lambda0 = NULL, null = FALSE, upper.integration.limit = 10, ...){
  
  lower1 = likelihood_ratio_p_chi2(alpha, tuning.power, alpha, log = TRUE, upper.bound.ncp = upper.bound.ncp)
  
  if(p < 0 | p > 1 ) return(NA)
  if(p == 0) return(2*lower1)
  if(p == 1) return(Inf)
  
  optimize(function(z, ...) (p - sum2llr_cdf_chi2(z, tuning.power, alpha, lambda0, null = null, upper.integration.limit, ...))^2,
           c(2*lower1, upper.integration.limit))$minimum
}, "p")

#' Derivative of the log likelihood ratio for a set of p values
#'
#' @param p p values for the Z P curve analysis
#' @param tuning.power The probability that the chi^2 random variable takes on a value greater than the alpha-level critical value of the central chi^2 distribution
#' @param alpha The probability that a central chi^2 value takes on a value greater than the critical value
#' @param upper.bound.ncp The upper bound for the search for the noncentrality parameter
#'
#' @return Returns the derivative of the log likelihood ratio for \eqn{\lambda = \lambda_0} to 
#' \eqn{\lambda = 0}, computed at p.
#' @export
#'
#' @examples
d_likelihood_ratio_p_chi2 = function(p, tuning.power, alpha, upper.bound.ncp = 30){
  
  lambda0 = find_ncp_chi2(tuning.power, alpha, upper.bound = upper.bound.ncp)
  
  t = erfci(p)
  t2 = t^2
  
  result = log(lambda0) + .5 * log(pi) + log(hypergeo::genhypergeo(NULL,1+1/2,lambda0*t2/2)) - log(hypergeo::genhypergeo(NULL,1/2,lambda0*t2/2)) + log(t) + t2 
  return(-exp(result))
}

#' CDF of sum of two log p values 
#'
#' @param sumlogp Value at which to compute the density
#' @param tuning.power The probability that the chi^2 random variable takes on a value greater than the alpha-level critical value of the central chi^2 distribution
#' @param alpha The probability that a central chi^2 value takes on a value greater than the critical value
#' @param lambda0 If this is not null, then the function will not search for the noncentrality parameter (saves time)
#' @param log If TRUE, return the log density
#' @param lower.integration.limit The lower bound for quadrature algorithm
#' @param ... Values to pass to logp_density_chi2
#'
#' @return
#' @export
#'
#' @examples
#' lambda0 = find_ncp_chi2(1/3, 0.05)
#' 
#' Z2 = qchisq(runif(2000, 2/3, 1), 1, lambda0)
#' dim(Z2) = c(2, 1000)
#' 
#' ps = 1 - pchisq(Z2, 1)
#' sums = colSums(log(ps))
#' plot(ecdf(sums))
#' 
#' curve(sum2logp_cdf_chi2(x, 1/3, .05, lambda0),
#'  min(sums), 2*log(.05)-.0001, col = "red", add = TRUE)
sum2logp_cdf_chi2 = Vectorize(function(sumlogp, tuning.power, alpha, lambda0 = NULL, null = FALSE, lower.integration.limit = -30, ...){
  
  if(null) return(1 - pchisq(-2*(sumlogp - 2*log(alpha)), 4))
  if(sumlogp < lower.integration.limit) return(0)
  
  integrate(sum2logp_density_chi2, lower.integration.limit, sumlogp, tuning.power = tuning.power, alpha = alpha, lambda0 = lambda0, ...)[[1]]  
}, "sumlogp")

#' Density of sum of two log p values
#'
#' @param sumlogp Value at which to compute the density
#' @param tuning.power The probability that the chi^2 random variable takes on a value greater than the alpha-level critical value of the central chi^2 distribution
#' @param alpha The probability that a central chi^2 value takes on a value greater than the critical value
#' @param lambda0 If this is not null, then the function will not search for the noncentrality parameter (saves time)
#' @param log If TRUE, return the log density
#' @param upper.bound.ncp The upper bound for the search for the noncentrality parameter
#' @param ... Values to pass to logp_density_chi2
#'
#' @return Returns the (log) density at sum(log(p))
#' @export
#'
#' @examples
#' lambda0 = find_ncp_chi2(1/3, 0.05)
#' 
#' Z2 = qchisq(runif(2000, 2/3, 1), 1, lambda0)
#' dim(Z2) = c(2, 1000)
#' 
#' ps = 1 - pchisq(Z2, 1)
#' sums = colSums(log(ps))
#' hist(sums, freq = FALSE)
#' 
#' curve(sum2logp_density_chi2(x, 1/3, .05, lambda0),
#'  min(sums), 2*log(.05)-.0001, col = "red", add = TRUE)
sum2logp_density_chi2 = Vectorize(function(sumlogp, tuning.power, alpha, lambda0 = NULL, null = FALSE, log = FALSE, upper.bound.ncp = 30, ...){
  
  if(sumlogp > 2*log(alpha)) return(ifelse(log, -Inf, 0))
  
  if(null & !is.null(tuning.power))
    warning("null was true, but tuning.power was not set to NULL. Null assumed.")
  
  if(null){
    res = dchisq(-2 * (sumlogp - 2*log(alpha)), 4, log = TRUE) + log(2)
  }else{
    
    if(is.null(lambda0))
      lambda0 = find_ncp_chi2(tuning.power, alpha, upper.bound.ncp)
    
    # numerical convolution
    fv_convo = function(y, z, ...){
      x = z - y
      exp(
        logp_density_chi2(x, tuning.power, alpha, lambda0, log = TRUE, ...) +
          logp_density_chi2(y, tuning.power, alpha, lambda0, log = TRUE, ...)
      )
    }
    
    res = log(integrate(fv_convo, lower = sumlogp - log(alpha), upper = log(alpha), z = sumlogp, ...)[[1]])
  }
  
  if(log){
    return(res)
  }else{
    return(exp(res))
  }
}, "sumlogp")



#' Quantile function of sum of two log p values
#'
#' @param p 
#' @param tuning.power The probability that the chi^2 random variable takes on a value greater than the alpha-level critical value of the central chi^2 distribution
#' @param alpha The probability that a central chi^2 value takes on a value greater than the critical value
#' @param lambda0 If this is not null, then the function will not search for the noncentrality parameter (saves time)
#' @param log If TRUE, return the log density
#' @param lower.integration.limit The lower bound for quadrature algorithm
#' @param ... Values to pass to logp_density_chi2
#'
#' @return
#' @export
#'
#' @examples
#' lambda0 = find_ncp_chi2(1/3, 0.05)
#' 
#' Z2 = qchisq(runif(2000, 2/3, 1), 1, lambda0)
#' dim(Z2) = c(2, 1000)
#' 
#' ps = 1 - pchisq(Z2, 1)
#' sums = colSums(log(ps))
#' plot(ecdf(sums))
#' 
#' median.value = sum2logp_quantile_chi2(.5, 1/3, .05, lambda0)
#' 
#' abline(h = .5, col = "red")
#' abline(v = median.value, col = "red")
#' points(median.value, .5, col = "red", pch = 19)
sum2logp_quantile_chi2 = Vectorize(function(p, tuning.power, alpha, lambda0 = NULL, null = FALSE, lower.integration.limit = -30, ...){
  
  if(p < 0 | p > 1 ) return(NA)
  if(p == 0) return(-Inf)
  if(p == 1) return(2*log(alpha))
  
  if(null) return( -pchisq(1 - p, 4)/2 + 2*log(alpha) )
  
  optimize(function(z, ...) (p - sum2logp_cdf_chi2(z, tuning.power, alpha, lambda0, null = FALSE, lower.integration.limit, ...))^2,
           c(lower.integration.limit, 2*log(alpha)))$minimum
}, "p")

log_likelihood_ratio_density = function(llr, tuning.power, alpha, lambda0 = NULL, null = FALSE, upper.bound.ncp = 30, log = FALSE, ...){
  
  if( is.null(lambda0) )
    lambda0 = find_ncp_chi2(tuning.power, alpha, upper.bound.ncp)
  
  pp = log_likelihood_ratio_inverse(llr, tuning.power, alpha, lambda0 = lambda0, ...)
  
  res = -log(alpha) - log(abs(d_likelihood_ratio_p_chi2(pp, tuning.power, alpha)))
  
  if(!null) res = res + llr
  
  if(log){
    return(res)
  }else{
    return(exp(res))
  }
}
