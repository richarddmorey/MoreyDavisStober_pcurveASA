teststat_string = function(vals, types = 'z'){
  glue::glue('{types}={vals}') |> paste(collapse='\n')
}

pcurve_functions0 = function(
    url = 'https://p-curve.com/app4/pcurve_app4.10.r',
    path = here::here('R_utility/Simonsohn/pcurve_app4.10.r'),
    md5sum = 'b7b8cb6e5aee8fc25263b5c0493146a9',
    check_md5 = TRUE, from_web = FALSE)
{
  
  if(from_web){
    path = tempfile(fileext = '.R')
    download.file(url, destfile = path, mode = "wb")
  }
  
  ## Error if the MD5 hash of the file does not equal the one that
  ## I obtained when I downloaded the file on 1 Jul 2024
  if(check_md5)
    stopifnot(tools::md5sum(path) == md5sum)
  
  pkg_list = renv::dependencies(path)
  pkg_list$Package |>
    purrr::walk(\(pkg){
      if(!require(pkg, character.only = TRUE, quietly = TRUE)){
        stop('Package ', pkg, ' not installed.')
      }
    })
  
  readLines(path) |>
    paste(collapse = '\n')
}

pcurve_functions = memoise::memoise(pcurve_functions0)

do_online_pcurve_analysis = function(string){
  
  new_env = new.env()
  source(textConnection(pcurve_functions()),local = new_env)
  
  td = tempfile(pattern = 'dir')
  stopifnot(dir.create(td))
  tf_stats = tempfile(tmpdir = td, pattern = 'teststats',fileext = '.txt')
  cat(string, file = tf_stats)
  
  withr::with_dir(td, new_env$pcurve_app(basename(tf_stats), td))
  
  analysis_file = dir(td, pattern = '^STOUFFER_.*\\.txt$')
  analysis_results = file.path(td, analysis_file) |> readLines() |> as.numeric()
  
  names(analysis_results) = c("ktot", "ksig", "khalf", "Zppr", "p.Zppr", "Zpp33", "p.Zpp33", "Zppr.half", "p.Zppr.half", "Zpp33.half", "p.Zpp33.half")
  
  return(list(dir = td, results = analysis_results))
}

#' Find the noncentrality parameter such that the chi-squared RV has 
#' some probability of being above the alpha critical value.
#'
#' This function uses optimize(), which allows the penalizing of solutions
#' that would lead to numerical problems.
#'
#' @param tuning_power probability desired for being above the lower bound
#' @param alpha alpha used to determining lower bound (under the null)
#' @param df Degrees of freedom for the chi-squared RVs
#' @param lb lower bound for the ncp search 
#' @param ub upper bound for the ncp search 
#' @param ... Extra arguments for optimize()
find_ncp_chi2_opt = Vectorize(function(tuning_power = 1/3, alpha = .05, df = 1, lb = 0, ub, ...){
  if(tuning_power <  alpha) stop('tuning_power must be >= alpha.')
  if(tuning_power == alpha) return(0)
  crit.z2 = qchisq( alpha, df, lower.tail = FALSE )
  f = function( ncp0 ){
    a = pchisq( crit.z2, df, ncp0, lower.tail = FALSE) 
    if(a>tuning_power) return(.Machine$double.xmax)
    (a-tuning_power)^2
  }
  ncp0 = optimize(f,c(lb,ub),tol = .Machine$double.eps^0.8,...)
  return(ncp0$minimum)
},"df")

#' Find the noncentrality parameter such that the chi-squared RV has 
#' some probability of being above the alpha critical value.
#'
#' This function will subsequently call find_ncp_chi2_optim (which uses 
#' optimize instead of uniroot) if the
#' search produces a value that would lead to a probability greater than
#' tuning_power (which would cause numerical problems). However, I've 
#' generally found uniroot to give better solutions so we do a first pass
#' with uniroot.
#'
#' @param tuning_power probability desired for being above the lower bound
#' @param alpha alpha used to determining lower bound (under the null)
#' @param df Degrees of freedom for the chi-squared RVs
#' @param ub upper-bound for the ncp search 
#' @param eps value to subtract off the solution for the lower bound used for the 
#'            subsequent find_ncp_chi2_optim call
#' @param ... Extra arguments for uniroot
find_ncp_chi2_uniroot = Vectorize(function(tuning_power = 1/3, alpha = .05, df = 1, ub = 100, eps = .01, ...){
  if(tuning_power <  alpha) stop('tuning_power must be >= alpha.')
  if(tuning_power == alpha) return(0) # Obviously the NCP for alpha is 0
  crit.z2 = qchisq( alpha, df, lower.tail = FALSE )
  f = function( ncp0 ){
    a = pchisq( crit.z2, df, ncp0, lower.tail = FALSE) 
    a-tuning_power
  }
  r = uniroot(f,c(0,ub),...)
  lambda0 = r$root
  if(r$f.root>0)
    lambda0 = find_ncp_chi2_opt(tuning_power = tuning_power, alpha = alpha, df = df, lb = lambda0 - eps, ub = lambda0)
  return(lambda0)
},"df")


#' Compute the 2014 and 2015 p curve tests for z values.
#'
#' @param z vector of z values
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
pcurve_Z <- function(z, tuning_power = 1/3, alpha_bound = .05, log.p = FALSE, recompute.power = FALSE){
  pcurve_X2(
    z^2,
    df = 1,
    tuning_power = tuning_power,
    alpha_bound = alpha_bound,
    log.p = log.p,
    recompute.power = recompute.power)
}

#' Compute the 2014 and 2015 p curve tests for chi-squared values.
#'
#' @param x2 vector of chisq values
#' @param df degrees of freedom for the chi-squared values
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
pcurve_X2 <- function(x2, df=1, tuning_power = 1/3, alpha_bound = .05, log.p = FALSE, recompute.power = FALSE){
  log.p = isTRUE(log.p)
  crit.x2 = qchisq(alpha_bound,df,lower.tail = FALSE)
  lambda0 = find_ncp_chi2_uniroot(tuning_power, alpha_bound, df = df)
  if(isTRUE(recompute.power)){
    tuning_power2 = pchisq(crit.x2, df, lambda0, lower.tail = FALSE)
  }else{
    tuning_power2 = tuning_power
  }
  log_p = pchisq(x2, df, log.p = TRUE, lower.tail = FALSE)
  if(any(log_p >= log(alpha_bound)))
    stop('Some p values were nonsignificant. Please remove them.')
  n = length(x2)
  # Log p values for test 1
  log_pp = log_p - log(alpha_bound)
  # Log p values for test 2
  log_qq = -log(tuning_power2) +
    log(1-tuning_power2) + 
    log( expm1(pchisq(x2, df, lambda0, log.p = TRUE) - log(1-tuning_power2))) 
  
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

