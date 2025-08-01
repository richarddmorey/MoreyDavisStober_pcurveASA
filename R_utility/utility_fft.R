moments = function(f1, interval = c(-Inf,Inf)){
  Ex = integrate(\(x) x*f1(x), lower = min(interval), upper = max(interval))[[1]] 
  Ex2= integrate(\(x) x^2*f1(x), lower = min(interval), upper = max(interval))[[1]] 
  Sx = sqrt(Ex2 - Ex^2)
  c(Ex = Ex, Sx = Sx)
}

convolve_rep = function (x, k) 
{
  n <- length(x)
  p <- c(x, rep.int(0, (k-1)*(n-1)))
  s <- sum(p)
  d <- fft(fft(p/s)^k, inverse = TRUE)
  Re(d)*s/n
}

## The functions below need to be calibrated to the distributions passed to it, 
## so that they capture enough of the space for accurate estimates. In particular,
## a (lower bound), b (upper bound) and nx (number of points) -- and maybe
## sc, which is a scale factor that attempts to estimate good bounds from the mean
## and the variance -- need to be set well.

## https://stackoverflow.com/a/77482400/1129889
## https://stackoverflow.com/a/50324164/1129889
dsumf1 <- function(f1, k1, nx = 2^10, a=NULL, b=NULL, sc = 7, recenter = TRUE, left_tr = FALSE, moments_interval = c(-Inf,Inf)){
  m1 = moments(f1, moments_interval)
  if(is.null(b))
    b = m1[1] + sc*m1[2]
  if(is.null(a))
    a = m1[1] - sc*m1[2]
  x1 = seq(a, b, length.out = nx)
  cv = convolve_rep(f1(x1), k1)
  cs = cumsum(cv)/sum(cv)
  if(!recenter) return(list(x = 1:length(cv), cv = cv, cs = cs))
  i = 1:length(cv)
  Ei = sum(i*cv)/sum(cv)
  Ei2 = sum(i^2*cv)/sum(cv)
  Si = sqrt(Ei2 - Ei^2)
  x0 = (i - Ei)/Si * m1[2]*sqrt(k1) + k1*m1[1]
  if(!is.null(left_tr)){
    cv[x0<a] = 0
  }
  cv = cv/sum(cv)/(x0[2]-x0[1])
  list(
    x = x0,
    cs = cs,
    cv = cv
  )
}


dsumf2 <- function(f1, f2, k1, k2, nx = 2^10, sc = 7, left_tr = NULL){
  m1 = moments(f1)
  m2 = moments(f2)
  sum_Ex = m1[1] + m2[1]
  sum_Sx = sqrt(m1[2]^2*k1 + m2[2]^2*k2)
  b = sum_Ex + sc*sum_Sx
  a = ifelse(is.null(left_tr), sum_Ex - sc*sum_Sx, left_tr)
  cv = convolve(
    dsumf1(f1, k1, nx = nx, sc = sc, a = a, b = b, left_tr = !is.null(left_tr), recenter = FALSE)$cv,
    rev(dsumf1(f2, k2, nx = nx, sc = sc, a = a, b = b, left_tr = !is.null(left_tr), recenter = FALSE)$cv),
    type = 'open'
  )
  i = 1:length(cv)
  Ei = sum(i*cv)/sum(cv)
  Ei2 = sum(i^2*cv)/sum(cv)
  Si = sqrt(Ei2 - Ei^2)
  x0 = (i - Ei)/Si * sqrt(m1[2]^2*k1 + m2[2]^2*k2) + k1*m1[1] + k2*m2[1]
  if(!is.null(left_tr)){
    cv[x0<left_tr] = 0
  }
  cv = cv/sum(cv)/(x0[2]-x0[1])
  cs = cumsum(cv)/sum(cv)
  list(
    x = x0,
    cs = cs,
    cv = cv
  )
}
