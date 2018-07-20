#' Approximation of Bessel's functions
#'
#' @param x argumet of Bessel's function
#'
#' @return value
#' @export
#'
#' @seealso \code{\link{BesselI1}} and \code{\link{BesselK0}} and \code{\link{BesselK1}}
#' @examples
#' x=0.5
#' BesselI0(x)
BesselI0 = function(x){
  absx=abs(x)
  ans=c()
  if (absx < 3.75) {
    y = x/3.75
    y = y*y
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))))
  } else {
    y=3.75/absx
    ans=(exp(absx)/sqrt(absx))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))))))))
  }
  return(ans)
}
#' Approximation of Bessel's functions
#'
#' @param x argumet of Bessel's function
#'
#' @return value
#' @export
#'
#' @seealso \code{\link{BesselI0}} and \code{\link{BesselK0}} and \code{\link{BesselK1}}
#' @examples x=0.5
#' BesselI1(x)
BesselI1 = function(x){
  absx = abs(x)
  ans=c()
  if (absx<3.75) {
    y=x/3.75
    y=y*y
    ans=absx*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))))
  } else {
    y=3.75/absx
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2))
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))))
    ans= ans * (exp(absx)/sqrt(absx))
  }
  if(x<0)  ans = -ans
  return(ans)
}

#' Approximation of Bessel's functions
#'
#' @param x argumet of Bessel's function
#'
#' @return value
#' @export
#'
#' @seealso \code{\link{BesselI1}} and \code{\link{BesselI0}} and \code{\link{BesselK1}}
#' @examples x=0.5
#' BesselK0(x)
BesselK0 = function(x){
  ans=c()
  if (x <= 2.0) {
    y=x*x/4.0
    ans=(-log(x/2.0)*BesselI0(x))-0.57721566+0.42278420*y+0.23069756*y*y+0.03488590*y*y*y+0.00262698*y*y*y*y+0.00010750*y*y*y*y*y+0.00000740*y*y*y*y*y*y
  } else {
    y=2.0/x
    ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2+y*(-0.251540e-2+y*0.53208e-3))))))
  }
  return(ans)
}


#' Approximation of Bessel's functions
#'
#' @param x argumet of Bessel's function
#'
#' @return value
#' @export
#'
#' @seealso \code{\link{BesselI1}} and \code{\link{BesselI0}} and \code{\link{BesselK0}}
#' @examples
#' x=0.5
#' BesselI0(x)
BesselK1 = function(x){
  ans=c()
  if (x <= 2.0) {
    y=x*x/4.0
    ans=(log(x/2.0)*BesselI1(x))+(1.0/x)*(1.0+y*(0.15443144+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1+y*(-0.110404e-2+y*(-0.4686e-4)))))))
  } else {
    y=2.0/x
    ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2+y*(0.325614e-2+y*(-0.68245e-3)))))))
  }
  return(ans);
}
