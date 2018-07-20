#' The dimensional drawdown of pumping test calcuated using the stehfest algorithm and aquifer parameters.
#'
#' @title The dimensionles drawdown of well pumping test
#'
#' @param Gwinput list of \code{TD} vector of dimensionless time calculated as aquifer Transmisivity[m2/s] * time[s] / well radius [m] / storativity [-], \code{Cbar} dimensionless welbore storage constant,,\code{Skin} skin effect factor
#' @param ... another agruments
#'
#' @return The dimensionless drawdown calculated for each time of \code{TD}.
#' @export
#'
#' @examples
#' # Artifitial exmaple of unconfined aquifer with following parameters of pumping test
#' # Pumped well discharge rate [m3/s]
#' Q=0.0045
#' # Aquifer Transmisivity [m2/s]
#' Transm=0.012
#' # Well radius [m]
#' Rv=1.2
#' # Storativity [-]
#' STRR=0.0001
#' # Time [s]
#' T = c( 1.88,2.36,2.97,3.74,4.71,5.93,7.46,9.40,11.83,14.89)
#' # Dimesionles time [-]
#' TD=T*Transm/Rv/Rv/STRR
#' # Well drawdown during pumping test [m]
#' S =c(0.006,0.008,0.010,0.012,0.015,0.019,0.024,0.03,0.038,0.047)
#' # Dimensionless drawdown [-]
#' dlS = (S*2*3.14*Transm)/Q
#' # Aquifer parameters
#' Cbar= 1562.858133
#' Skin = 4.976677
#' #calcutating the dimensionless drawdown using stehfest algorithm
#' Kytlice_2018_15_05 = PumpingTestDataLoad(Q=Q,Storativity = STRR,
#' Transmisivity = Transm,Radius = Rv,Time = T,DrawDown = S,WellTestType = "IsotropicMediumBesselR")
#' Gwinput = SetStehfestData(Kytlice_2018_15_05$DimlesTime,10,Cbar,Skin,WellTestType= "IsotropicMediumBesselR")
#' simDD <- DimensionlessDrawDown(Gwinput)
#' # Comparison plot
#' plot(TD,dlS, pch=19,ylab="Dimensionless Drawdown [-]",xlab = "Dimensionless Time [-]")
#' lines(TD,simDD,col="red")
#'
#' @references
#' R. G. Agarwal, R. Al-Hussainy and H. J. Ramey Jr.,
#' “An Investigation of Wellbore Storage and Skin Effect in Unsteady Liquid Flow: I. Analytical Treatment,”
#' SPE Journal, Vol. 10, No. 3, 1970, pp. 279-290.
#'
#' Stehfest, H. 1970. Algorithm 368 numerical inversion of Laplace transforms. Communication of ACM 13, no. 1: 47–49.
#'
#' Walton, W.C. 2006. Aquifer Test Modeling. Boca Raton, Florida: CRC Press, Taylor & Francis Group.
#'
DimensionlessDrawDown <-  function(Gwinput){
  UseMethod("DimensionlessDrawDown")
}

#' @export
DimensionlessDrawDown.default <- function(Gwinput){
  warning(paste("DimensionlessDrawDown does not know how to handle object of class ",
                class(Gwinput),
                "and can only be used on classes IsotropicMediumBesselR, IsotropicMediumBesselApr."))
}

#' @export
DimensionlessDrawDown.IsotropicMediumBesselR <- function(Gwinput){
  N = length(Gwinput$Vi)
  # Walton aproximation of Vi
  # Vi = c(0.08333333333333333,-32.08333333333334,1279.0,-15623.66666666667,84244.1666666666,-236957.5,375911.66666666667,-340071.6666666667,164062.5,-32812.5)
  res=c()
  for(j  in 1:length(Gwinput$TD)){
    timeLT = Gwinput$TD[j];
    sum = 0;
    for(i in 1:N){
      k = (i + i)/2;
      # sum = sum + Vi[i] * AgarwalFunBESSELapr((i * (log(2)/t)),Cbar,Skin);
      sum = sum + Gwinput$Vi[i] * AgarwalFunBesselapr((i * (log(2)/timeLT)),Gwinput$Cbar,Gwinput$Skin);
    }
    res = c(res, log(2) / timeLT * sum);
  }
  res
}

#' @export
DimensionlessDrawDown.IsotropicMediumBesselApr <- function(Gwinput){
  N = length(Gwinput$Vi)
  # Walton aproximation of Vi
  # Vi = c(0.08333333333333333,-32.08333333333334,1279.0,-15623.66666666667,84244.1666666666,-236957.5,375911.66666666667,-340071.6666666667,164062.5,-32812.5)
  res=c()
  for(j  in 1:length(Gwinput$TD)){
    timeLT = Gwinput$TD[j];
    sum = 0;
    for(i in 1:N){
      k = (i + i)/2;
      sum = sum + Vi[i] * AgarwalFunBESSELapr((i * (log(2)/t)),Cbar,Skin);
    }
    res = c(res, log(2) / timeLT * sum);
  }
  res
}


#' Laplace transforms of Well function following the Agarwal  et. al. 1970, based on Bessels functions built in R.
#' Ideal isotropic medium, slightly compressible fluid and ideal radial flow, wellbore storage and skin effect to Transient flow.
#'
#' @title The Laplace transforms of Well function
#'
#' @param par variable of Laplace tramsform
#' @param Cbar dimensionless welbore storage constant
#' @param Skin skin effect factor
#'
#' @return discrete value of Laplace transforms of Well function
#'
#' @references
#' R. G. Agarwal, R. Al-Hussainy and H. J. Ramey Jr.,
#' “An Investigation of Wellbore Storage and Skin Effect in Unsteady Liquid Flow: I. Analytical Treatment,”
#' SPE Journal, Vol. 10, No. 3, 1970, pp. 279-290.
#'
#' Walton, W.C. 2006. Aquifer Test Modeling. Boca Raton, Florida: CRC Press, Taylor & Francis Group.
AgarwalFunBesselR = function(par,Cbar,Skin){
  denom = par*(sqrt(par)*besselK(sqrt(par),1) + Cbar * par * (besselK(sqrt(par),0) + Skin * sqrt(par) * besselK(sqrt(par),1)))
  gwfun = (besselK(sqrt(par),0) + (Skin * sqrt(par) * besselK(sqrt(par),1)))/denom;
  return(gwfun);
}

#' Laplace transforms of Well function following the Agarwal  et. al. 1970, based on approximation of Bessels functions built in package.
#' Ideal isotropic medium, slightly compressible fluid and ideal radial flow, wellbore storage and skin effect to Transient flow.
#'
#' @title The Laplace transforms of Well function
#'
#' @param par variable of Laplace tramsform
#' @param Cbar dimensionless welbore storage constant
#' @param Skin skin effect factor
#'
#' @return discrete value of Laplace transforms of Well function
#'
#' @references
#' R. G. Agarwal, R. Al-Hussainy and H. J. Ramey Jr.,
#' “An Investigation of Wellbore Storage and Skin Effect in Unsteady Liquid Flow: I. Analytical Treatment,”
#' SPE Journal, Vol. 10, No. 3, 1970, pp. 279-290.
#'
#' Walton, W.C. 2006. Aquifer Test Modeling. Boca Raton, Florida: CRC Press, Taylor & Francis Group.
AgarwalFunBesselapr = function(par,Cbar,Skin){
  # skin -- skin effect
  denom = par*(sqrt(par)*BesselK1(sqrt(par)) + Cbar * par * (BesselK0(sqrt(par)) + Skin * sqrt(par) * BesselK1(sqrt(par))));
  gwfun = (BesselK0(sqrt(par)) + (Skin * sqrt(par) * BesselK1(sqrt(par)))) / denom;
  return(gwfun);
}

#' #' @export
#' dimHd <- function(x) {
#'   UseMethod("dimHd")
#' }
#' #' @export
#' dimHd.character <- function(x) {paste("char", x)}
#' #' @export
#' dimHd.numeric <- function(x) {paste("num", x)}
#' #' @export
#' dimHd.ups <- function(x) {
#'   paste("list",x$h)
#' }
#' #' @export
#' dimHd.list <- function(x) {
#'   paste("list")
#'   N = length(x$Vi)
#'   # Walton aproximation of Vi
#'   # Vi = c(0.08333333333333333,-32.08333333333334,1279.0,-15623.66666666667,84244.1666666666,-236957.5,375911.66666666667,-340071.6666666667,164062.5,-32812.5)
#'   res=c()
#'   for(j  in 1:length(x$TD)){
#'     timeLT = x$TD[j];
#'     sum = 0;
#'     for(i in 1:N){
#'       k = (i + i)/2;
#'       # sum = sum + Vi[i] * AgarwalFunBESSELapr((i * (log(2)/t)),Cbar,Skin);
#'       sum = sum + x$Vi[i] * AgarwalFunBesselR((i * (log(2)/timeLT)),x$Cbar,x$Skin);
#'     }
#'     res = c(res, log(2) / timeLT * sum);
#'   }
#'   res
#' }
#'
