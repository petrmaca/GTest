#' Setting the information about the evaluated pumping test into basic list. Current  Well test seted for stehfest algorithm on Laplace tranform for isotropic medium aquifers with radial flow.
#'
#' @title Pupming Test Data
#'
#' @param Q pumping test rate  [m3/s]
#' @param Storativity aquifer storativity [-]
#' @param Transmisivity aquifer tranmisivity [m2/s]
#' @param Radius well radius [m]
#' @param Time numeric vector time intervals of puming test [s]
#' @param DrawDown numeric vector drawdown recorded for each time of pupmping test [m]
#' @param WellTestType the identifier of type of well test , \code{"IsotropicMediumBesselR"} stands for transient radial flow and skin effect to isotropic aquifer with R bessels functions
#'
#' @return list of pumpiming test data
#' @export
#'
#' @examples
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
#' Kytlice_2018_15_05 = PumpingTestDataLoad(Q=Q,Storativity = STRR,
#' Transmisivity = Transm,Radius = Rv,Time = T,DrawDown = S,WellTestType = "IsotropicMediumBesselR")
PumpingTestDataLoad = function(Q=NA,Storativity=NA,Transmisivity=NA,Radius=NA,Time=NA,DrawDown=NA,WellTestType = "IsotropicMediumBesselR"){

  if(!is.numeric(Q) | (Q<0)) stop("Pumpimg rate Q is not numeric or is outside the proper bounds.")
  if(!is.numeric(Storativity) | (Storativity<0)) stop("Pumpimg rate Q is not numeric or is outside the proper bounds.")
  if(!is.numeric(Transmisivity) | (Transmisivity<0)) stop("Pumpimg rate Q is not numeric or is outside the proper bounds.")
  if(!is.numeric(Radius) | (Radius<0)) stop("Pumpimg rate Q is not numeric or is outside the proper bounds.")
  if(anyNA(Time)) stop("The time had NA values.")
  if(anyNA(DrawDown)) stop("The time had NA values.")
  if(!is.character(WellTestType)) stop("The WellTestType should be a character.")


  DimTime = Time*Transmisivity/Radius/Radius/Storativity
  # # Dimensionless drawdown [-]
  # dlS = (S*2*3.14*Transm)/Q
  DimDrawDown = (DrawDown*2*3.14*Transmisivity)/Q

  res = list(Q=Q,
       Storativity=Storativity,
       Transmisivity=Transmisivity,
       Radius=Radius,
       Time=Time,
       DrawDown=DrawDown,
       DimlesTime = DimTime,
       DimDrawDown = DimDrawDown,
       WellTestType = WellTestType)

  class(res) <-WellTestType

  return(res)
}


#' Object for calculating dimensionles drawdaown
#'
#' @title Object for calculating dimensionles drawdaown
#'
#' @param PumpingTestDataList data
#' @param N number
#' @param Cbar c bar
#' @param Skin skin factor
#' @param WellTestType aquifer type
#'
#' @return object for calculating Dimensionless drawdown using stehfest algorthm of aquifer \code{IsotropicMediumBesselR} type
#' @export
#'
#' @examples
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
#' Kytlice_2018_15_05 = PumpingTestDataLoad(Q=Q,Storativity = STRR,
#' Transmisivity = Transm,Radius = Rv,Time = T,DrawDown = S,WellTestType = "IsotropicMediumBesselR")
#' Gwinput = SetStehfestData(Kytlice_2018_15_05$DimlesTime,
#'                          N=10,
#'                          Cbar= 1562.858133,
#'                          Skin = 4.976677,
#'                          WellTestType= "IsotropicMediumBesselR")
SetStehfestData <- function(Time=NA, N=10, Cbar = NA, Skin =NA, WellTestType = "IsotropicMediumBesselR"){
  if(is.na(N) | (N!=10)) stop("The number of Vi in stehfest aproximation should be used.")

  Vi = c(0.08333333333333333,-32.08333333333334,1279.0,-15623.66666666667,
         84244.1666666666,-236957.5,375911.66666666667,
         -340071.6666666667,164062.5,-32812.5)

  InpList = list(Vi = Vi,
                 TD=Time,
                 Cbar = Cbar,
                 Skin = Skin)

  class(InpList) <-as.character(WellTestType)
  return(InpList)
}
