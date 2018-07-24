#' The identification of aquifer parameters using differential evolution
#'
#' @title The identification of aquifer parameters using inverse problem solution
#'
#' @param obs observed dimensionless drwandown
#' @param DimlesTime dimensionless time
#' @param OPtFunType type of inverse problem minimization based on MSE - \code{MSE}, MAE -  \code{MAE}, MAx error  \code{MAXER}
#' @param ... additional parameters of \code{RcppDE}
#'
#' @return list,  well pumping test parameters Cbar and Skin, simulated dimensionles drawdown
#' @export
#'
#'
#' @examples
#'
#' DimlesTime =c (1.88, 3.74, 11.83,47.1,118.31,470.98,1489.37,3741.12, 14893.66)
#' obsSnormDra = c( 0.09899987,0.19699975,0.60699923,2.18499723,4.58399418,9.32698816,11.05998596,11.64,12.38)
#' dd=GwOptim(obsSnormDra,DimlesTime,"MSE",lower=c(0,0),upper=c(200000,200),control = DEoptim.control(itermax=100))
#' plot(log10(DimlesTime), obsSnormDra,pch=19,ylab="Dimensionless Drawdown [-]",xlab = "log10(Dimensionless Time) [-]")
#' lines(log10(DimlesTime), dd$dimDrawdown,col="red")
#'
GwOptim = function(obs=NA,DimlesTime=NA,OPtFunType = "MSE",...){

  obsSnorm <<- obs
  Cbar=0.5
  Skin=-0.5

  # obsSnorm=dat$S

  # Gwinput  = SetStehfestData(dat$T,10,Cbar,Skin,WellTestType= "IsotropicMediumBesselR")

  Gwinput <<- SetStehfestData(DimlesTime,10,Cbar,Skin,WellTestType= "IsotropicMediumBesselR")

  # print(Gwinput)

  # lower = c(00,0),upper = c(200000,100),DEoptim.control(itermax=100)

  switch(OPtFunType,
         MSE={
           fit=DEoptim(Msefit,...)
         },
         MAE={
           fit=DEoptim(Maefit,...)
         },
         MAXER={
           fit=DEoptim(MaxErrfit,...)
         })

  Gwinput$Cbar=fit$optim$bestmem[1]
  Gwinput$Skin=fit$optim$bestmem[2]
  simDDD = DimensionlessDrawDown(Gwinput)

  return(
    list(Cbar = fit$optim$bestmem[1], Skin = fit$optim$bestmem[2], dimDrawdown=simDDD)
  )
}

#' @export
Msefit = function(par){
  Gwinput$Cbar=par[1]
  Gwinput$Skin=par[2]

  simSnorm = DimensionlessDrawDown(Gwinput)
  # N=length(simSnorm)
  # N=20
  if(is.na(simSnorm[1])) simSnorm[1]=0
  res=obsSnorm - simSnorm
  Mse = sum((res)^2)
  return(Mse)
}

#' @export
MaxErrfit = function(par){
  Cbar=par[1]
  Skin=par[2]
  simSnorm= DimensionlessDrawDown(Gwinput)
  # N=length(simSnorm)
  # N=20
  if(is.na(simSnorm[1])) simSnorm[1]=0
  res=obsSnorm - simSnorm
  MaxErr = max(abs(res))
  return(MaxErr)
}

#' @export
Maefit = function(par){
  Cbar=par[1]
  Skin=par[2]
  simSnorm= DimensionlessDrawDown(Gwinput)
  # N=length(simSnorm)
  # N=20
  if(is.na(simSnorm[1])) simSnorm[1]=0
  res=obsSnorm - simSnorm
  Mae = sum(abs(res))
  return(Mae)
}

