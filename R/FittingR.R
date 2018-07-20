#' @export
GwOptim = function(obs=NA,DimlesTime=NA, OPtFunType = "MSE"){

  obsSnorm = obs
  Cbar=0.5
  Skin=-0.5

  switch(OPtFunType,
         MSE={
           fit=DEoptim(Msefit,lower = c(00,0),upper = c(200000,100),DEoptim.control(itermax=100))
         },
         MAE={
           fit=DEoptim(Maefit,lower = c(00,0),upper = c(200000,100),DEoptim.control(itermax=100))
         },
         MAXER={
           fit=DEoptim(MaxErrfit,lower = c(00,0),upper = c(200000,100),DEoptim.control(itermax=100))
         })

  Gwinput$Cbar=fit$optim$bestmem[1]
  Gwinput$Skin=fit$optim$bestmem[2]
  simDDD = DimensionlessDrawDown(Gwinput)
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

