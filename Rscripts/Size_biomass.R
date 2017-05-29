#Extract mean size, size variance, and the relationships with total biomass from surface
nameit = 'npacific'
pwd1   = 'npacific'
setwd(paste0('~/Roms_tools/Run/',pwd1))

#Load total biomass:
load('PHY0m.Rdata')

#Load mean size data:
load(paste0('LNVann',pwd1,'.Rdata'))

#Load size variance data:
load(paste0('VARann',pwd1,'.Rdata'))

#Calculate the microphytoplankton fraction:
Nlon = dim(LNV.ann)[1]
Nlat = dim(LNV.ann)[2]

microp <- matrix(NA, nrow=Nlon, ncol=Nlat)
for (i in 1:Nlon){
  for (j in 1:Nlat){
      if (!is.na(LNV.ann[i,j] && !is.na(VAR.ann[i,j]))){
          mean_ <- log(

      }
  }
}
