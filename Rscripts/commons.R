nameit <- 'NPacS1_0.1'
setwd(paste0('~/Roms_tools/Run/',nameit))
source('~/Roms_tools/Rscripts/get_roms_data.R')
#install.packages('plotrix',repos='https://cran.ism.ac.jp/')
library(plot3D)
library(ggplot2)
library(plotrix) #For Taylor diagram

#Coordinates
lon1 = seq(100,280,by=20)
lon2 = lon1
lon2[lon2>180]=lon2[lon2>180]-360

#Read file
nameit  <- 'npacS'
avgfile <- paste0(nameit,'_avg.nc')
biofile <- paste0(nameit,'_dbio_avg.nc')

#Get 'surface area'
pm      <- ncread(avgfile,'pm')  #Unit: m-1
pn      <- ncread(avgfile,'pn')

#Get depth of ROMS
CFF     <- readnc('NO3_roms', avgfile, nameit='npacS',ROMS=T)
Hz      <- CFF$Hz
Depth   <- CFF$depth
Nroms   <- dim(Depth)[3]

#Calculate total volume (m^3) of each computational grid
gvol    <- Hz
for (i in 1:Nroms){
    gvol[,,i] <- Hz[,,i]/(pm*pn)
}

time    <- ncread(avgfile, 'time_step')
time    <- time[1, ] #Time step from initialization (0 year, Jan 1)
dt      <- 1200      #a time step
days    <- as.integer(time*dt/86400)
year    <- days/365

#Get final year
NMo     <- length(days)
NMo     <- (NMo-11):NMo

#DOY of final year
DOYf    <- days[NMo] %% 365

#Get mask:
mask    <- ncread(avgfile,'mask_rho')
NO3_s   <- get_sur('NO3_roms', file = avgfile, nameit = 'npacS')
Lon     <- NO3_s$Lon
Lat     <- NO3_s$Lat
L       <- length(Lon)
M       <- length(Lat)

#Get depth integrated NPP (model)
#ROMS outputs of NPP
#Vertical integration
#Get obs. data:
library(Rcpp)
sourceCpp('~/Roms_tools/Rscripts/integNPP.cpp')
MASK <- function(dat){
   NT <- dim(dat)[3]
   for (i in 1:NT){
       y          <- dat[,,i]
       y[mask==0] <- NA
       dat[,,i]   <- y
   }
   return(dat)
}


vinteg <- function(biofile, var = 'oPPt'){
   NPPm <- ncread(biofile, var, 
                  start = c(1,1,1,NMo[1]), 
                  count = c(L,M,Nroms,length(NMo)))  
   NPPI <- integC(L, M, length(NMo), Nroms, Depth, Hz, NPPm)
   #NPPI <- array(0, dim = c(L, M, 12)) #Integrated NPP to be calculated (unit: mgC m-2 d-1)
   #
   ##Integrate:
   #for (i in 1:L){
   #    for (j in 1:M){
   #         dep <- Depth[i,j,]   #Depth profile
   #         w   <- dep > -260
   #        for (k in 1:12){
   #            npp         <- NPPm[i,j,w,k]  #NPP depth profile
   #            NPPI[i,j,k] <- NPPI[i,j,k] +
   #                           sum(npp*Hz[i,j,w])*mask[i,j] #Unit: mgC m-2 d-1
   #        }
   #    }
   #}
   NPPI      <- array(NPPI, c(L, M, length(NMo)))
   NPPI[,1,] <- NA  #Remove data from the south boundary
   dat       <- MASK(NPPI)
   return(dat)
}

#Calculate the annual mean NPP integrated over the vertical water column
integann <- function(biofile, var = 'oPPt'){
   NPP_               <- vinteg(biofile)
   NPP_ann            <- apply(NPP_, c(1,2), mean)
   NPP_ann[mask == 0] <- NA 
   return(NPP_ann)
}
integT <- function(biofile, var = 'oPPt'){
   dat <- vinteg(biofile, var)
   for (i in 1:dim(dat)[3]){
       dat[,,i] <- dat[,,i]/(pm*pn)
   }
   dat <- apply(dat, c(1,2), mean)
   dat <- sum(dat,na.rm=T)*365   #Unit: mgC y-1
   dat <- dat/1E18    #Unit: PetagC y-1 = 10^15 gC y-1
   return(dat)
}
