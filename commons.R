nameit <- 'NPacS1_0.1'
setwd(paste0('~/Roms_tools/Run/',nameit))
source('~/Roms_tools/Rscripts/get_roms_data.R')
#install.packages('plotrix',repos='https://cran.ism.ac.jp/')
library(plot3D)
library(ggplot2)
library(plotrix) #For Taylor diagram

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


