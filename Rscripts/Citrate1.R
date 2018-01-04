nameit <- 'NPacS1_0.1'
setwd(paste0('~/Roms_tools/Run/',nameit))
source('~/Roms_tools/Rscripts/get_roms_data.R')
#install.packages('devtools',repos='https://cran.ism.ac.jp/')
library(plot3D)
library(ggplot2)
library(plotrix) #For Taylor diagram

#Read file
nameit  <- 'npacS'
avgfile <- paste0(nameit,'_avg1.nc')
biofile <- paste0(nameit,'_dbio_avg1.nc')

#Get 'surface area'
pm      <- ncread(avgfile,'pm')  #Unit: m-1
pn      <- ncread(avgfile,'pn')

#Get depth of ROMS
CFF     <- readnc('NO3_roms', avgfile, nameit='npacS',ROMS=T)
Hz      <- CFF$Hz
Depth   <- CFF$depth
Nroms   <- dim(Depth)[3]
time    <- ncread(avgfile, 'time_step')
time    <- time[1,]  #Time step from initialization (0 year, Jan 1)
dt      <- 1200      #a time step
days    <- as.integer(time*dt/86400)
year    <- days/365
# ROMS avg. model data
#Nf  <- system('ls npacific_avg_*.nc | wc -l')
#pat = paste0(nameit,"_avg_*.nc")
#f   = list.files(pattern = pat)
#Calculate total N:

NO3     <- ncread(avgfile, 'NO3')
PHY     <- ncread(avgfile, 'PHYTO')
ZOO     <- ncread(avgfile, 'ZOO')
DET     <- ncread(avgfile, 'DET')
DFe     <- ncread(avgfile, 'DFE')
DETFe   <- ncread(avgfile, 'DETFe')

#Get mask:
mask    <- ncread(avgfile,'mask_rho')
NO3_s   <- get_sur('NO3_roms', file = avgfile, nameit = 'npacS')
Lon     <- NO3_s$Lon
Lat     <- NO3_s$Lat
L       <- length(Lon)
M       <- length(Lat)
source('~/Roms_tools/Rscripts/inventory.R')
#Retrieve annual surface mean data:
Sur_mean <- function(ncfile, VAR){
   DAT  <- ncread(ncfile, VAR)
   DAT  <- DAT[,,Nroms,NMo]
   
   #Annual mean DAT:
   DAT.a <- apply(DAT, c(1,2),mean)
   DAT.a[mask==0] <- NA
   DAT.a[,1]      <- NA
   return(DAT.a)
}
#Get surface data of modeled data:
MASK <- function(dat){
   NT <- dim(dat)[3]
   for (i in 1:NT){
       y          <- dat[,,i]
       y[mask==0] <- NA
       dat[,,i]   <- y
   }
   return(dat)
}

#Get final year
NMo     <- length(days)
NMo     <- (NMo-11):NMo

#Get surface DFe:
DFe     <- DFe[,,Nroms,NMo]

#Annual mean surface DFe
DFe.a   <- apply(DFe, c(1,2),mean)
DFe.a[mask == 0] <- NA

#Annual mean NO3:
NO3m_ann <- Sur_mean(avgfile, 'NO3')
#Annual mean CHL:
CHL1 <- Sur_mean(avgfile, 'CHL')

#longitude for plotting
lon1 = seq(100,280,by=20)
lon2 = lon1
lon2[lon2>180]=lon2[lon2>180]-360

source('~/Roms_tools/Rscripts/mod_dat_paper.R')
#source('~/Roms_tools/Rscripts/mod_dat_month.R')
source('~/Roms_tools/Rscripts/plot_bio.R')
source('~/Roms_tools/Rscripts/Size_Biomass.R')
source('~/Roms_tools/Rscripts/NPP_VAR.R')
source('~/Roms_tools/Rscripts/Plot_diff.R')

#Get time series data for the two boxes:
source('~/Roms_tools/Rscripts/get_BOX.R')
#Plot out:
source('~/Roms_tools/Rscripts/TS_2box.R')
#Make a video:
#CHL.m 
for (i in 1:dim(CHL.m)[3]){
   tmp <- CHL.m[,,i]
   tmp[tmp > CHLmax] <- CHLmax
   tmp[,1]           <- NA
   fname             <- sprintf("%02d", i)
   fname             <- paste0('CHL',fname,'.png')
   png(fname)
   image2D(tmp, Lon, Lat,    #Model CHL
             col = jet.colors(18),   zlim = c(0,CHLmax), 
            xaxt = 'n',frame = F,
            xlab = "Longitude (ºE)", ylab = "Latitude (ºN)")
   axis(side=1, at = lon1, labels=lon2)
   dev.off()
}

system('ffmpeg -framerate 4 -i CHL%03d.png output.mp4')
#Plot out physics variables:
#Annual mean temperature:
Temp <- Sur_mean(avgfile, 'temp')
image2D(Temp, Lon, Lat, 
          col = jet.colors(18),
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('Modeled SST (ºC)',adj=0)

#Get mixed layer depth:
MLD  <- ncread(avgfile,'hbl')
MLD  <- MLD[,,NMo]#Last year

#Average MLD:
MLD.mean  <- apply(MLD, c(1,2), mean)
image2D(MLD.mean, Lon, Lat, 
          col = jet.colors(18),
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('Annual mean MLD (m)',adj=0)

#Variation of MLD (CV):
MLD.SD  <- apply(MLD, c(1,2), sd)
MLD.CV  <- MLD.SD/MLD.mean
image2D(MLD.CV, Lon, Lat, 
          col = jet.colors(18),
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('CV of MLD (m)',adj=0)



#Test the hypothesis whether environmental variability affects the effect of diversity on productivity
#Plot productivity difference vs. variability of NO3:
  #Calculate variability of NO3:
  NO3   <- get_sur('NO3_roms',  file = avgfile)
  NO3_v <- apply(NO3$data, c(1,2), CV)

pdf('NO3_var_mu_diff.pdf',width=5, height=4,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(4,4,1,1),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             pch     = 16,
             mfrow   = c(1,2))
  #NPP.dff1 <- log(muAvg.ann/muAvg1.ann)
  Z_lim    <- c(-.2,.2)
  zones=matrix(c(1,2), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5))
  par(mar=c(4,4,1,.1))
  plot(as.vector(log(NO3m_ann)), as.vector(NPP.dff),pch=16, cex=.6, ylim=Z_lim,
       xlab='Variability of MLD',ylab='Effect of diversity on productivity')
  abline(0,0, lty=2)
  ds    = dat_lim(as.vector(muAvg.dff), Z_lim)
  yhist = hist(ds, plot=FALSE, breaks=25)
  par(mar=c(4,0,1,.1))
  barplot(yhist$counts, axes=FALSE, space=0, horiz=TRUE)

  #Plot productivity difference with productivity:
#  plot(as.vector(log(NPP_s_ann$data)), as.vector(muAvg.dff),pch=16, cex=.6, ylim=Z_lim,
#       xlab='Annual mean NPP',ylab='Effect of diversity on productivity')
dev.off()

