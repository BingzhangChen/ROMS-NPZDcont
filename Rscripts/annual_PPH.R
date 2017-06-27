pwd1   = 'npacific_0.001'
setwd(paste0('~/Roms_tools/Run/',pwd1))
source('~/Roms_tools/Rscripts/get_roms_data.R')
library(plot3D)
library(ggplot2)

# ROMS avg. model data
#Nf        <- system('ls npacific_avg_*.nc | wc -l')
#pat = paste0(nameit,"_avg_*.nc")
#f   = list.files(pattern = pat)

#Get mask:
nameit  <- 'npacific'
avgfile <- paste0(nameit,'_avg_1.nc')
NO3_s   <- get_sur('NO3_roms', file = avgfile, nameit = 'npacific')
mask    <- readnc('mask_rho', avgfile, nameit = nameit)
Lon     <- NO3_s$Lon
Lat     <- NO3_s$Lat
latlon  <- list(Lon=Lon, Lat=Lat)
save(latlon, file = 'latlon.Rdata')
mask    <- mask$data

Gets.list = function(VAR, nameit = 'npacific', BIOfile = F){
   N_per_y    = 365
   N_per_file = 5
   NFile      = N_per_y/N_per_file
   Days       = numeric(N_per_y)
   Dat        = array(NA, dim = c(length(Lon), length(Lat), length(Days)))
   
   for (i in 1:NFile){
     if (BIOfile){
      avgfile <- paste0(nameit,'_dbio_avg_',i,'.nc') 
     }else {
      avgfile <- paste0(nameit,'_avg_',i,'.nc')
     }
     dat      <- get_sur(VAR, file = avgfile, nameit = nameit)
     w        <- ((i-1)*N_per_file+1):(i*N_per_file)
     Days[w]  <- dat$days
     Dat[,,w] <- dat$data
   }
   return(list(Days=Days,Dat=Dat))
}

PHY0m = Gets.list('PHYTO')
save(PHY0m, file = 'PHY0m.Rdata')

load('PHY0m.Rdata')

Days    = PHY0m$Days
save(Days, file = 'Days.Rdata')

PHY.ann = apply(PHY0m$Dat,c(1,2),mean)
PHY.ann[PHY.ann <= 0] <- NA
save(PHY.ann, file = 'PHYann.Rdata')

NO30m = Gets.list('NO3_roms')
save(NO30m, file = 'NO30m.Rdata')

load('NO30m.Rdata')
NO3.ann = apply(NO30m$Dat,c(1,2),mean)
NO3.ann[NO3.ann <= 0] <- NA
save(NO3.ann, file = paste0('NO3ann.Rdata'))

Chl0m = Gets.list('Chl_roms')
save(Chl0m, file = 'Chl0m.Rdata')
Chl.ann = apply(Chl0m$Dat,c(1,2),mean)
Chl.ann[Chl.ann <= 0] <- NA
save(Chl.ann, file = paste0('Chlann.Rdata'))

LNV0m = Gets.list('LNV')
save(LNV0m, file = 'LNV0m.Rdata')

load('LNV0m.Rdata')
LNV.ann = apply(LNV0m$Dat,c(1,2),mean)
LNV.ann[LNV.ann <= 0] <- NA
save(LNV.ann, file = 'LNVann.Rdata')

VAR0m = Gets.list('VAR')
save(VAR0m, file = 'VAR0m.Rdata')
load('VAR0m.Rdata')
VAR.ann = apply(VAR0m$Dat,c(1,2),mean)
VAR.ann[VAR.ann <= 0] <- NA
save(VAR.ann, file = 'VARann.Rdata')

NPP0m = Gets.list('oPPt', BIOfile = T)
save(NPP0m, file = 'NPP0m.Rdata')

load('NPP0m.Rdata')
NPP.ann               <- apply(NPP0m$Dat,c(1,2),mean)
NPP.ann[NPP.ann <= 0] <- NA
save(NPP.ann, file = 'NPPann.Rdata')

mu0m = Gets.list('omuNet', BIOfile = T)
save(mu0m, file = 'mu0m.Rdata')

d2mu0m = Gets.list('od2mudl2', BIOfile = T)
save(d2mu0m, file = 'd2mu0m.Rdata')

load('mu0m.Rdata')
load('d2mu0m.Rdata')
muAvg  = mu0m$Dat + .5*VAR0m$Dat*d2mu0m$Dat
save(muAvg, file = 'muAvg.Rdata')

load('muAvg.Rdata')
muAvg.ann <- apply(muAvg, c(1,2), mean)
save(muAvg.ann, file = 'muAvg.ann.Rdata')

source('~/Roms_tools/Rscripts/get_BOXH.R')
source('~/Roms_tools/Rscripts/TS_2box.R')

#Test the hypothesis whether environmental variability affects the effect of diversity on productivity
#Plot productivity difference vs. variability of NO3:

#Calculate variability of NO3:
load('~/Roms_tools/Run/npacific/NO30m.Rdata')
#NO3   <- get_sur('NO3_roms',  file = avgfile)
CV    <- function(x)sd(x,na.rm=T)/mean(x,na.rm=T)
NO3_v <- apply(NO30m$Dat, c(1,2), CV)

#Calculate kernel densities of NO3_v vs. muAvg.dff:
f2   <- data.frame(x=as.vector(NO3_v), y=as.vector(muAvg.dff))
f2   <- na.omit(f2)
#f2   <- kde2d(f2$x, f2$y, n = 50, lims = c(0, 1.5, -0.18, .18))
cor.test(f2$x, f2$y)
N    <- 5000
f3   <- sample(1:nrow(f2), N) 
f3   <- f2[f3,]
pdf('NO3_var_mu_diff_H.pdf',width=5, height=4,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(4,4,1,1),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             pch     = 16,
             mfrow   = c(1,2))
  #NPP.dff1 <- log(muAvg.ann/muAvg1.ann)
  Z_lim    <- c(-.18,.18)
  zones=matrix(c(1,2), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5))
  par(mar=c(4,4,1,.1))

  #image2D(f2, col = jet.colors(20), colkey=F, # zlim = c(0,.05), 
  #     xlab='Variability of DIN',ylab='Effect of diversity on productivity')
  plot(f3$x, f3$y,pch=16, cex=.6, ylim=Z_lim,
       xlab='Variability of DIN',ylab='Effect of diversity on productivity')

  abline(0,0, lty=2)
  abline(lm(y ~ x, f2))
  ds    = dat_lim(as.vector(muAvg.dff), Z_lim)
  yhist = hist(ds, plot=FALSE, breaks=25)
  par(mar=c(4,0,1,.1))
  barplot(yhist$counts, axes=FALSE, space=0, horiz=TRUE)

  #Plot productivity difference with productivity:
#  plot(as.vector(log(NPP_s_ann$data)), as.vector(muAvg.dff),pch=16, cex=.6, ylim=Z_lim,
#       xlab='Annual mean NPP',ylab='Effect of diversity on productivity')
dev.off()

