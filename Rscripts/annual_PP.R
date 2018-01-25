nameit = 'npacS1'
setwd(paste0('~/Roms_tools/Run/',nameit))
source('~/Roms_tools/Rscripts/get_roms_data.R')
library(plot3D)
library(ggplot2)
library(plotrix) #For Taylor diagram

# ROMS avg. model data
#Nf  <- system('ls npacific_avg_*.nc | wc -l')
#pat = paste0(nameit,"_avg_*.nc")
#f   = list.files(pattern = pat)

#Get mask:
nameit  <- 'npacS'
avgfile <- paste0(nameit,'_avg.nc')
NO3_s   <- get_sur('NO3_roms', file = avgfile, nameit = 'npacS')
mask    <- readnc('mask_rho',         avgfile, nameit = nameit)
mask    <- mask$data
Lon     <- NO3_s$Lon
Lat     <- NO3_s$Lat
MASK <- function(dat){
   NT <- dim(dat)[3]
   for (i in 1:NT){
       y          <- dat[,,i]
       y[mask==0] <- NA
       dat[,,i]   <- y
   }
   return(dat)
}


#Get obs. data:
clmname <- paste0(nameit, '-clim.nc') 
NO3clm  <- ncread(clmname, 'NO3')
Nz      <- dim(NO3clm)[3]
NO3s    <- NO3clm[,,Nz,]  #Surface NO3 from WOA13
NO3s    <- MASK(NO3s)

NO3a    <- apply(NO3s, c(1,2), mean)  #Seasonal average
NO3a[mask==0] <- NA

image2D(NO3a, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,30), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")


#Get final year
NMo     <- length(NO3_s$days)
NMo     <- (NMo-11):NMo
NO3m    <- NO3_s$data[,,NMo]   #Final year of NO3 from model
NO3m    <- MASK(NO3m)


#Plot Taylorgram for NO3:
NO3m    <- as.numeric(NO3m)
NO3s    <- as.numeric(NO3s)
no3     <- data.frame(obs=NO3s, mod=NO3m)
no3[no3 <= 0] <- NA
no3     <- na.omit(no3)

#Annual mean NO3:
NO3_s_ann          = apply(NO3, c(1,2),mean)
NO3_s_ann[mask==0] = NA
NO3_s_ann[,1]      = NA

#Get Chl:
CHL   = get_sur('Chl_roms',  file = avgfile, nameit = 'npacS')
CHL   = CHL$data[,,NMo]
CHL.m = MASK(CHL)

#Annual mean CHL:
CHL1          = apply(CHL, c(1,2),mean)
CHL1[mask==0] = NA
CHL1[,1] = NA

#Get seawifs data:
CHLclm  <- ncread(clmname, 'CHLA')
CHLs    <- CHLclm[,,Nz,]  #Surface CHL from WOA13
CHLs    <- MASK(CHLs)
CHLm    <- as.numeric(CHL.m)  #Model
CHLs    <- as.numeric(CHLs)   #Observation
chl     <- data.frame(obs=CHLs, mod=CHLm)
chl[chl <= 0] <- NA
chl     <- na.omit(chl)
taylor.diagram(log(no3$obs), log(no3$mod), normalize=T)
taylor.diagram(log(chl$obs), log(chl$mod), col=3, add=T, normalize=T)

#

#Get mean size:
LNV   = get_sur('LNV',  file = avgfile, nameit = 'npacS')
LNV   = LNV$data[,,NMo]
LNV[LNV <= 0] = NA
#Annual mean LNV:
LNV1   = apply(LNV, c(1,2),mean)
LNV1[mask==0] = NA
LNV1[,1] = NA

#Get size variance:
VAR   = get_sur('VAR',  file = avgfile, nameit = 'npacS')
V1    = VAR$data[wi,wj,]

VAR   = VAR$data[,,NMo]
VAR[VAR <= 0] = NA
#Annual mean VAR:
VAR1   = apply(VAR, c(1,2),mean)
VAR1[mask==0] = NA
VAR1[,1] = NA

#winter
VAR1          = VAR[,,3]
VAR1[mask==0] = NA
VAR1[,1]      = NA
image2D(VAR1, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,10), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")

#Summer:
VAR2          = VAR[,,9]
VAR2[mask==0] = NA
VAR2[,1]      = NA
image2D(VAR2, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,10), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")


ln20   <- log(20**3*pi/6)
ln2    <- log(2 **3*pi/6)
mean_  <- log(LNV**3*pi/6)
sd_    <- sqrt(VAR)
microp <- 1-pnorm(ln20, mean_,sd_)
 picop <-   pnorm(ln2,  mean_,sd_)
 nanop <- 1-microp-picop

#Get data in the last time:
NO3 = NO3_s$data
NO3 = NO3[,,dim(NO3)[3]]
NO3 = NO3[,,1]
NO3[mask==0] = NA
NO3max = quantile(NO3, probs=0.99, na.rm=T)
image2D(NO3, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,NO3max), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")

PHY   = get_sur('PHYTO',  file = avgfile, nameit = 'npacS')
PHY   = PHY$data
PHY1  = PHY[,,dim(PHY)[3]]
PHY1[PHY1 <= 0] = NA
PHYmax = quantile(PHY1, probs=0.99, na.rm=T)
PHY1[PHY1 >= PHYmax] = PHYmax
image2D(PHY1, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,PHYmax), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")

LNV   = get_sur('LNV',    file = avgfile, nameit = 'npacS')
LNV   = LNV$data
LNV1  = LNV[,,dim(LNV)[3]]
LNV2  = LNV[,,1]
LNVmax = quantile(LNV2, probs=0.99, na.rm=T)
LNV1[LNV1 >= LNVmax] = LNVmax
image2D(LNV1, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,15), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")

MTo   = get_sur('MTo',    file = avgfile, nameit = 'npacS')
MTo   = MTo$data
MTo   = MTo/PHY
MTo1  = MTo[,,dim(MTo)[3]]
MTomax = quantile(MTo1, probs=0.99, na.rm=T)
image2D(MTo1, Lon, Lat, 
          col = jet.colors(18),   zlim = c(-5,35), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")

VTo   = get_sur('VTo',    file = avgfile, nameit = 'npacS')
VTo   = VTo$data
VTo   = VTo/PHY - MTo^2
VTo1  = VTo[,,dim(VTo)[3]]
VTomax = quantile(VTo1, probs=0.99, na.rm=T)
image2D(VTo1, Lon, Lat, 
          col = jet.colors(18),   zlim = c(-500,500), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")


MIo   = get_sur('MIo',    file = avgfile, nameit = 'npacS')
MIo   = MIo$data
MIo   = MIo/PHY
MIo1  = MIo[,,dim(MIo)[3]]
MIomax = quantile(MIo1, probs=0.99, na.rm=T)

image2D(MIo1, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,log(1000)), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")

VIo   = get_sur('VIo',    file = avgfile, nameit = 'npacS')
VIo   = VIo$data
VIo   = VIo/PHY - MIo^2
VIo1  = VIo[,,dim(VIo)[3]]
VIomax = quantile(VIo1, probs=0.99, na.rm=T)
image2D(VIo1, Lon, Lat, 
          col = jet.colors(18),   zlim = c(-5,10), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")


mask    <- mask$data

VAR        = 'Chl_roms'

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

NO30m = Gets.list('NO3_roms')
save(NO30m, file = 'NO30m.Rdata')

Chl0m = Gets.list('Chl_roms')
save(Chl0m, file = 'Chl0m.Rdata')

LNV0m = Gets.list('LNV_roms')
save(LNV0m, file = 'LNV0m.Rdata')

avgfile2  <- '~/Roms_tools/run/NpacS_0.001/npacS_avg.nc'
avgfile3  <- '~/Roms_tools/run/NpacS_NODIV/npacS_avg2.nc'

#Check file without temp. dependency:
avgfile4  <- '~/Roms_tools/run/NpacS_NOTEMP/npacS_avg11.nc' 

# ROMS biological diagnostic file:
biofile   <- '~/Roms_tools/run/NpacS/npacS_dbio_avg.nc'
mu   = get_sur('omuNet',    file = biofile, nameit = 'npacS')
mu   = mu$data
mu1  = mu[,,dim(mu)[3]]
mumax = quantile(mu1, probs=0.99, na.rm=T)

image2D(mu1, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,4), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")

dmudl   = get_sur('odmuNetdl',    file = biofile, nameit = 'npacS')
dmudl   = dmudl$data
dmudl   = dmudl[,,1]
dmudl[mask==0] = NA
image2D(dmudl, Lon, Lat, 
          col = jet.colors(18),
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")

dmudZ   = get_sur('odmudZ',    file = biofile, nameit = 'npacS')
dmudZ   = dmudZ$data
dmudZ   = dmudZ[,,dim(dmudZ)[3]]
dmudZ[mask==0] = NA
image2D(dmudZ, Lon, Lat, 
          col = jet.colors(18),
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")

d2mudZ2   = get_sur('od2mudZ2',    file = biofile, nameit = 'npacS')
d2mudZ2   = d2mudZ2$data
d2mudZ2   = d2mudZ2[,,dim(d2mudZ2)[3]]
d2mudZ2[mask==0] = NA
image2D(d2mudZ2, Lon, Lat, 
          col = jet.colors(18),
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")


biofile2  <- '~/Roms_tools/run/NpacS_0.001/npacS_dbio_avg.nc'
biofile3  <- '~/Roms_tools/run/NpacS_NODIV/npacS_dbio_avg2.nc'
biofile4  <- '~/Roms_tools/run/NpacS_NOTEMP/npacS_dbio_avg11.nc'


#Get NO3:
NO3_s  = get_sur('NO3_roms', file = avgfile, nameit = 'npacific')

NO3_s_ann  = get_ann_sur('NO3_roms', file = avgfile, nameit = 'npacS')
NO3_s_ann2 = get_ann_sur('NO3_roms', file = avgfile3)
NO3_s_ann4 = get_ann_sur('NO3_roms', file = avgfile4)

#Get temp:
temp_s_ann  = get_ann_sur('temp_roms', file = avgfile)

#Get PAR:
PAR_s_ann  = get_ann_sur('radsw', file = avgfile)
PAR_s_ann$data = PAR_s_ann$data * .43 * mask$data
PAR_s_ann$data[PAR_s_ann$data <= 0] <- NA

#Get mixed layer depth:
MLD_s_ann  = get_ann_sur('hbl', file = avgfile)

#Get Chl a:
Chl_s_ann  = get_ann_sur('Chl_roms', file = avgfile)

#Get Zooplankton:
ZOO_s_ann  = get_ann_sur('ZOO', file = avgfile)
ZOO_s_ann2 = get_ann_sur('ZOO', file = avgfile3)

#Get PHYTO biomass:
PHYTO_s_ann = get_ann_sur('PHYTO', file = avgfile)
PHYTO_s_ann2= get_ann_sur('PHYTO', file = avgfile3)

#Get fer:
DFE_s_ann  = get_ann_sur('DFE', file = avgfile)
DFE_s_ann3 = get_ann_sur('DFE', file = avgfile3)

#Get NPP:
NPP_s_ann  = get_ann_sur('oPPt', file = biofile)
NPP_s_ann1 = get_ann_sur('oPPt', file = biofile2)
NPP_s_ann2 = get_ann_sur('oPPt', file = biofile3)
NPP_s_ann4 = get_ann_sur('oPPt', file = biofile4)

#Get mean size:
LNV_s_ann  = get_ann_sur('LNV',  file = avgfile)
LNV_s_ann1 = get_ann_sur('LNV',  file = avgfile2)
#LNV_s_ann4 = get_ann_sur('LNV',  file = avgfile4)

#Get size variance 
VAR_s_ann  = get_ann_sur('VAR',  file = avgfile)
VAR_s_ann1 = get_ann_sur('VAR',  file = avgfile2)
#VAR_s_ann4 = get_ann_sur('VAR',  file = avgfile4)

#Get growth rate at the mean size:
muNet  = get_sur('omuNet', file = biofile)
muNet1 = get_sur('omuNet', file = biofile2)
#muNet2 = get_sur('omuNet', file = biofile3)  #Growth rate from fixed size
#muNet4 = get_sur('omuNet', file = biofile4)  #Growth rate from NO temperature dependency

#Get d2mudl2:
d2mudl2  = get_sur('od2mudl2', file = biofile)
d2mudl21 = get_sur('od2mudl2', file = biofile2)
#d2mudl24 = get_sur('od2mudl2', file = biofile4)

#Get size variance:
VAR_s  = get_sur('VAR',  file = avgfile)
VAR_s1 = get_sur('VAR',  file = avgfile2)
#VAR_s4 = get_sur('VAR',  file = avgfile4)

#Calculate avg. growth rate:
muAvg  =  muNet$data + .5*VAR_s$data*d2mudl2$data
muAvg1 = muNet1$data + .5*VAR_s1$data*d2mudl21$data
#muAvg4 = muNet4$data + .5*VAR_s4$data*d2mudl24$data

#Calculate CV of muAvg:

CV  <- function(x)sd(x,na.rm=T)/mean(x,na.rm=T)
muAvg.CV  = apply(muAvg,  c(1,2), CV)
muAvg1.CV = apply(muAvg1, c(1,2), CV)

#Calculate mean of muAvg:
 muAvg.ann  = apply(muAvg,  c(1,2), mean)
 muAvg.ann[muAvg.ann <= 0] = NA
muAvg1.ann  = apply(muAvg1, c(1,2), mean)
 muAvg1.ann[muAvg1.ann <= 0] = NA
#muAvg2.ann  = apply(muNet2$data, c(1,2), mean)
# muAvg2.ann[muAvg2.ann <= 0] = NA
#muAvg4.ann  = apply(muNet4$data, c(1,2), mean)
# muAvg4.ann[muAvg4.ann <= 0] = NA



#Plot environmental conditions:
source('~/Roms_tools/Rscripts/Annu_mean_Env.R')

#Plot NPP and productivity for no temperature dependency:
source('~/Roms_tools/Rscripts/NO_TEMP.R')



#Plot annual patterns of three different treatments:
source('~/Roms_tools/Rscripts/ann_NPP_VAR.R')

#Plot productivity differences:
source('~/Roms_tools/Rscripts/Plot_diff.R')

#Get time series data for the two boxes:
source('~/Roms_tools/Rscripts/get_BOX.R')
#Plot out:
source('~/Roms_tools/Rscripts/TS_2box.R')

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
  plot(as.vector(NO3_v), as.vector(muAvg.dff),pch=16, cex=.6, ylim=Z_lim,
       xlab='Variability of DIN',ylab='Effect of diversity on productivity')
  abline(0,0, lty=2)
  ds    = dat_lim(as.vector(muAvg.dff), Z_lim)
  yhist = hist(ds, plot=FALSE, breaks=25)
  par(mar=c(4,0,1,.1))
  barplot(yhist$counts, axes=FALSE, space=0, horiz=TRUE)

  #Plot productivity difference with productivity:
#  plot(as.vector(log(NPP_s_ann$data)), as.vector(muAvg.dff),pch=16, cex=.6, ylim=Z_lim,
#       xlab='Annual mean NPP',ylab='Effect of diversity on productivity')
dev.off()


#Test the effect of diversity on stability:


#Plot density plot between size diversity and NPP:
source('~/Working/FlexEFT1D/Rscripts/LO_theme.R')
Y   <- data.frame(VAR=as.vector(VAR_s_ann$data), 
                  NPP=as.vector(muAvg.ann))
Y   <- na.omit(Y)
Y   <- Y[Y$VAR > 0,]

Y1  <- data.frame(VAR=as.vector(VAR_s_ann1$data), 
                  NPP=as.vector(muAvg1.ann))
Y1  <- na.omit(Y1)
Y1  <- Y1[Y1$VAR > 0,]

op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(4,4,1,1),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             cex     = .6,
             pch     = 16,
             mfrow   = c(1,1))

plot(Y$NPP, Y$VAR, pch=16, xlim=c(0,1.5),ylim=c(0,10),
     xlab='Productivity',ylab='Size diversity')
points(Y1$NPP, Y1$VAR, pch=16, col=2)
#At annual mean scales: 
xx  <- ggplot(Y, aes(VAR, NPP)) +
          stat_density2d(aes(fill = ..level..), geom="polygon")+
          scale_fill_gradientn(colours=jet.colors(20))+  
          xlab('Size diversity') +
          ylab('Productivity') +
          LO_theme 

