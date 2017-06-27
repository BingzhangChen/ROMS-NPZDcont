nameit = 'npacific'
setwd(paste0('~/Roms_tools/Run/',nameit))
source('~/Roms_tools/Rscripts/get_roms_data.R')
library(plot3D)
library(ggplot2)

# ROMS avg. model data
#Nf        <- system('ls npacific_avg_*.nc | wc -l')
#pat = paste0(nameit,"_avg_*.nc")
#f   = list.files(pattern = pat)

#Get mask:
avgfile <- paste0(nameit,'_avg_1.nc')
NO3_s   <- get_sur('NO3_roms', file = avgfile, nameit = 'npacific')
mask    <- readnc('mask_rho', avgfile, nameit = nameit)
Lon     <- NO3_s$Lon
Lat     <- NO3_s$Lat

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
biofile   <- '~/Roms_tools/run/NpacS/npacS_dbio_avg11.nc'
biofile2  <- '~/Roms_tools/run/NpacS_0.001/npacS_dbio_avg.nc'
biofile3  <- '~/Roms_tools/run/NpacS_NODIV/npacS_dbio_avg2.nc'
biofile4  <- '~/Roms_tools/run/NpacS_NOTEMP/npacS_dbio_avg11.nc'


#Get NO3:
NO3_s  = get_sur('NO3_roms', file = avgfile, nameit = 'npacific')

NO3_s_ann  = get_ann_sur('NO3_roms', file = avgfile, nameit = 'npacific')
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

