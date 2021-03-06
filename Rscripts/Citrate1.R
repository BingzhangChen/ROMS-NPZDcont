source('~/Roms_tools/Rscripts/commons.R')
# ROMS avg. model data
#Nf  <- system('ls npacific_avg_*.nc | wc -l')
#pat = paste0(nameit,"_avg_*.nc")
#f   = list.files(pattern = pat)
#Calculate total N:

NM      <- length(NMo)  #Number of months
NO3     <- ncread(avgfile, 'NO3', 
                  start = c(1,1,    1,NMo[1]),
                  count = c(L,M,Nroms,length(NMo)))
CHL     <- ncread(avgfile, 'CHL', 
                  start = c(1,1,    1,NMo[1]),
                  count = c(L,M,Nroms,length(NMo)))
PHY     <- ncread(avgfile, 'PHYTO',
                  start = c(1,1,    1,NMo[1]),
                  count = c(L,M,Nroms,length(NMo)))
ZOO     <- ncread(avgfile, 'ZOO')
DET     <- ncread(avgfile, 'DET')
DFe     <- ncread(avgfile, 'DFE')
DETFe   <- ncread(avgfile, 'DETFe')

#Compare model outputs with time series observational data of JGOFS
source('~/Roms_tools/Rscripts/jgofs.R')

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
             oma     = c(2,2,2,2),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             cex.axis= 1.2,
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

