#Get obs. data:
library(Rcpp)
sourceCpp('~/Roms_tools/Rscripts/integNPP.cpp')
clmname <- paste0(nameit, '-clim.nc') 
NO3clm  <- ncread(clmname, 'NO3')
Nz      <- dim(NO3clm)[3]
NO3s    <- NO3clm[,,Nz,]  #Surface NO3 from WOA13
NO3s    <- MASK(NO3s)
NO3a    <- apply(NO3s, c(1,2), mean)  #Seasonal average
NO3a[mask==0] <- NA

#Get seawifs data:
CHLclm  <- ncread(clmname, 'CHLA')
CHLs    <- CHLclm[,,Nz,]  #Surface CHL from SEAWIFS
CHLs    <- MASK(CHLs)

#Get SST from clim file:
source('~/Roms_tools/Rscripts/Picof_Taylor.R')
#library(Rcpp)
#sourceCpp('~/Roms_tools/Rscripts/Ward2015.cpp')
SST       <- ncread(clmname,"temp")
SST       <- SST[,,Nz,]
picof.clm <- Picofrac(CHLs, SST)    #Picophytoplankton fractions from climatology
picof.clm.ann  <- apply(picof.clm, c(1,2), mean)
nanof.clm      <- Nanofrac(CHLs, SST)    #Nanophytoplankton fractions from climatology
microf.clm     <- 1-picof.clm-nanof.clm
nanof.clm.ann  <- apply(nanof.clm, c(1,2), mean)
microf.clm.ann <- apply(microf.clm, c(1,2), mean)

CHL2    <- apply(CHLs, c(1,2),mean)
CHL2[mask==0] <- NA   #Annual mean CHL from SEAWIFS

#Get NPP data from vgpm:
vgpm    <- ncread(clmname,'NPP_vgpm')
cbpm    <- ncread(clmname,'NPP_cbpm')

vgpm.a  <- apply(vgpm, c(1,2), mean)
vgpm.a[mask == 0] <- NA
cbpm.a  <- apply(cbpm, c(1,2), mean)
cbpm.a[mask == 0] <- NA

#Get depth integrated NPP (model)
#ROMS outputs of NPP
#Vertical integration
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
NPPI  <- vinteg(biofile)
NPP.m <- NPPI
vgpm  <- MASK(vgpm)
cbpm  <- MASK(cbpm)
vgpm  <- as.numeric(vgpm)
cbpm  <- as.numeric(cbpm)
npp   <- data.frame( mod =as.numeric(NPP.m), 
                     vgpm=vgpm, 
                     cbpm=cbpm)
npp   <- na.omit(npp)

NPP_ann            <- apply(NPPI, c(1,2), mean)
NPP_ann[mask == 0] <- NA 

NO3m    <- NO3_s$data
NO3m    <- NO3m[,,NMo]
NO3m    <- MASK(NO3m)
CHL     <- ncread(avgfile, 'CHL')
CHL     <- CHL[,,Nroms,NMo]
CHL.m   <- MASK(CHL)

#Plot Taylorgram for NO3:
NO3m    <- as.numeric(NO3m)
NO3s    <- as.numeric(NO3s)
no3     <- data.frame(obs=as.numeric(NO3s), mod=NO3m)
no3[no3 <= 0] <- NA
no3     <- na.omit(no3)

CHLm    <- as.numeric(CHL.m)  #Model
CHLs    <- as.numeric(CHLs)   #Observation
chl     <- data.frame(obs=CHLs, mod=CHLm)
chl[chl <= 0] <- NA
chl     <- na.omit(chl)

#Load mean size data:
LNV <- ncread(avgfile, 'LNV')
LNV <- LNV[,,Nroms,NMo]
LNV <- MASK(LNV)
PHY <- ncread(avgfile, 'PHYTO')
PHY <- PHY[,,Nroms,NMo]
PHY <- MASK(PHY)
LNV <- LNV/PHY
VAR <- ncread(avgfile, 'VAR')
VAR <- VAR[,,Nroms,NMo]
VAR <- MASK(VAR)
VAR <- VAR/PHY - LNV^2

#Calculate the microphytoplankton fraction:
logV            <- function(ESD)log(ESD^3*pi/6)
mean_           <- LNV-log(10)
sd_             <- sqrt(VAR)
picof.m         <- pnorm(logV(2),  mean_, sd_)
nanof.m         <- pnorm(logV(20), mean_, sd_)
microf.m        <- 1 - nanof.m
nanof.m         <- nanof.m - picof.m
picof.m.ann     <- apply(picof.m,c(1,2),function(x)mean(x,na.rm=T))
nanof.m.ann     <- apply(nanof.m, c(1,2), mean)
microf.m.ann    <- apply(microf.m, c(1,2), mean)
picof.m.ann[,1] <- NA
nanof.m.ann[,1] <- NA
microf.m.ann[,1]<- NA

mod_dat_file <- paste0(nameit,'-mod-dat-paper.pdf')

pdf(mod_dat_file,width=6, height=11,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(2.5,3.5,2,2),
             mgp     = c(2.3,1,0),
             oma     = c(4,4,0,0),
             pch     = 16,
             mfrow   = c(5,2))

image2D(NO3m_ann, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,30), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('Modeled NO3 (µM)',adj=0)

image2D(NO3a, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,30), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('WOA13 NO3 (µM)',adj=0)

CHLmax <- 2
CHL1[CHL1 > CHLmax] <- CHLmax
image2D(CHL1, Lon, Lat,    #Model CHL
          col = jet.colors(18),   zlim = c(0,CHLmax), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('Modeled Chl (µg/L)',adj=0)

CHL2[CHL2 > CHLmax] <- CHLmax
image2D(CHL2, Lon, Lat,    #Seawifs CHL
          col = jet.colors(18),   zlim = c(0,CHLmax), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('SeaWIFS Chl (µg/L)',adj=0)

cff <- picof.m.ann
image2D(cff, Lon, Lat, zlim =c(0,1),
          col = jet.colors(18), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- "Modeled pico. fractions"
mtext(Varname,adj=0)

cff <- picof.clm.ann
image2D(cff, Lon, Lat, zlim =c(0,1),
          col = jet.colors(18), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- "Pico. fraction based on Ward (2015)"
mtext(Varname,adj=0,cex=.8)

NPPmax <- 1.6E3
vgpm.a[vgpm.a > NPPmax] <- NPPmax
image2D(vgpm.a, Lon, Lat, zlim =c(0,NPPmax),
          col = jet.colors(18), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- bquote('NPP, vgpm (mg C '*m^-2*' '*d^-1*')')
mtext(Varname,adj=0)

#image2D(eppl.a, Lon, Lat, zlim =c(0,1000),
#          col = jet.colors(18), 
#         xaxt = 'n',frame = F,
#         xlab = "", ylab = "")
#axis(side=1, at = lon1, labels=lon2)
#Varname  <- bquote('NPP, VGPM-Eppley (mg C '*m^-2*' '*d^-1*')')
#mtext(Varname,adj=0)

cbpm.a[cbpm.a > NPPmax] <- NPPmax
image2D(cbpm.a, Lon, Lat, zlim =c(0,NPPmax),
          col = jet.colors(18), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- bquote('NPP, cbpm (mg C '*m^-2*' '*d^-1*')')
mtext(Varname,adj=0)

NPP_ann[NPP_ann > NPPmax] <- NPPmax
NPP_ann[,1]               <- NA
image2D(NPP_ann, Lon, Lat, zlim =c(0,NPPmax),
          col = jet.colors(18), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- bquote('Modeled NPP (mg C '*m^-2*' '*d^-1*')')
mtext(Varname,adj=0)

pico <- data.frame(obs=as.numeric(picof.clm), mod=as.numeric(picof.m))
pico <- na.omit(pico)
nano <- data.frame(obs=as.numeric(nanof.clm), mod=as.numeric(nanof.m))
nano <- na.omit(nano)
micro<- data.frame(obs=as.numeric(microf.clm), mod=as.numeric(microf.m))
micro<- na.omit(micro)
taylor.diagram(log(npp$vgpm), log(npp$mod),
               col=4, add=F, normalize=T, mar=c(2.5,3.5,2,2))
taylor.diagram(log(no3$obs), log(no3$mod), add=T,normalize=T, )
taylor.diagram(log(chl$obs), log(chl$mod), col=3, add=T, normalize=T)
taylor.diagram(log(npp$cbpm),log(npp$mod), col=5, add=T, normalize=T)
taylor.diagram(pico$obs,     pico$mod,     col=6, add=T, normalize=T)
#taylor.diagram(nano$obs,     nano$mod,     col=7, add=T, normalize=T)
#taylor.diagram(micro$obs,    micro$mod,    col=8, add=T, normalize=T)
legend('topright', 
       legend = c('NO3','Chl','NPP,vgpm','NPP,cbpm',
                  '%Pico'),
       col = 2:6, pch = 16, cex=0.6)

mtext('Latitude (ºN)', side=2, outer=T)
mtext('Longitude (ºE)',side=1, outer=T, adj=0.25)
dev.off()

#Another Taylorgram:
TNO3 <- data.frame(MOD=as.numeric(NO3a),OBS=as.numeric(NO3m_ann))
TNO3 <- na.omit(TNO3)
TCHL <- data.frame(MOD=as.numeric(CHL1),OBS=as.numeric(CHL2))
TCHL <- na.omit(TCHL)
TNPP <- data.frame(MOD=as.numeric(NPP_ann),OBS=as.numeric(cbpm.a))
TNPP <- na.omit(TNPP)
TPico<- data.frame(MOD=as.numeric(picof.m.ann),
                   OBS=as.numeric(picof.clm.ann))
TPico<- na.omit(TPico)
Tnano<- data.frame(MOD=as.numeric(nanof.m.ann),
                   OBS=as.numeric(nanof.clm.ann))
Tnano<- na.omit(Tnano)
Tmicro<- data.frame(MOD=as.numeric(microf.m.ann),
                    OBS=as.numeric(microf.clm.ann))
Tmicro<- na.omit(Tmicro)


file1='NewTaylor.pdf'
pdf(file1,width=6, height=6,paper='a4')

op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(2.5,3.5,2,2),
             mgp     = c(2.3,1,0),
             oma     = c(4,4,0,0),
             pch     = 16,
             cex.axis= 1.4,
             cex.lab = 1.4,
             mfrow   = c(1,1))

taylor.diagram(TNPP$OBS, TNPP$MOD, col=2, pch=1,cex=1.5,
               grad.corr.lines=c(0.6,0.9),sd.arcs=c(0.5,1,1.5),
               add=F, normalize=T, mar=c(2.5,3.5,2,2))
taylor.diagram(log(TNO3$OBS), log(TNO3$MOD),pch=1, cex=1.5,
               add=T,normalize=T,col='blue')
taylor.diagram(log(TCHL$OBS), log(TCHL$MOD),pch=1,cex=1.5,
               col='cyan', add=T, normalize=T)
taylor.diagram(TPico$OBS,     TPico$MOD,cex=1.5, col=6, add=T, normalize=T)
taylor.diagram(Tnano$OBS,     Tnano$MOD,cex=1.5, col=7, add=T, normalize=T)
taylor.diagram(micro$obs,    micro$mod,cex=1.5,col=8, add=T, normalize=T)
legend('topright', 
       legend = c('NO3','Chl','NPP',
                  '%Pico','%Micro'),
       col = c('blue','cyan',2,6,8), pch = c(1,1,1,16,16), 
       cex = 0.8)
dev.off()

file1 <- paste0(nameit,'-mod-dat-bulk.pdf')
pdf(file1,width=6, height=7,paper='a4')

op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(2.5,3.5,2,2),
             mgp     = c(2.3,1,0),
             oma     = c(4,4,0,0),
             pch     = 16,
             cex.axis= 1.4,
             cex.lab = 1.4,
             cex     = 1.4,
             mfrow   = c(3,2))

image2D(NO3m_ann, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,30), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('Modeled NO3 (µM)',adj=0, cex=1.2)

image2D(NO3a, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,30), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('WOA13 NO3 (µM)',adj=0,cex=1.2)

CHLmax <- 2
CHL1[CHL1 > CHLmax] <- CHLmax
image2D(CHL1, Lon, Lat,    #Model CHL
          col = jet.colors(18),   zlim = c(0,CHLmax), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('Modeled Chl (µg/L)',adj=0,cex=1.2)

CHL2[CHL2 > CHLmax] <- CHLmax
image2D(CHL2, Lon, Lat,    #Seawifs CHL
          col = jet.colors(18),   zlim = c(0,CHLmax), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('SeaWIFS Chl (µg/L)',adj=0,cex=1.2)

cff <- picof.m.ann
image2D(cff, Lon, Lat, zlim =c(0,1),
          col = jet.colors(18), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- "Modeled pico. fractions"
mtext(Varname,adj=0,cex=1.2)

cff <- picof.clm.ann
image2D(cff, Lon, Lat, zlim =c(0,1),
          col = jet.colors(18), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- "Pico. fraction based on Ward (2015)"
mtext(Varname,adj=0,cex=.8)

mtext('Latitude (ºN)', side=2, outer=T, cex=1.4)
mtext('Longitude (ºE)',side=1, outer=T, adj=0.5,cex=1.4)
dev.off()

#Draw a separate figure with NPP only
file2 <- paste0(nameit,'-mod-dat-NPP.pdf')
pdf(file2,width=8, height=4,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(2.5,3.5,2,3),
             mgp     = c(2.3,1,0),
             oma     = c(4,4,0,0),
             pch     = 16,
             cex     = 1.4,
             cex.axis= 1.4,
             cex.lab = 1.4,
             mfrow   = c(1,2))

NPP_ann[NPP_ann > NPPmax] <- NPPmax
NPP_ann[,1]               <- NA
image2D(NPP_ann, Lon, Lat, zlim =c(0,NPPmax),
          col = jet.colors(18), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- bquote('Modeled NPP (mg C '*m^-2*' '*d^-1*')')
mtext(Varname,adj=0,cex=1.2)

NPPmax <- 1.6E3
#vgpm.a[vgpm.a > NPPmax] <- NPPmax
#image2D(vgpm.a, Lon, Lat, zlim =c(0,NPPmax),
#          col = jet.colors(18), 
#         xaxt = 'n',frame = F,
#         xlab = "", ylab = "")
#axis(side=1, at = lon1, labels=lon2)
#Varname  <- bquote('NPP, vgpm (mg C '*m^-2*' '*d^-1*')')
#mtext(Varname,adj=0,cex=1.2)
cbpm.a[cbpm.a > NPPmax] <- NPPmax
image2D(cbpm.a, Lon, Lat, zlim =c(0,NPPmax),
          col = jet.colors(18), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- bquote('Satellite NPP')
mtext(Varname,adj=0, cex=1.2)

mtext('Latitude (ºN)', side=2, cex=1.4,outer=T)
mtext('Longitude (ºE)',side=1, cex=1.4,outer=T, adj=0.5)
dev.off()

#=====================================
#Calculate NPP from discrete model:
vtr      <- c(0, 0.01, 0.03, 0.05, 0.07, 0.1)
KTW      <- c(0, 0.01, 0.03, 0.05, 0.07, 0.1)
pre      <- '~/Roms_tools/Run/NPacS1_DIS_'
preKTW   <- '~/Roms_tools/Run/NPacS1_DIS_KTW_'
biofile  <- paste0(pre,vtr,'/','npacS_dbio_avg.nc')
avgfile  <- paste0(pre,vtr,'/','npacS_avg.nc')
biof_KTW <- paste0(preKTW,KTW,'/','npacS_dbio_avg.nc')
avgf_KTW <- paste0(preKTW,KTW,'/','npacS_avg.nc')
biof_KTW[1] <- biofile[1]
avgf_KTW[1] <- avgfile[1]

system.time(NPP_VTR  <- sapply(1:length(vtr), function(i)integT(biofile[i])))

system.time(NPP_KTW  <- sapply(1:length(KTW), function(i)integT(biof_KTW[i])))

pdf('TNPP_VTR_DISC.pdf',width=5, height=5,paper = 'a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(4,4,1,1),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             mfrow   = c(1,1)  )
plot(vtr, NPP_VTR, type = 'b',
     xlab = 'Coefficient of TD or KTW',
     ylab = bquote('Annual primary production ('*'PgC '*y^-1*')'))
points(KTW,NPP_KTW, type = 'b', pch=2)
legend('topright', legend = c('TD', 'KTW'), pch=1:2)
dev.off()


#Obtain NO3, CHL, NPP from model files
Mod_sur <- function(avgfile, biofile){
   CHL    <- ncread(avgfile, 'CHL', 
                    start = c(1,1,Nroms,NMo[1]),
                    count = c(L,M,    1,length(NMo)))
   CHL.m  <- MASK(CHL)
   NO3    <- ncread(avgfile, 'NO3', 
                    start = c(1,1,Nroms,NMo[1]),
                    count = c(L,M,    1,length(NMo)))
   NO3.m  <- MASK(NO3)
   NPP.m  <- integann(biofile)
   return(list(CHL=CHL.m, NO3=NO3.m, NPP=NPP.m))
}

dat1  <- Mod_sur(avgfile[4], biofile[4])
file2 <- paste0(nameit,'-mod-dat-bulk-disc.pdf')
pdf(file2,width=6, height=7,paper='a4')

op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(2.5,3.5,2,2),
             mgp     = c(2.3,1,0),
             oma     = c(4,4,0,0),
             pch     = 16,
             cex.axis= 1.4,
             cex.lab = 1.4,
             cex     = 1.4,
             mfrow   = c(3,2))

image2D(NO3m_ann, Lon, Lat, 
          col = jet.colors(18),   zlim = c(0,30), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('Modeled NO3 (µM)',adj=0, cex=1.2)

cff <- apply(dat1$NO3, c(1,2), mean)
image2D(cff, Lon, Lat,  #Modeled NO3 from discrete models
          col = jet.colors(18),   zlim = c(0,30), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('WOA13 NO3 (µM)',adj=0,cex=1.2)

CHLmax <- 2
CHL1   <- apply(dat1$CHL, c(1,2), mean)
CHL1[CHL1 > CHLmax] <- CHLmax
image2D(CHL1, Lon, Lat,    #Model CHL
          col = jet.colors(18),   zlim = c(0,CHLmax), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('Modeled Chl (µg/L)',adj=0,cex=1.2)

CHL2[CHL2 > CHLmax] <- CHLmax
image2D(CHL2, Lon, Lat,    #Seawifs CHL
          col = jet.colors(18),   zlim = c(0,CHLmax), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('SeaWIFS Chl (µg/L)',adj=0,cex=1.2)

#Calculate the pico fraction of each community:
ESD2V   <- function(x) log(pi/6*x^3)
PMU_max <- ESD2V(60)
PMU_min <- ESD2V(.6)
NPHY    <- 20L
PMU_    <- seq(PMU_min, PMU_max, length.out=NPHY)
ip      <- which(PMU_ <= ESD2V(2))

#Read PHYTO from surface:
PHY_res <- function(avgfile){
    Y   <- array(, c(L, M, length(NMo), NPHY))
    for (i in 1:NPHY){
        vname   <- paste0('PHYTO', i)
        Y[,,,i] <- ncread(avgfile, vname,
                        start = c(1,1,Nroms, NMo[1]),
                        count = c(L,M,1,     length(NMo)))
    }
    PHYt <- apply(Y, c(1,2,3), sum)
    pico <- Y[,,,ip]
    pico <- apply(pico,c(1,2,3),sum)/PHYt
    return(list(PHY = PHYt, pico=pico))
}

pico1   <- PHY_res(avgfile[4])$pico
cff     <- apply(pico1, c(1,2), mean)
image2D(cff, Lon, Lat, zlim =c(0,1),
          col = jet.colors(18L), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- "Modeled pico. fractions"
mtext(Varname,adj=0,cex=1.2)

cff <- picof.clm.ann
image2D(cff, Lon, Lat, zlim =c(0,1),
          col = jet.colors(18), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- "Pico. fraction based on Ward (2015)"
mtext(Varname,adj=0,cex=.8)

cbpm.a[cbpm.a > NPPmax] <- NPPmax
image2D(cbpm.a, Lon, Lat, zlim =c(0,NPPmax),
          col = jet.colors(18), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- bquote('NPP, cbpm (mg C '*m^-2*' '*d^-1*')')
mtext(Varname,adj=0)

NPP_ann <- apply(dat1$NPP, c(1,2), mean)
NPP_ann[NPP_ann > NPPmax] <- NPPmax
NPP_ann[,1]               <- NA
image2D(NPP_ann, Lon, Lat, zlim =c(0,NPPmax),
          col = jet.colors(18), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- bquote('Modeled NPP (mg C '*m^-2*' '*d^-1*')')
mtext(Varname,adj=0)

mtext('Latitude (ºN)', side=2, outer=T, cex=1.4)
mtext('Longitude (ºE)',side=1, outer=T, adj=0.5,cex=1.4)
dev.off()

#Compare discrete models with and without flexibility

#Fixed N:C ratio
setwd('~/Roms_tools/Run/NPacS1_DIS_KTW_0.1')
biof1 <- '~/Roms_tools/Run/NPacS1_DIS_KTW_0.1_INFLEX/npacS_dbio_avg.nc'

#Flexible N:C ratio
biof2 <- '~/Roms_tools/Run/NPacS1_DIS_KTW_0.1/npacS_dbio_avg.nc'

NPP1  <- integann(biof1)
NPP2  <- integann(biof2)

NPP_  <- list(NPP1 = NPP1, NPP2 = NPP2)
N     <- length(NPP_)

pdf('Compare_flex.pdf',width=6, height=7,paper='a4')

op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(2.5,3.5,2,3),
             mgp     = c(2.3,1,0),
             oma     = c(4,4,0,0),
             pch     = 16,
             cex.axis= 1.2,
             cex.lab = 1.2,
             cex     = 1.2,
             mfrow   = c(2,1))

   NPPmax <- 1200
   for (i in 1:N){
     cff <- NPP_[[i]]
     cff[cff > NPPmax] <- NPPmax
     cff[,1]           <- NA
     image2D(cff, Lon, Lat, zlim =c(0,NPPmax),
               col = jet.colors(18), 
              xaxt = 'n',frame = F,
              xlab = "", ylab = "")
     axis(side=1, at = lon1, labels=lon2)
     if (i == 1){
       Varname  <- bquote('Modeled NPP without flexibility (mg C '*m^-2*' '*d^-1*')')
     }else{
       Varname  <- bquote('Modeled NPP with flexibility (mg C '*m^-2*' '*d^-1*')')
     }
     mtext(Varname,adj=0)
   }

mtext('Latitude (ºN)', side=2, outer=T, cex=1.4)
mtext('Longitude (ºE)',side=1, outer=T, adj=0.5,cex=1.4)

dev.off()

#Plot vertical profiles of Chl:C
Lon_HOT <- 

