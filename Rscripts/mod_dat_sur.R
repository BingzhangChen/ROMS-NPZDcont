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

#Get obs. data:
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
picof.clm.ann <- apply(picof.clm, c(1,2), mean)
nanof.clm <- Nanofrac(CHLs, SST)    #Picophytoplankton fractions from climatology

CHL2    <- apply(CHLs, c(1,2),mean)
CHL2[mask==0] <- NA   #Annual mean CHL from SEAWIFS

#Get NPP data from vgpm:
vgpm    <- ncread(clmname,'NPP_vgpm')
eppl    <- ncread(clmname,'NPP_eppl')
cbpm    <- ncread(clmname,'NPP_cbpm')

vgpm.a  <- apply(vgpm, c(1,2), mean)
vgpm.a[mask == 0] <- NA
eppl.a  <- apply(eppl, c(1,2), mean)
eppl.a[mask == 0] <- NA
cbpm.a  <- apply(cbpm, c(1,2), mean)
cbpm.a[mask == 0] <- NA

#Get depth integrated NPP (model)
#ROMS outputs of NPP
integ <- function(biofile, var = 'oPPt'){
   NPPm <- ncread(biofile, var)  
   NPPm <- NPPm[,,,NMo]
   NPPI <- array(0, dim = c(L, M, 12)) #Integrated NPP to be calculated (unit: mgC m-2 d-1)
   
   #Integrate:
   for (i in 1:L){
       for (j in 1:M){
            dep <- Depth[i,j,]   #Depth profile
            w   <- dep > -260
           for (k in 1:12){
               npp         <- NPPm[i,j,w,k]  #NPP depth profile
               NPPI[i,j,k] <- NPPI[i,j,k] +
                              sum(npp*Hz[i,j,w])*mask[i,j] #Unit: mgC m-2 d-1
           }
       }
   }
   
   NPPI[,1,] <- NA  #Remove data from the south boundary
   dat       <- MASK(NPPI)
   return(dat)
}

integann <- function(biofile, var = 'oPPt'){
   NPP_               <- integ(biofile)
   NPP_ann            <- apply(NPP_, c(1,2), mean)
   NPP_ann[mask == 0] <- NA 
   return(NPP_ann)
}
integT <- function(biofile, var = 'oPPt'){
   dat <- integ(biofile, var)
   for (i in 1:dim(dat)[3]){
       dat[,,i] <- dat[,,i]/(pm*pn)
   }
   dat <- apply(dat, c(1,2), mean)
   dat <- sum(dat,na.rm=T)*365   #Unit: mgC y-1
   dat <- dat/1E18    #Unit: PetagC y-1 = 10^15 gC y-1
   return(dat)
}
NPPI  <- integ(biofile)
NPP.m <- NPPI
vgpm  <- MASK(vgpm)
cbpm  <- MASK(cbpm)
eppl  <- MASK(eppl)
vgpm  <- as.numeric(vgpm)
cbpm  <- as.numeric(cbpm)
eppl  <- as.numeric(eppl)
npp   <- data.frame( mod =as.numeric(NPP.m), 
                     vgpm=vgpm, 
                     cbpm=cbpm)
 #                    eppl=eppl )
npp   <- na.omit(npp)

NPP_ann            <- apply(NPPI, c(1,2), mean)
NPP_ann[mask == 0] <- NA 

NO3m <- NO3_s$data
NO3m <- NO3m[,,NMo]
NO3m <- MASK(NO3m)
CHL  <- ncread(avgfile, 'CHL')
CHL  <- CHL[,,Nroms,NMo]
CHL.m<- MASK(CHL)

#Plot Taylorgram for NO3:
NO3m    <- as.numeric(NO3m)
NO3s    <- as.numeric(NO3s)
no3     <- data.frame(obs=NO3s, mod=NO3m)
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
LNV <- LNV-log(10)

#Calculate the microphytoplankton fraction:
ln20   <- log(20**3*pi/6)
ln2    <- log(2 **3*pi/6)
mean_  <- LNV
sd_    <- sqrt(VAR)
picof.m         <- pnorm(ln2, mean_, sd_)
picof.m.ann     <- apply(picof.m, c(1,2), mean)
picof.m.ann[,1] <- NA

mod_dat_file <- paste0(nameit,'-mod-dat.pdf')
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
mtext(Varname,adj=0)

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
taylor.diagram(log(npp$vgpm), log(npp$mod), col=4, add=F, normalize=T, mar=c(2.5,3.5,2,2))
taylor.diagram(log(no3$obs), log(no3$mod), add=T,normalize=T, )
taylor.diagram(log(chl$obs), log(chl$mod), col=3, add=T, normalize=T)
taylor.diagram(log(npp$cbpm),log(npp$mod), col=5, add=T, normalize=T)
taylor.diagram(pico$obs,     pico$mod,     col=6, add=T, normalize=T)
legend(2,1.5, 
       legend = c('NO3','Chl','NPP,vgpm','NPP,cbpm','%Pico'),
       col = 2:6, pch = 16, cex=0.96)

mtext('Latitude (ºN)', side=2, outer=T)
mtext('Longitude (ºE)',side=1, outer=T, adj=0.25)
dev.off()


