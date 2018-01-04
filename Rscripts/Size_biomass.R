#Extract mean size, size variance, and the relationships with total biomass from surface
library(MASS)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))  #A good color

CHLsize <- function(avgfile, biofile){
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
   
   #Load NPP:
   NPPm  <- ncread(biofile,'oPPt')
   NPP0m <- NPPm[,,Nroms,NMo]
   NPP0m <- MASK(NPP0m)
   
   #Calculate the microphytoplankton fraction:
   ln20   <- log(20**3*pi/6)
   ln2    <- log(2 **3*pi/6)
   mean_  <- LNV
   sd_    <- sqrt(VAR)
   microp <- 1-pnorm(ln20, mean_,sd_)
    picop <-   pnorm(ln2,  mean_,sd_)
    nanop <- 1-microp-picop
   
   #Use daily data for size-fractionated Chl:
   CHL <- data.frame(Tchl  = as.numeric(CHL.m),
                     LNV   = as.numeric(LNV),
                     VAR   = as.numeric(VAR),
                    microp = as.numeric(microp),
                     nanop = as.numeric(nanop),
                     picop = as.numeric(picop))
   CHL <- CHL[CHL$VAR  > 0,]
   CHL <- CHL[CHL$Tchl > 0,]
   CHL <- na.omit(CHL)
   return(CHL)
}
source('~/Roms_tools/Rscripts/Maranon_data.R')

avgfile1 <- '~/Roms_tools/Run/NPacS1_0.1/npacS_avg.nc'
avgfile2 <- '~/Roms_tools/Run/NPacS1_0/npacS_avg.nc'
biofile1 <- '~/Roms_tools/Run/NPacS1_0.1/npacS_dbio_avg.nc'
biofile2 <- '~/Roms_tools/Run/NPacS1_0/npacS_dbio_avg.nc'

CHL1_ <- CHLsize(avgfile1, biofile1)
CHL2_ <- CHLsize(avgfile2, biofile2)

pdf('Size_fraction_mod_obs2.pdf',width=6,height=7.5,paper='a4')
op <- par(font.lab  = 1,
             family ="serif",
             mar    = c(2,4.3,2,2.2),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,2,0),
             mfrow  = c(3,2), pch=16,
             cex.lab=1.8,cex.axis=1.8) 

f2 <- kde2d(log10(CHL1_$Tchl), CHL1_$picop, n = 50, 
            lims = c(-2, log10(25), 0, 1))
image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
           xaxt = 'n',
           xlab = '', 
           ylab = '%Pico')
mtext('a',adj=0,cex=1.2)
x2 <- c(0.01, 0.1,  1, 5,25)
x1 <- log10(x2)
axis(side=1, at = x1, labels=x2)
points(log10(SFchl$chltot), SFchl$chl0.2p, col=2,cex=.5)

f2 <- kde2d(log10(CHL2_$Tchl), CHL2_$picop, n = 50, 
            lims = c(-2, log10(25), 0, 1))
image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
           xaxt = 'n',
           xlab = '', 
           ylab = '')
mtext('b',adj=0,cex=1.2)
x2 <- c(0.01, 0.1,  1, 5,25)
x1 <- log10(x2)
axis(side=1, at = x1, labels=x2)
points(log10(SFchl$chltot), SFchl$chl0.2p, col=2,cex=.5)

f2 <- kde2d(log10(CHL1_$Tchl), CHL1_$LNV, n = 50, 
            lims = c(-2, log10(25), 0, 8))
image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
           xaxt = 'n',
           xlab = '', 
           ylab = expression(paste("Log mean volume "*'(ln '*µm^3*')')))
axis(side=1, at = x1, labels=x2)
mtext('c',adj=0,cex=1.2)
points(log10(dat$Tchl), dat$PMU, col=2, cex=.5)

f2 <- kde2d(log10(CHL2_$Tchl), CHL2_$LNV, n = 50, 
            lims = c(-2, log10(25), 0, 8))
image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
           xaxt = 'n',
           xlab = '', 
           ylab = '')
axis(side=1, at = x1, labels=x2)
points(log10(dat$Tchl), dat$PMU, col=2, cex=.5)
mtext('d',adj=0,cex=1.2)

#plot(as.vector(CHL$Tchl), as.vector(CHL$LNV), log='x',
#     xlab= '',
#     ylab= expression(paste("Log mean volume "*'(ln '*µm^3*')')),
#     xlim= c(0.02,25), ylim=c(0,8),
#     pch = 16, cex = .2)
f2 <- kde2d(log10(CHL1_$Tchl), CHL1_$VAR, n = 50, 
            lims = c(-2, log10(25), 0, 8))
image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
           xaxt = 'n',
           xlab = '', 
           ylab = expression(paste("Size variance "*'(ln '*µm^3*')'^2)))
axis(side=1, at = x1, labels=x2)
points(log10(dat$Tchl), dat$VAR, col=2,cex=.5)
mtext('e',adj=0,cex=1.2)

f2 <- kde2d(log10(CHL2_$Tchl), CHL2_$VAR, n = 50, 
            lims = c(-2, log10(25), 0, 8))
image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
           xaxt = 'n',
           xlab = '', 
           ylab = '')
axis(side=1, at = x1, labels=x2)
points(log10(dat$Tchl), dat$VAR, col=2,cex=.5)
mtext('f',adj=0,cex=1.2)

#plot(as.vector(CHL$Tchl), as.vector(CHL$VAR), log='x',
#     xlab= '',
#     ylab= expression(paste("Size variance "*'(ln '*µm^3*')'^2)),
#     xlim= c(0.02,25), ylim=c(0,8),
#     pch = 16, cex = .2)
#points(dat$Tchl, dat$VAR, col=2)
mtext('Chl (µg/L)',side=1, outer=T,line=.7, adj=0.5,cex=1.6)
mtext('High diversity', side=3,outer=T,adj=.25, cex=1.6)
mtext('Low  diversity', side=3,outer=T,adj=.8, cex=1.6)
#points(K2size$chlt,K2size$picop,col=3)
#points(S1size$chlt,S1size$picop,col=4)
#legend('topright',legend=c('Model','Global','K2','S1'),
#        cex=1.2,
#        pch=c(16,1,2,2),col=1:4)

#Plot size variance vs. total Chl:

#Calculate the Global NPP~size diversity relationship:

#plot(as.vector(CHL), as.vector(nanop), log='x',
#     xlab= '',
#     ylab= '%Nano',
#     xlim= c(0.02,25), ylim=c(0,1),
#     pch = 16, cex = .2)
#points(dat$chltot, dat$chl2p, col=2)
#points(K2size$chlt,K2size$nanop,pch=2,col=3)
#points(S1size$chlt,S1size$nanop,pch=2,col=4)
#
#par(mar    = c(4,4,1.5,0.5))
#plot(as.vector(CHL), as.vector(microp), log='x',
#     xlab= 'Total Chl a (µg/L)',
#     ylab= '%Micro',
#     xlim= c(0.02,25), ylim=c(0,1),
#     pch = 16, cex = .2)
#points(dat$chltot, dat$chl20p, col=2)
#points(K2size$chlt,K2size$microp,pch=2,col=3)
#points(S1size$chlt,S1size$microp,pch=2,col=4)
dev.off()

#Plot daily NPP (log scale) against VAR:
DAT2    <- data.frame(VAR = as.vector(VAR),
                      NPP = as.vector(NPP0m))
DAT2    <- na.omit(DAT2)
DAT2    <- DAT2[DAT2$VAR > 0, ]
DAT2    <- DAT2[DAT2$NPP > 0, ]
f2      <- kde2d(DAT2$VAR, DAT2$NPP, n = 50, lims = c(1, 5, .1, 36))

#Plot for old simulation runs
setwd('~/Roms_tools/Run/NPacS/NPacS')
avgfile <- 'npacS_avg11.nc'
biofile <- 'npacS_dbio_avg11.nc' 
PHY     <- ncread(avgfile, 'PHYTO')
N       <- dim(PHY)[3]
VAR     <- ncread(avgfile, 'VAR')
NPP     <- ncread(biofile, 'oPPt') 
VAR     <- VAR/PHY^2
VAR0m   <- VAR[,,N,]
NPP0m   <- NPP[,,N,]
DAT2    <- data.frame(VAR = as.vector(VAR0m),
                      NPP = as.vector(NPP0m))
DAT2    <- na.omit(DAT2)
DAT2    <- DAT2[DAT2$VAR > 0, ]
DAT2    <- DAT2[DAT2$NPP > 0, ]

f2  <- kde2d(DAT2$NPP, DAT2$VAR, n = 50, lims = c(.1, 36, 1, 8))
f3  <- kde2d(dat$NPP,   dat$VAR, n = 50, lims = c(.1, 36, 1, 8))
pdf('SizeDiversity_NPP_old.pdf',width=5,height=5,paper='a4')
op <- par(font.lab  = 1,
             family ="serif",
             mar    = c(4,4,1.5,3),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(1,1), pch=16, cex=.8, 
             cex.lab=1.4,cex.axis=1.2 ) 

image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
           ylab = expression("Size diversity ((Ln "*µm^3*')'^2*")"), 
           xlab = expression("Primary production (µgC"*' '*L^-1*' '*d^-1*")"))
#plot(DAT2$NPP, DAT2$VAR, cex=.3,pch=16,col=1)
points(dat$NPP, dat$VAR, cex=.9,pch=2, col=2)
dev.off()

#Compare biomass fraction and Chl fraction in Hong Kong waters:
#setwd('~/Working/RGC-monthly/')
#file <- 'spectra.csv'
#dat  <- read.csv(file)
#dat$micChl_frac <- dat$microChl/dat$Chla
#dat$micB_frac   <- dat$microB/dat$Btot
#
#pdf('Micro_Btot.pdf',width=4,height=6,paper='a4')
#op <- par(font.lab  = 1,
#             family ="serif",
#             mar    = c(4,4,1.5,3),
#             mgp    = c(2.3,1,0),
#             oma    = c(4,4,0,0),
#             mfcol  = c(2,1), cex=.8, 
#             cex.lab=1.2,cex.axis=1.2 ) 
#
#x = dat[dat$Stn == 'PM7',]
#y = dat[dat$Stn == 'NM3',]
#plot(x$Chla,   x$micChl_frac, pch=16,
#     xlim=range(dat$Chla, na.rm=T),ylim=c(0,1), xlab='Total Chl (µg/L)',ylab='%Micro')
#points(y$Chla, y$micChl_frac)
#legend('bottomright', legend = c('PM7','NM3'), pch=c(16,1), cex=.8)
#
#plot(x$Btot,   x$micB_frac,   pch=16,
#     xlim=range(dat$Btot, na.rm=T),ylim=c(0,1), xlab='Total Biomass (µgC/L)',ylab='%Micro')
#points(y$Btot, y$micB_frac)
#
#dev.off()
#
