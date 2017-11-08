#Extract mean size, size variance, and the relationships with total biomass from surface
library(ggplot2)
library(MASS)
library(plot3D)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))  #A good color

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
source('~/Roms_tools/Rscripts/Maranon_data.R')


pdf('Size_fraction_mod_obs.pdf',width=4,height=9,paper='a4')
op <- par(font.lab  = 1,
             family ="serif",
             mar    = c(2,4,1.5,2),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(3,1), pch=16,
             cex.lab=1.2,cex.axis=1.2 ) 

f2 <- kde2d(log10(CHL$Tchl), CHL$picop, n = 50, 
            lims = c(-2, log10(25), 0, 1))
image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
           xaxt = 'n',
           xlab = '', 
           ylab = '%Pico')
x2 <- c(0.01, 0.1,  1, 5,25)
x1 <- log10(x2)
axis(side=1, at = x1, labels=x2)
points(log10(SFchl$chltot), SFchl$chl0.2p, col=2,cex=.5)


#plot(as.vector(CHL$Tchl), as.vector(CHL$picop), log='x',
#     xlab= '',
#     ylab= '%Pico',
#     xlim= c(0.02,25), ylim=c(0,1),
#     pch = 16, cex = .2)
#legend('topright',legend=c('Model','Observation'),
#        cex=1.2,
#        pch=c(16,1),col=1:2)
f2 <- kde2d(log10(CHL$Tchl), CHL$LNV, n = 50, 
            lims = c(-2, log10(25), 0, 8))
image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
           xaxt = 'n',
           xlab = '', 
           ylab = expression(paste("Log mean volume "*'(ln '*µm^3*')')))
axis(side=1, at = x1, labels=x2)
points(log10(dat$Tchl), dat$PMU, col=2, cex=.5)

#plot(as.vector(CHL$Tchl), as.vector(CHL$LNV), log='x',
#     xlab= '',
#     ylab= expression(paste("Log mean volume "*'(ln '*µm^3*')')),
#     xlim= c(0.02,25), ylim=c(0,8),
#     pch = 16, cex = .2)
f2 <- kde2d(log10(CHL$Tchl), CHL$VAR, n = 50, 
            lims = c(-2, log10(25), 0, 8))
image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
           xaxt = 'n',
           xlab = '', 
           ylab = expression(paste("Size variance "*'(ln '*µm^3*')'^2)))
axis(side=1, at = x1, labels=x2)
points(log10(dat$Tchl), dat$VAR, col=2,cex=.5)

#plot(as.vector(CHL$Tchl), as.vector(CHL$VAR), log='x',
#     xlab= '',
#     ylab= expression(paste("Size variance "*'(ln '*µm^3*')'^2)),
#     xlim= c(0.02,25), ylim=c(0,8),
#     pch = 16, cex = .2)
#points(dat$Tchl, dat$VAR, col=2)
mtext('Chl (µg/L)',side=1, outer=T,line=.7, adj=0.5)
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

pdf('SizeDiversity_NPP.pdf',width=6,height=6,paper='a4')
op <- par(font.lab  = 1,
             family ="serif",
             mar    = c(4,4,1.5,3),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(1,1), pch=16, cex=.8, 
             cex.lab=1.2,cex.axis=1.2 ) 

image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
           xlab = expression("Size diversity ((Ln "*µm^3*')'^2*")"), 
           ylab = expression("Primary production (µgC"*' '*L^-1*' '*d^-1*")"))

points(dat$VAR, dat$NPP, cex=.6,pch=16, col=2)
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
