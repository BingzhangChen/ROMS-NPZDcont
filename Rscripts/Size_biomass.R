#Extract mean size, size variance, and the relationships with total biomass from surface
nameit = 'npacific'
pwd1   = 'npacific'
library(ggplot2)
source('~/Working/FlexEFT1D/Rscripts/LO_theme.R')
setwd(paste0('~/Roms_tools/Run/',pwd1))

#Load total biomass:
load('PHY0m.Rdata')

#Load total Chl:
load('Chl0m.Rdata')
Chl0m$Dat[Chl0m$Dat <= 0] <- NA
Chl0m$Dat[!is.finite(Chl0m$Dat)] <- NA

#Load mean size data:
load('LNV0m.Rdata')
LNV0m$Dat[LNV0m$Dat <= 0] <- NA
LNV0m$Dat[!is.finite(LNV0m$Dat)] <- NA

#Load size variance data:
load('VAR0m.Rdata')
VAR0m$Dat[VAR0m$Dat <= 0] <- NA
VAR0m$Dat[!is.finite(VAR0m$Dat)] <- NA

#Calculate the microphytoplankton fraction:
Nlon   <- dim(LNV0m$Dat)[1]
Nlat   <- dim(LNV0m$Dat)[2]
NDay   <- dim(LNV0m$Dat)[3]
ln20   <- log(20**3*pi/6)
ln2    <- log(2 **3*pi/6)
mean_  <- log(LNV0m$Dat**3*pi/6)
sd_    <- sqrt(VAR0m$Dat)
microp <- 1-pnorm(ln20, mean_,sd_)
 picop <-   pnorm(ln2,  mean_,sd_)
 nanop <- 1-microp-picop

#Calculate annual mean PHY:
#PHY.ann <- apply(PHY0m$Dat, c(1,2), mean)
#PHY.ann[PHY.ann<=0] <- NA

#Calculate annual mean Chl:
Chl.ann <- apply(Chl0m$Dat, c(1,2), mean)

#Calculate annual mean picop:

microp.ann <- apply(microp, c(1,2), mean)
nanop.ann  <- apply(nanop, c(1,2), mean)
picop.ann  <- apply(picop, c(1,2), mean)
 picop.ann[picop.ann  <= 0] <- NA
 nanop.ann[nanop.ann  <= 0] <- NA
microp.ann[microp.ann <= 0] <- NA

pdf('Size_fraction_model_obs.pdf',width=5,height=8,paper='a4')
op <- par(font.lab  = 1,
             family ="serif",
             mar    = c(2,4,1.5,0.5),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(3,1), pch=16, cex=.8, 
             cex.lab=1.2,cex.axis=1.2 ) 

plot(as.vector(Chl.ann), as.vector(picop.ann), log='x',
     xlab= '',
     ylab= '%Pico',
     xlim= c(0.02,25),
     pch = 16, cex = .2)
points(dat$chltot, dat$chl0.2p, col=2)
points(K2size$chlt,K2size$picop,pch=2,col=3)
points(S1size$chlt,S1size$picop,pch=2,col=4)
legend('topright',legend=c('Model','Global','K2','S1'),
        cex=1.2,
        pch=c(16,1,2,2),col=1:4)

plot(as.vector(Chl.ann), as.vector(nanop.ann), log='x',
     xlab= '',
     ylab= '%Nano',
     xlim= c(0.02,25),
     pch = 16, cex = .2)
points(dat$chltot, dat$chl2p, col=2)
points(K2size$chlt,K2size$nanop,pch=2,col=3)
points(S1size$chlt,S1size$nanop,pch=2,col=4)

par(mar    = c(4,4,1.5,0.5))
plot(as.vector(Chl.ann), as.vector(microp.ann), log='x',
     xlab= 'Total Chl a (µg/L)',
     ylab= '%Micro',
     xlim= c(0.02,25),
     pch = 16, cex = .2)
points(dat$chltot, dat$chl20p, col=2)
points(K2size$chlt,K2size$microp,pch=2,col=3)
points(S1size$chlt,S1size$microp,pch=2,col=4)
dev.off()
#df <- data.frame(Chl   = as.vector(log(Chl.ann)),
#                 micro = as.vector(microp.ann),
#                 nano  = as.vector(nanop.ann),
#                 pico  = as.vector(picop.ann))
#
#df <- na.omit(df)
#
##Plot in contours:
#xx  <- ggplot(df, aes(x=Chl, y=pico)) +
#         stat_density2d(aes(fill = ..level..), geom="polygon")+
#         scale_x_continuous(limits=c(log(0.02),log(5)))+
#         scale_fill_gradientn(colours=rev(rainbow(100, start = 0.01, end=1)))+
#         xlab('Log Chl a (µg/L)') +
#         ylab('%Pico') +
#         #annotate("text",
#         # x      = median(enspar$x,na.rm=T),
#         # y      = median(enspar$y,na.rm=T),
#         # label  = paste('r =',round(corp,2)),
#         # family = "serif", size=4)+
#         LO_theme
#
##Add global obs. data
#xx  <- xx +
#       geom_point(data=dat, mapping=aes(x=log(chltot),y=chl0.2p))

  image2D(LNV0m$Dat[,,1], Lon, Lat, 
             col = jet.colors(18),   zlim = c(0,8), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "Latitude (ºN)")

