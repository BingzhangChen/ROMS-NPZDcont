#Surface NO3:
NO30    <- NO3[,,Nroms,]
CHL0    <- CHL[,,Nroms,]

obs_0m <- function(stn = 'HOT', type = 'TIN'){
    file <- paste0('~/Working/FlexEFT1D/',stn,'/',stn,'_',type,'.dat')
    dat  <- read.table(file, header = T)
    dat  <- dat[dat$Depth <= 5, ]
    return(dat)
}
HOT_CHL     <- obs_0m('HOT', 'CHL')
HOT_TIN     <- obs_0m('HOT', 'TIN')
S1_CHL      <- obs_0m('S1', 'CHL')
S1_TIN      <- obs_0m('S1', 'TIN')
K2_CHL      <- obs_0m('K2', 'CHL')
K2_TIN      <- obs_0m('K2', 'TIN')
P_CHL       <- obs_0m('P', 'CHL')
P_TIN       <- obs_0m('P', 'TIN')
EQPAC_CHL   <- obs_0m('EQPAC', 'CHL')
EQPAC_TIN   <- obs_0m('EQPAC', 'TIN')

#Interpolate data at individual stations
var_0m <- function(Stn_name = 'HOT', type = 'NO3'){
    if (Stn_name == 'S1'){
        stn_lon  = 145
        stn_lat  = 30
    } else if(Stn_name == 'K2'){
        stn_lon  = 160 
        stn_lat  = 47
    } else if(Stn_name == 'HOT'){
        stn_lon  = -158
        stn_lat  = 22.75
    } else if(Stn_name == 'P'){
        stn_lon  = -145
        stn_lat  = 50
    } else if(Stn_name == 'EQPAC'){
        stn_lon  = -140
        stn_lat  = 0
    }
    stn_lon[stn_lon<0]  <- stn_lon[stn_lon<0] + 360
    gridlist <- list(x=stn_lon,y=stn_lat)
    dat <- numeric(NM) 
    if (type == 'NO3'){
        DAT0 <- NO30
    }else if (type == 'CHL'){
        DAT0 <- CHL0
    }
    for (k in 1:NM){
        d      <- list(x=Lon, y=Lat, z=DAT0[,,k] )
        dat[k] <- interp.surface.grid(d,gridlist)$z
    }
    return(dat)
}
NO3.HOT <- var_0m('HOT', 'NO3')
CHL.HOT <- var_0m('HOT', 'CHL')
NO3.S1  <- var_0m('S1',  'NO3')
CHL.S1  <- var_0m('S1',  'CHL')
NO3.K2  <- var_0m('K2',  'NO3')
CHL.K2  <- var_0m('K2',  'CHL')
NO3.P   <- var_0m('P',   'NO3')
CHL.P   <- var_0m('P',  'CHL')
NO3.EQPAC   <- var_0m('EQPAC',   'NO3')
CHL.EQPAC   <- var_0m('EQPAC',  'CHL')

pdf('jgofs_obs_mod.pdf',width=5, height=8)
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(3,2,1,2),
             mgp     = c(3,1,0),
             oma     = c(3,3,3,3),
             cex.lab = 1.2,
             cex.axis= 1.2,
             pch     = 16,
             mfrow   = c(5,2))

plot(DOYf, NO3.HOT, type = 'l', xaxt = 'n', las = 1,
    ylim = c(0, 0.2),
    xlab = '', ylab = '')
points(HOT_TIN$DOY, HOT_TIN$TIN, pch=16, cex=.5)
mtext('ALOHA, NO3',adj=0)
axis(side=1, at = DOYf, labels=1:12)

plot(DOYf, CHL.HOT, type = 'l', ylim = c(0, 0.5), xaxt = 'n', yaxt = 'n',
    xlab='', ylab='')
points(HOT_CHL$DOY, HOT_CHL$CHL, pch=16, cex=.5)
mtext('ALOHA, Chl',adj=0)
axis(side=4, las=1)
axis(side=1, at = DOYf, labels=1:12)

plot(DOYf, NO3.S1, type = 'l', xaxt = 'n', las = 1,
    ylim = c(0, 6),
    xlab = '', ylab = '')
axis(side=1, at = DOYf, labels=1:12)
points(S1_TIN$DOY, S1_TIN$TIN, pch=16, cex=.5)
mtext('S1, NO3',adj=0)
axis(side=1, at = DOYf, labels=1:12)

plot(DOYf, CHL.S1, type = 'l', ylim = c(0,2),  xaxt = 'n', yaxt = 'n',
    xlab='',ylab='')
points(S1_CHL$DOY, S1_CHL[,3], pch=16, cex=.5)
mtext('S1, Chl',adj=0)
axis(side=4, las=1)
axis(side=1, at = DOYf, labels=1:12)

plot(DOYf, NO3.K2, type = 'l',las=1,xaxt='n',
    ylim = c(0,30),
    xlab = '', ylab = '')
points(K2_TIN$DOY, K2_TIN$TIN, pch=16, cex=.5)
mtext('K2, NO3',adj=0)
axis(side=1, at = DOYf, labels=1:12)

plot(DOYf, CHL.K2, type = 'l', xlab='',ylab='',ylim=c(0,3), 
    xaxt='n',yaxt='n')
points(K2_CHL$DOY, K2_CHL[,3], pch=16, cex=.5)
mtext('K2, Chl',adj=0)
axis(side=4, las=1)
axis(side=1, at = DOYf, labels=1:12)

plot(DOYf, NO3.P, type = 'l', xaxt='n',las=1,
    ylim = c(0, max(c(NO3.P,P_TIN[,3]), na.rm=T)),
    xlab = '', ylab = '')
points(P_TIN$DOY, P_TIN$TIN, pch=16, cex=.5)
mtext('Station P, NO3',adj=0)
axis(side=1, at = DOYf, labels=1:12)

plot(DOYf, CHL.P, type = 'l',  
     xaxt ='n',yaxt='n',
     ylim = c(0, max(c(CHL.P,P_CHL[,3]), na.rm=T)))
points(P_CHL$DOY, P_CHL[,3], pch=16, cex=.5)
mtext('Station P, Chl',adj=0)
axis(side=4, las=1)
axis(side=1, at = DOYf, labels=1:12)

plot(DOYf, NO3.EQPAC, type = 'l', xaxt='n', las=1,
    ylim = c(0, max(c(NO3.EQPAC,EQPAC_TIN[,3]), na.rm=T)),
    xlab = '', ylab = '')
points(EQPAC_TIN$DOY, EQPAC_TIN$TIN, pch=16, cex=.5)
mtext('Eq. Pacific, NO3',adj=0)
axis(side=1, at = DOYf, labels=1:12)

plot(DOYf, CHL.EQPAC, type = 'l', xaxt='n',las=1,yaxt='n',xlab='',ylab='',
     ylim = c(0, max(c(CHL.EQPAC,EQPAC_CHL[,3]), na.rm=T)))
points(EQPAC_CHL$DOY, EQPAC_CHL[,3], pch=16, cex=.5)
mtext('Eq. Pacific, Chl',adj=0)
axis(side=4, las=1)
axis(side=1, at = DOYf, labels=1:12)

mtext('Month', side=1,outer=T, adj=.5)
mtext('Nitrate (µM)', side=2,outer=T, adj=.5, line=1)
mtext('Chl (µg/L)',   side=4,outer=T, adj=.5, line=1)
dev.off()


