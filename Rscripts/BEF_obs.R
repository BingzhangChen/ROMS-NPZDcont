#Data from Maranon about size fractionated Chl and NPP:
setwd('~/Working/Global_PP')
dat = '~/Working/Global_PP/size_Chl.csv'  #Maranon data
dat = read.csv(dat)  #For Seed I, Large fraction is > 10 um
dat = dat[!is.na(dat$chl0.2),]
dat = dat[!is.na(dat$chl20),]
dat$depth = abs(dat$depth)
dat$date  = as.character(dat$date)

#First, standardize date:
for (i in 1:nrow(dat)){

    cff = dat[i,'project']
    if (substr(cff, start=1, stop=3) != 'AMT'    && 
                                 cff != 'HK2007' &&
                                 cff != 'SCS2007'){
        dat[i,'date'] = as.Date(dat[i,'date'])
    }
}

dat0 = dat[dat$depth <=5,]

dat$chl20p  = dat$chl20/dat$chltot
dat$chl2p   = dat$chl2/dat$chltot
dat$chl0.2p = dat$chl0.2/dat$chltot

#Get data from K2 and S1:
sizeStn = function(Stn){
   file = paste0('~/Working/FlexEFT1D/',Stn,'/',Stn,'_size.dat')
   Dat  = read.table(file, header=T)
   Dat$chlt = apply(Dat[,3:ncol(Dat)], 1, sum)
   Dat$microp = Dat$SIZE10/Dat$chlt
   Dat$nanop  = Dat$SIZE3 /Dat$chlt
   Dat$picop  = (Dat$SIZE1 + Dat$SIZE_L1)/Dat$chlt
   
   Dat0 = Dat[Dat$Depth <= 10, ]
   return(Dat0)
}

K2size = sizeStn('K2')
S1size = sizeStn('S1')

pdf('Size_fraction_wK2S1.pdf',width=3,height=7,paper='a4')
op <- par(font.lab  = 1,
             family ="serif",
             mar    = c(2,4,1.5,0.5),
             mgp    = c(2.3,1,0),
             oma    = c(4,4,0,0),
             mfcol  = c(3,1), pch=16, 
             cex.lab=1.2,cex.axis=1.2 ) 

plot(dat$chltot, dat$chl20p, log='x',cex=.5,
     xlab = '',
     ylab = '%Micro',
     xlim = c(0.02,25),
     ylim = c(0,1))
points(K2size$chlt,K2size$microp,col=2)
points(S1size$chlt,S1size$microp,col=3)
legend('bottomright',legend=c('Global','K2','S1'),
        pch=16,col=1:3)

plot(dat$chltot, dat$chl2p, log='x', cex=.5,
     xlab = '',
     ylab = '%Nano',
     xlim = c(0.02,25),
     ylim = c(0,1))
points(K2size$chlt,K2size$nanop,col=2)
points(S1size$chlt,S1size$nanop,col=3)

par(mar    = c(4,4,1.5,0.5))
plot(dat$chltot, dat$chl0.2p, log='x', cex=.5,
     xlab = 'Total Chl (µg/L)',
     ylab = '%pico',
     xlim = c(0.02,25),
     ylim = c(0,1))
points(K2size$chlt,K2size$picop,col=2)
points(S1size$chlt,S1size$picop,col=3)

dev.off()

plot(Dat0$chlt, Dat0$microp, log = 'x',
     xlab = 'Total Chl (µg/L)',
     ylab = 'Percentage',
     xlim = c(0.02,25),
     ylim = c(0,1.2), pch=16, cex=.5)
points(Dat0$chlt, Dat0$chl2p,   pch=16, cex=.5, col=2)
points(dat$chltot, dat$chl0.2p, pch=16, cex=.5, col=3)
legend('topright',legend=c('%Micro','%Nano','%Pico'),
        pch=16,col=1:3)
dev.off()

