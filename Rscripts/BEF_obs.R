#Calculate size diversity:
dat$PMU      <- NA
dat$VAR      <- NA
dat$p20      <- 1 - dat$chl20p
dat$p20      <- sapply(1:nrow(dat), function(i)min(dat$p20[i],.9999))
dat$chl0.2p  <- sapply(1:nrow(dat), function(i)max(dat$chl0.2p[i],1e-3))
dat$PMU  <- sapply(1:nrow(dat), function(i)MUsig(Y = c(dat$chl0.2p[i], dat$p20[i]), X)$PMU)
dat$VAR  <- sapply(1:nrow(dat), function(i)MUsig(Y = c(dat$chl0.2p[i], dat$p20[i]), X)$VAR)

#Calculate shannon-weaver index:
# H = -sum(PlnP)
shannon  <- function(x) {
  #x: a vector of length N
  N <- length(x)
  H <- -sum(x*log(x))
  return(H)
}
dat$SWH  <- sapply(1:nrow(dat), 
              function(i)shannon(c(dat$chl20p[i], dat$chl2p[i], dat$chl0.2p[i])))


#Combine three datasets:
sizeDat <- matrix(NA, nrow = nrow(dat) + nrow(K2size) + nrow(S1size), ncol = 4)
sizeDat <- as.data.frame(sizeDat)
colnames(sizeDat) <- c('Tchl', 'chl20p', 'chl2p', 'chl0.2p')
sizeDat[1:nrow(dat),'Tchl']                                <- dat$chltot
sizeDat[1:nrow(dat),'chl20p']                              <- dat$chl20p
sizeDat[1:nrow(dat),'chl2p']                               <- dat$chl2p
sizeDat[1:nrow(dat),'chl0.2p']                             <- dat$chl0.2p

sizeDat[(nrow(dat)+1):(nrow(dat)+nrow(K2size)),  'Tchl']   <- K2size$chlt
sizeDat[(nrow(dat)+1):(nrow(dat)+nrow(K2size)),  'chl20p'] <- K2size$microp
sizeDat[(nrow(dat)+1):(nrow(dat)+nrow(K2size)),  'chl2p']  <- K2size$nanop
sizeDat[(nrow(dat)+1):(nrow(dat)+nrow(K2size)),  'chl0.2p']<- K2size$picop

sizeDat[(nrow(dat)+nrow(K2size)+1):nrow(sizeDat),'Tchl']   <- S1size$chlt
sizeDat[(nrow(dat)+nrow(K2size)+1):nrow(sizeDat),'chl20p'] <- S1size$microp
sizeDat[(nrow(dat)+nrow(K2size)+1):nrow(sizeDat),'chl2p']  <- S1size$nanop
sizeDat[(nrow(dat)+nrow(K2size)+1):nrow(sizeDat),'chl0.2p']<- S1size$picop



#Check size diversity with Shannon-Weaver index:


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

