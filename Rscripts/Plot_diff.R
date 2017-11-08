dat_lim <- function(dat, ZLIM){
   dat[dat < ZLIM[1]] = ZLIM[1]
   dat[dat > ZLIM[2]] = ZLIM[2]
   return(dat)
}
library(MASS)
#Calculate production log ratios between high and low diversity treatments:

avgFiles  <- c(avgf_TD[1], avgf_TD[N])
bioFiles  <- c(biof_TD[1], biof_TD[N])
NPP_H     <- integann(bioFiles[2])
NPP_L     <- integann(bioFiles[1])
NPP.dff   <- log(NPP_H/NPP_L)  #Difference between high and low diversity
x1 <--165+360; x2 <- -155+360; y1 <- 22; y2 <- 28
X1 <- 190; X2 <- 225; Y1 <- -5; Y2 <- 5

DRAW <- function(){
  #Northern area
  DrawRect(x1, y1, x2, y2)
  x <- (x1+x2)/2
  y <- (y1+y2)/2
  text(x,y,'A')
  
  #Southern area
  DrawRect(X1, Y1, X2, Y2)
  X <- (X1+X2)/2
  Y <- (Y1+Y2)/2
  text(X,Y,'B')
}

muAvg <- function(avgfile, biofile){
   PHY.a   <- Sur_mean(avgfile, 'PHYTO')
    mu.a   <- Sur_mean(biofile, 'omuNet')
   d2mudL2 <- Sur_mean(biofile, 'od2mudl2')
   LNV.a   <- Sur_mean(avgfile, 'LNV')
    VAR.a  <- Sur_mean(avgfile, 'VAR')
    LNV.a  <- LNV.a/PHY.a
    VAR.a  <- VAR.a/PHY.a - LNV.a^2
    muAvg  <- mu.a + d2mudL2*VAR.a/2
    return(muAvg)
}

muAvg_L   <- muAvg(avgFiles[1], bioFiles[1])
muAvg_H   <- muAvg(avgFiles[2], bioFiles[2])
muAvg.dff <- log(muAvg_H/muAvg_L)

#Calculate kernel densities of mu vs. NPP:
f2   <- data.frame(x=as.numeric(muAvg.dff), y=as.numeric(NPP.dff))
f2   <- na.omit(f2)
f2   <- kde2d(f2$x, f2$y, n = 50, lims = c(-0.1, .1, -.1, .1))

pdf('NPP_VAR_diff.pdf',width=6, height=8,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(4,4,1.5,2),
             mgp     = c(2.3,1,0),
             cex.lab = 1.4,
             lwd     = 1.5,
             mfcol   = c(2,1),
             cex.axis=1) 

image2D(NPP.dff, Lon, Lat, 
          col = jet.colors(18), # zlim = ZLIM,
         xaxt = 'n',frame = F,
         xlab = "Longitude (ºE)", ylab = "Latitude (ºN)")
DRAW()
mtext('a) NPP log ratio of high vs. low diversity', adj=0, line = .5)
lon1 = seq(100,280,by=20)
lon2 = lon1
lon2[lon2>180]=lon2[lon2>180]-360
axis(side=1, at = lon1, labels=lon2)

image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
        xlab='Productivity log ratios',
        ylab='NPP log ratios')

#plot(as.vector(muAvg.dff), as.vector(NPP.dff), pch=16, cex=.5, 
abline(0,1)
mtext('b)', adj=0, line = .5)
dev.off()

