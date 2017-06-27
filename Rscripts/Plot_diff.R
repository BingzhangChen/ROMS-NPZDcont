dat_lim <- function(dat, ZLIM){
   dat[dat < ZLIM[1]] = ZLIM[1]
   dat[dat > ZLIM[2]] = ZLIM[2]
   return(dat)
}
library(plot3D)
library(MASS)
#Calculate productivity log ratios between high and low diversity treatments:
pwd1      <- 'npacific_0.001'
setwd(paste0('~/Roms_tools/Run/',pwd1))
load('muAvg.ann.Rdata')
muAvg2    <- muAvg.ann

load('NPPann.Rdata')
NPP2      <- NPP.ann

pwd1      <- 'npacific'
setwd(paste0('~/Roms_tools/Run/',pwd1))
load('muAvg.ann.Rdata')
muAvg1    <- muAvg.ann
load('NPPann.Rdata')
NPP1      <- NPP.ann

muAvg.dff <- log(muAvg1/muAvg2)
NPP.dff   <- log(NPP1/NPP2)  #Difference between Low diversity and no diversity

#Load latlon:
load('latlon.Rdata')

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))  #A good color

ZLIM      <- c(-0.2,0.2)
muAvg.dff <- dat_lim(muAvg.dff, ZLIM)
NPP.dff   <- dat_lim(NPP.dff,   ZLIM)
#NPPI.dff  <- dat_lim(NPPI.dff, ZLIM)
#Plot:
x1 <- 175; x2 <- 185; y1 <- 20; y2 <- 25
X1 <- 180; X2 <- 200; Y1 <- -5; Y2 <- 5

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

#Calculate kernel densities of mu vs. NPP:
f2   <- data.frame(x=as.vector(muAvg.dff), y=as.vector(NPP.dff))
f2   <- na.omit(f2)
f2   <- kde2d(f2$x, f2$y, n = 50, lims = c(-0.18, .18, -.18, .18))

#plot(as.vector(muAvg.dff), as.vector(NPP.dff))
pdf('NPP_VAR_diff.pdf',width=6, height=8,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(4,4,2,2),
             mgp     = c(2.3,1,0),
             cex.lab = 1.4,
             lwd     = 1.5,
             mfcol   = c(2,1),
             cex.axis=1) 

image2D(muAvg.dff, latlon$Lon, latlon$Lat, 
          col = jet.colors(18),  zlim = ZLIM,
         xaxt = 'n',frame = F,
         xlab = "Longitude (ºE)", ylab = "Latitude (ºN)")
DRAW()
mtext('a)', adj=0, line = .5)
#par(mar=c(4,4,1,2))
#image2D(NPP.dff2, NPP_s_ann$Lon, NPP_s_ann$Lat, 
#          col = jet.colors(18),  zlim = ZLIM,
#         xaxt = 'n',frame = F,
#         xlab = "Longitude (ºE)", ylab = "Latitude (ºN)")
#mtext('b) Difference between low and no diversity', adj=0, line = .5)
##Select two areas for comparison:
#DRAW()
lon1 = seq(100,280,by=20)
lon2 = lon1
lon2[lon2>180]=lon2[lon2>180]-360
axis(side=1, at = lon1, labels=lon2)

axis(side=1, at = lon1, labels=lon2)

image2D(f2, col = jet.colors(20),  # zlim = c(0,.05), 
        xlab='Productivity log ratios',
        ylab='NPP log ratios')

#plot(as.vector(muAvg.dff), as.vector(NPP.dff), pch=16, cex=.5, 
abline(0,1)
mtext('b)', adj=0, line = .5)
dev.off()

