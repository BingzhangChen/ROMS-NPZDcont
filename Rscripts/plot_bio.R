#Plot PHY, mean size, size variance,N:C ratio, theta, and mu:

#Annual mean surface PHY
#ZOO.a   <- Sur_mean(avgfile, 'ZOO')
PHY.a   <- Sur_mean(avgfile, 'PHYTO')
QN.a    <- Sur_mean(avgfile, 'Qbulk')
LNV.a   <- Sur_mean(avgfile, 'LNV')
VAR.a   <- Sur_mean(avgfile, 'VAR')
mu.a    <- Sur_mean(biofile, 'omuNet')
the.a   <- Sur_mean(biofile, 'otheta')
Felim.a <- Sur_mean(biofile, 'oFe_lim')
SI.a    <- Sur_mean(biofile, 'oSI')
fN.a    <- Sur_mean(biofile, 'ofN')

bio_out_file <- paste0(nameit,'-bio-out.pdf')
pdf(bio_out_file,width=5, height=8,paper='a4')

op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(1.5,3,2,2),
             mgp     = c(2.3,1,0),
             oma     = c(4,4,0,0),
             pch     = 16,
             mfrow   = c(5,2))

image2D(DFe.a, Lon, Lat,    #Model dissolved iron
          col = jet.colors(18),   zlim = c(0,1.2), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('Dissolved iron (nM)',adj = 0)

image2D(PHY.a, Lon, Lat,    #Model phyto. biomass
          col = jet.colors(18),   zlim = c(0,1.2), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('Phyto. biomass (µM N)',adj = 0)

#image2D(ZOO.a, Lon, Lat,    #Model zoo. biomass
#          col = jet.colors(18),   zlim = c(0,1.2), 
#         xaxt = 'n',frame = F,
#         xlab = "", ylab = "")
#axis(side=1, at = lon1, labels=lon2)
#mtext('Modeled ZOO biomass (µM N)',adj = 0)

#Annual mean surface QN
image2D(QN.a, Lon, Lat,    #Modeled N:C ratio
          col = jet.colors(18),   zlim = c(0,20/106), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('N:C ratio',adj = 0)

#Annual mean surface theta
thetamax <- 0.62
the.a[the.a > thetamax] <- thetamax
image2D(the.a, Lon, Lat,    #Modeled Chl:C ratio at mean size
          col = jet.colors(18), 
         zlim = c(0,thetamax), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- expression(paste("Chl:C "*' (gChl '*molC^-1*')'))
mtext(Varname,adj = 0)

#Get mean size:
#Convert to µm:
LNV.a  <- LNV.a/PHY.a
LNV.b  <- exp(LNV.a - log(10))
LNV.b  <- (LNV.b*6/pi)**0.33333  #Unit: µm

image2D(LNV.b, Lon, Lat,    #Modeled mean size
          col = jet.colors(18),   zlim = c(0,6), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
mtext('Mean size (µm)',adj = 0)

#Get size variance:
#Convert to (log(µm^3))^2:
VAR.a  <- VAR.a/PHY.a - LNV.a^2
image2D(VAR.a, Lon, Lat,    #Modeled size variance
          col = jet.colors(18),   zlim = c(0,8), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
axis(side=1, at = lon1, labels=lon2)
Varname  <- expression(paste("Size variance "*'(ln '*µm^3*')'^2))
mtext(Varname,adj = 0)

mumax <- 2
mu.a[mu.a > mumax] <- mumax
image2D(mu.a, Lon, Lat,    #Modeled growth rate at mean size
          col = jet.colors(18),   zlim = c(0,mumax), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
Varname  <- bquote('µ at mean size ( '*d^-1*')')
axis(side=1, at = lon1, labels=lon2)
mtext(Varname,adj = 0)

image2D(Felim.a, Lon, Lat,  
          col = jet.colors(18),   zlim = c(0,1), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
Varname  <- bquote('Iron limitation')
axis(side=1, at = lon1, labels=lon2)
mtext(Varname,adj = 0)

image2D(SI.a, Lon, Lat,    #Light limitation
          col = jet.colors(18),   zlim = c(0,1), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
Varname  <- bquote('Light limitation')
axis(side=1, at = lon1, labels=lon2)
mtext(Varname,adj = 0)

image2D(fN.a, Lon, Lat,    #Nutrient (N or Fe) limitation
          col = jet.colors(18),   zlim = c(0,1), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "")
Varname  <- bquote('Nutrient limitation')
axis(side=1, at = lon1, labels=lon2)
mtext(Varname,adj = 0)

mtext('Latitude (ºN)', side=2, outer=T)
mtext('Longitude (ºE)',side=1, outer=T,line=.7, adj=0.5)

dev.off()
