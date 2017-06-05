#Compare spatial patterns of NPP and VAR:

setwd('~/Roms_tools/Run/npacific')
pdf('Annual_mean_NPP_VAR_H.pdf',width=8, height=10,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(1,4,1.6,2),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             pch     = 16,
             mfcol   = c(4,2),
             cex.axis=1) 

oi <- 0
for (pwd1 in c('npacific','npacific_0.001')) {

   oi <- oi + 1
   if (oi == 1) {
      index <- 'high'
   } else{
      index <- 'low'
   }
   setwd(paste0('~/Roms_tools/Run/',pwd1))
   Lon <- latlon$Lon
   Lat <- latlon$Lat
   load('NPPann.Rdata')
   image2D(NPP.ann, Lon, Lat, 
             col = jet.colors(18),   zlim = c(0,60), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "Latitude (ºN)")
   mtext(paste0(letters[oi],') NPP, ',index,' diversity'), line=.5,adj=0)
   
   load('muAvg.ann.Rdata')
   image2D(muAvg.ann, Lon, Lat, 
             col = jet.colors(18),   zlim = c(0,1.5), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "Latitude (ºN)")
   mtext(paste0(letters[oi],') Productivity, ',index,' diversity'), line=.5,adj=0)
   
   load('LNVann.Rdata')
   image2D(LNV.ann, Lon, Lat, 
             col = jet.colors(18),   zlim = c(0,5), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "Latitude (ºN)")
   mtext(paste0(letters[oi],') Mean size, ',index,' diversity'), line=.5,adj=0)
   
   par(mar     = c(4,4,1.6,2))
   load('VARann.Rdata')
   cff    = VAR.ann
   VARmax = 10
   cff[,1]           = NA
   #Remove high values at bottom
   cff[,Lat <= 0][cff[,Lat <= 0] > VARmax] <- NA
   cff[cff > VARmax] = VARmax

   image2D(cff, Lon, Lat, 
             col = jet.colors(18),   zlim = c(0,VARmax), 
            xaxt = 'n',frame = F,
            xlab = "Longitude (ºE)", ylab = "Latitude (ºN)")
   mtext(paste0(letters[oi],') Size diversity, ',index,' diversity'), line=.5,adj=0)
   
   lon1 = seq(100,280,by=20)
   lon2 = lon1
   lon2[lon2>180]=lon2[lon2>180]-360
   axis(side=1, at = lon1, labels=lon2)

   par(mar = c(1,4,1.6,2))
}
dev.off()

image2D(NPP_s_ann1$data, NPP_s_ann1$Lon, NPP_s_ann1$Lat, 
          col = jet.colors(18),   zlim = c(0,60), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")
mtext('e) NPP, low diversity', line=.5,adj=0)

image2D(muAvg1.ann, NPP_s_ann$Lon, NPP_s_ann$Lat, 
          col = jet.colors(18),   zlim = c(0,1.5), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")
mtext('f) Productivity, low diversity', line=.5,adj=0)


image2D(LNV_s_ann1$data, LNV_s_ann1$Lon, LNV_s_ann1$Lat, 
          col = jet.colors(18),   zlim = c(0,5), 
         xaxt = 'n',frame = F,
         xlab = "", ylab = "Latitude (ºN)")
mtext('g) Mean size, low diversity', line=.5,adj=0)

par(mar     = c(4,4,1.6,2))
image2D(VAR_s_ann1$data, VAR_s_ann1$Lon, VAR_s_ann1$Lat, 
          col = jet.colors(18),   zlim = c(0,10), 
         xaxt = 'n',frame = F,
         xlab = "Longitude (ºE)", ylab = "Latitude (ºN)")

mtext('h) Size diversity, low diversity', line=.5,adj=0)
axis(side=1, at = lon1, labels=lon2)

#par(mar = c(1,4,1.6,2))
#image2D(NPP_s_ann2$data, NPP_s_ann$Lon, NPP_s_ann$Lat, 
#          col = jet.colors(18),   zlim = c(0,60), 
#         xaxt = 'n',frame = F,
#         xlab = "", ylab = "Latitude (ºN)")
#mtext('i) NPP, no diversity', line=.5,adj=0)
#
#image2D(muAvg2.ann, NPP_s_ann$Lon, NPP_s_ann$Lat, 
#          col = jet.colors(18),   zlim = c(0,1.5), 
#         xaxt = 'n',frame = F,
#         xlab = "", ylab = "Latitude (ºN)")
#mtext('j) Productivity, no diversity', line=.5,adj=0)
#axis(side=1, at = lon1, labels=lon2)
dev.off()

