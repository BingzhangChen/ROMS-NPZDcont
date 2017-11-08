#Get NPP data from vgpm:
vgpm    <- ncread(clmname,'NPP_vgpm')
eppl    <- ncread(clmname,'NPP_eppl')
cbpm    <- ncread(clmname,'NPP_cbpm')

vgpm  <- MASK(vgpm)
cbpm  <- MASK(cbpm)
eppl  <- MASK(eppl)


NO3m <- NO3_s$data
NO3m <- NO3m[,,NMo]
NO3m <- MASK(NO3m)
CHL  <- ncread(avgfile, 'CHL')
CHL  <- CHL[,,Nroms,NMo]
CHL.m<- MASK(CHL)

NO3s <- NO3clm[,,Nz,]  #Surface NO3 from WOA13
NO3s <- MASK(NO3s)

CHLs <- CHLclm[,,Nz,]  #Surface CHL from SEAWIFS
CHLs <- MASK(CHLs)

for (i in 1:length(NMo)){
   #Plot Taylorgram for NO3:
   NO3m_   <- NO3m[,,i]
   NO3s_   <- NO3s[,,i]
   no3     <- data.frame(obs=as.numeric(NO3s_), mod=as.numeric(NO3m_))
   no3[no3 <= 0] <- NA
   no3     <- na.omit(no3)
   
   CHLm_   <- CHL.m[,,i]  #Model
   CHLs_   <-  CHLs[,,i]   #Observation
   chl     <- data.frame(obs=as.numeric(CHLs_), mod=as.numeric(CHLm_))
   chl[chl <= 0] <- NA
   chl     <- na.omit(chl)

   npp   <- data.frame( mod =as.numeric(NPPI[,,i]), 
                        vgpm=as.numeric(vgpm[,,i]), 
                        cbpm=as.numeric(cbpm[,,i]))
   npp   <- na.omit(npp)

   mod_dat_file <- paste0(nameit,'-mod-dat',i,'.pdf')
   pdf(mod_dat_file,width=6, height=10,paper='a4')
   
   op <- par( font.lab  = 1,
                family  ="serif",
                mar     = c(2.5,3.5,2,2),
                mgp     = c(2.3,1,0),
                oma     = c(4,4,0,0),
                pch     = 16,
                mfrow   = c(4,2))
   
   image2D(NO3m_, Lon, Lat, 
             col = jet.colors(18),   zlim = c(0,30), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)
   mtext('Modeled NO3 (µM)',adj=0)
   
   image2D(NO3s_, Lon, Lat, 
             col = jet.colors(18),   zlim = c(0,30), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)
   mtext('WOA13 NO3 (µM)',adj=0)
   
   CHLmax <- 2
   CHLm_[CHLm_ > CHLmax] <- CHLmax
   image2D(CHLm_, Lon, Lat,    #Model CHL
             col = jet.colors(18),   zlim = c(0,CHLmax), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)
   mtext('Modeled Chl (µg/L)',adj=0)
   
   CHLs_[CHLs_ > CHLmax] <- CHLmax
   image2D(CHLs_, Lon, Lat,    #Seawifs CHL
             col = jet.colors(18),   zlim = c(0,CHLmax), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)
   mtext('SeaWIFS Chl (µg/L)',adj=0)
   
   NPPmax <- 1.6E3
   cff    <- vgpm[,,i]
   cff[cff > NPPmax] <- NPPmax
   image2D(cff, Lon, Lat, zlim =c(0,NPPmax),
             col = jet.colors(18), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)
   Varname  <- bquote('NPP, vgpm (mg C '*m^-2*' '*d^-1*')')
   mtext(Varname,adj=0)
   
   cff <- cbpm[,,i]
   cff[cff > NPPmax] <- NPPmax
   image2D(cff, Lon, Lat, zlim =c(0,NPPmax),
             col = jet.colors(18), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)
   Varname  <- bquote('NPP, cbpm (mg C '*m^-2*' '*d^-1*')')
   mtext(Varname,adj=0)
   
   cff <- NPPI[,,i]
   cff[cff > NPPmax] <- NPPmax
   cff[,1]           <- NA
   image2D(cff, Lon, Lat, zlim =c(0,NPPmax),
             col = jet.colors(18), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)
   Varname  <- bquote('Modeled NPP (mg C '*m^-2*' '*d^-1*')')
   mtext(Varname,adj=0)
   
   if (sd(log(npp$vgpm)) > sd(log(npp$vgpm))){
      taylor.diagram(log(npp$vgpm), log(npp$mod),
                     col=4, add=F, normalize=T, mar=c(2.5,3.5,2,2))
      taylor.diagram(log(npp$cbpm), log(npp$mod), 
                     col=5, add=T, normalize=T)
   }else{
      taylor.diagram(log(npp$cbpm), log(npp$mod),
                     col=5, add=F, normalize=T, mar=c(2.5,3.5,2,2))
      taylor.diagram(log(npp$vgpm), log(npp$mod), 
                     col=4, add=T, normalize=T)
   }
   taylor.diagram(log(no3$obs), log(no3$mod), add=T,normalize=T, )
   taylor.diagram(log(chl$obs), log(chl$mod), col=3, add=T, normalize=T)
   legend(2,3, 
          legend = c('NO3','Chl','NPP,vgpm','NPP,cbpm'),
          col = 2:5, pch = 16, cex=0.96)
   
   mtext('Latitude (ºN)', side=2, outer=T)
   mtext('Longitude (ºE)',side=1, outer=T, adj=0.25)
   dev.off()
}

