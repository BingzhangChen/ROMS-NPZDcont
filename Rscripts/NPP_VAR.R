#Compare integrated NPP in different treatments:
prefix1  <- '~/Roms_tools/Run/NPacS1_'
prefix2  <- '~/Roms_tools/Run/NPacS1_KTW_'
nameit   <- 'npacS'
profix2  <- paste0('/',nameit,'_dbio_avg')
profix1  <- paste0('/',nameit,'_avg')
VTR      <- c(0,0.01,0.03,0.05,0.07,0.1)
alphaG   <- VTR
N        <- length(VTR)
biof_TD  <- matrix(NA, nr = 10, nc = N)
avgf_TD  <- biof_TD
biof_KTW <- biof_TD
avgf_KTW <- biof_TD
for (i in 1:nrow(biof_TD)){
   biof_TD[i,]  <- paste0(prefix1, VTR, profix2,i,'.nc')
   avgf_TD[i,]  <- paste0(prefix1, VTR, profix1,i,'.nc')
  biof_KTW[i,]  <- paste0(prefix2, alphaG, profix2,i,'.nc') 
  avgf_KTW[i,]  <- paste0(prefix2, alphaG, profix1,i,'.nc') 
}
biof_KTW[,1] <- biof_TD[,1]
avgf_KTW[,1] <- avgf_TD[,1]
NPP_TD       <- matrix(NA, nr = 10, nc = N)
NPP_KTW      <- NPP_TD
system.time( 
    for (i in 1:6){
       NPP_TD[,i] <- sapply(1:nrow(biof_TD), function(j)integT(biof_TD[j,i]))
    }
)

NPP_TD_mean <- apply(NPP_TD, 2, mean)
NPP_TD_sd   <- apply(NPP_TD, 2, sd)

system.time( 
    for (i in 2:6){
       try(NPP_KTW[,i] <- sapply(1:nrow(biof_KTW), function(j)integT(biof_KTW[j,i])))
    }
)


pdf('TNPP_VTR.pdf',width=5, height=5,paper = 'a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(4,4,1,1),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             mfrow   = c(1,1)  )
plot(VTR, NPP_TD_mean, type = 'b',
     xlab = 'Coefficient of TD or KTW',
     ylab = bquote('Annual primary production ('*'PgC '*y^-1*')'))
low  <- NPP_TD_mean - 2*NPP_TD_sd
high <- NPP_TD_mean + 2*NPP_TD_sd
segments(VTR, low, VTR, high) 
points(alphaG,NPP_KTW, type = 'b', pch=2)
legend('topright', legend = c('TD', 'KTW'), pch=1:2)
dev.off()

#Compare diversity and NPP at different coefficients:
meanVAR  <- function(file){
   PHY.a <- Sur_mean(file, 'PHYTO')
   LNV.a <- Sur_mean(file, 'LNV')
   LNV.a <- LNV.a/PHY.a
   VAR.a <- Sur_mean(file, 'VAR')
   VAR.a <- VAR.a/PHY.a - LNV.a^2
   LNV.b <- exp(LNV.a - log(10))
   LNV.b <- (LNV.b*6/pi)**0.33333  #Unit: µm
   return(list(VAR=VAR.a, LNV=LNV.b))
}

pdf('Compare_TD_KTW.pdf',width=8, height=6,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(2,1.7,1.4,2.3),
             mgp     = c(2.3,1,0),
             oma     = c(4,4,0,0),
             cex.lab = 1.2,
             pch     = 16,
             mfrow   = c(3,3)  )

Nrep     <- nrow(avgf_TD)
avgFiles <- matrix('', nr = Nrep, nc = 3)
bioFiles <- avgFiles
for (i in 1:Nrep){
  avgFiles[i,] <- c(avgf_TD[i,1], avgf_TD[i,N], avgf_KTW[i,N])
  bioFiles[i,] <- c(biof_TD[i,1], biof_TD[i,N], biof_KTW[i,N])
}
NPPname  <- bquote('NPP (mg C '*m^-2*' '*d^-1*')')
VARname  <- expression(paste("Size variance "*'(ln '*µm^3*')'^2))
for (i in 1:ncol(avgFiles)){

    #Calculate mean for 10 replicates
    cff <- array(, c(L, M, Nrep))
    for (j in 1:Nrep){
      cff[ , , j] <- meanVAR(avgFiles[j,i])$VAR
    }
    cff <- apply(cff, c(1,2), mean)
   image2D(cff, Lon, Lat,    #Modeled size variance
             col = jet.colors(18),   zlim = c(0,8), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)
   cff <- 0.1
   if (i == 1){
      cff      <- 0
      Varname  <- bquote(italic(u) ~ "=" ~.(cff) ~","~italic(alpha[g])~"=0")
   }else if (i == 2) {
      Varname  <- bquote(italic(u) ~ "=" ~.(cff))
   }else{
      Varname  <- bquote(italic(alpha[g]) ~ "=" ~ .(cff))
   }
   mtext(Varname,adj = -0.3, cex=.8)
   mtext(VARname,adj = 1.2, cex=.8)

   cff <- array(, c(L, M, Nrep))
   for (j in 1:Nrep){
      cff[ , , j] <- meanVAR(avgFiles[j,i])$LNV
   }
   cff <- apply(cff, c(1,2), mean)

   image2D(cff, Lon, Lat,    #Modeled size variance
             col = jet.colors(18),   zlim = c(0,8), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)
   mtext('Mean size (µm)',adj = 1, cex=.8)

   cff  <- meanVAR(avgf_TD[2,1])
   cff1 <- meanVAR(avgf_TD[2,6])
   par(mfrow=c(2,1))
   image2D(cff$LNV, Lon, Lat,    #Modeled size variance
             col = jet.colors(18),   zlim = c(0,8), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")
   image2D(cff1$LNV, Lon, Lat,    #Modeled size variance
             col = jet.colors(18),   zlim = c(0,8), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")

   axis(side=1, at = lon1, labels=lon2)
   mtext('Mean size (µm)',adj = 1, cex=.8)

   NPP_ann <- integann(bioFiles[i])
   image2D(NPP_ann, Lon, Lat,    #Modeled integrated NPP
             col = jet.colors(18),   zlim = c(0,1500), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)
   mtext(NPPname,adj = 1, cex=.8)
}
mtext('Longitude (ºE)',side=1, outer=T, line=1)
mtext('Latitude  (ºN)',side=2, outer=T, line=1)
dev.off()



