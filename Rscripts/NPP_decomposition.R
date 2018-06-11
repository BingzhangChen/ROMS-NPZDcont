source('~/Roms_tools/Rscripts/commons.R')
#This code decomposes dNPP/dVTR to PHY*QN*dmu, PHY*mu*dQN, dPHY * mu * QN at each grid (just compares the control run without TD/KTW and with maximal TD/KTW)
prefix1  <- '~/Roms_tools/Run/NPacS1_'
avgfile0 <- paste0(prefix1, '0.0/',nameit,'_avg.nc') 
avgfile1 <- paste0(prefix1, '0.1/',nameit,'_avg.nc') 

biofile0 <- paste0(prefix1, '0.0/',nameit,'_dbio_avg.nc') 
biofile1 <- paste0(prefix1, '0.1/',nameit,'_dbio_avg.nc') 

#Calculate 
PHY0     <- ncread(avgfile0, 'PHYTO',
                  start = c(1,1,    1,NMo[1]),
                  count = c(L,M,Nroms,length(NMo)))

PHY1     <- ncread(avgfile1, 'PHYTO',
                  start = c(1,1,    1,NMo[1]),
                  count = c(L,M,Nroms,length(NMo)))

QN0      <- ncread(avgfile0, 'Qbulk',
                  start = c(1,1,    1,NMo[1]),
                  count = c(L,M,Nroms,length(NMo)))

QN1      <- ncread(avgfile1, 'Qbulk',
                  start = c(1,1,    1,NMo[1]),
                  count = c(L,M,Nroms,length(NMo)))

LNV0     <- ncread(avgfile0, 'LNV',
                  start = c(1,1,    1,NMo[1]),
                  count = c(L,M,Nroms,length(NMo)))

LNV1     <- ncread(avgfile1, 'LNV',
                  start = c(1,1,    1,NMo[1]),
                  count = c(L,M,Nroms,length(NMo)))

#Convert to real mean size (log volume):
LNV0 <- LNV0/PHY0
LNV1 <- LNV1/PHY1

VAR0 <- ncread(avgfile0, 'VAR',
              start = c(1,1,    1,NMo[1]),
              count = c(L,M,Nroms,length(NMo)))

VAR1 <- ncread(avgfile1, 'VAR',
              start = c(1,1,    1,NMo[1]),
              count = c(L,M,Nroms,length(NMo)))

VAR0 <- VAR0/PHY0 - LNV0^2
VAR1 <- VAR1/PHY1 - LNV1^2

#NPP:
NPP0 <- ncread(biofile0, 'oPPt',
               start = c(1,1,    1,NMo[1]),
               count = c(L,M,Nroms,length(NMo)))

NPP1 <- ncread(biofile1, 'oPPt',
               start = c(1,1,    1,NMo[1]),
               count = c(L,M,Nroms,length(NMo)))

mu0  <- ncread(biofile0, 'omuNet',
               start = c(1,1,    1,NMo[1]),
               count = c(L,M,Nroms,length(NMo)))

mu1  <- ncread(biofile1, 'omuNet',
               start = c(1,1,    1,NMo[1]),
               count = c(L,M,Nroms,length(NMo)))

d2mu0 <- ncread(biofile0, 'od2mudl2',
               start = c(1,1,    1,NMo[1]),
               count = c(L,M,Nroms,length(NMo)))

d2mu1 <- ncread(biofile1,'od2mudl2',
               start = c(1,1,    1,NMo[1]),
               count = c(L,M,Nroms,length(NMo)))

d2mudQN0 <- ncread(biofile0, 'od2muNet_QNdl2',
               start = c(1,1,    1,NMo[1]),
               count = c(L,M,Nroms,length(NMo)))

d2mudQN1 <- ncread(biofile1,'od2muNet_QNdl2',
               start = c(1,1,    1,NMo[1]),
               count = c(L,M,Nroms,length(NMo)))

#volume for each grid (gvol)
x1 <- PHY0
x2 <- PHY0
x3 <- PHY0
x4 <- PHY0
x5 <- PHY0
x6 <- PHY0

#Calculate absolute amounts in each grid
for (i in 1:length(NMo)){
   x1[,,,i] <- PHY1[,,,i]/   QN1[,,,i]*(mu1[,,,i] - mu0[,,,i])*gvol
   x2[,,,i] <--PHY1[,,,i]*mu1[,,,i]/QN1[,,,i]**2*(QN1[,,,i] - QN0[,,,i])*gvol
   x3[,,,i] <-  mu1[,,,i] /QN1[,,,i]*(  PHY1[,,,i] -   PHY0[,,,i])*gvol
   x4[,,,i] <- PHY1[,,,i]*VAR1[,,,i]*(d2mudQN1[,,,i] - d2mudQN0[,,,i])*gvol/2
   x5[,,,i] <- VAR1[,,,i]*d2mudQN1[,,,i]*(PHY1[,,,i] -     PHY0[,,,i])*gvol/2
   x6[,,,i] <- PHY1[,,,i]*d2mudQN1[,,,i]*(VAR1[,,,i] -     VAR0[,,,i])*gvol/2
}

#Apply mask and remove the data along the southern boundary
for (i in 1:L){
    for (j in 1:M){
        if (mask[i,j] == 0){
            x1[i,j,,] <- 0
            x2[i,j,,] <- 0
            x3[i,j,,] <- 0
            x4[i,j,,] <- 0
            x5[i,j,,] <- 0
            x6[i,j,,] <- 0
        }
        #Integrate above 260 m
        w <- which(Depth[i,j,] < -260)
        x1[i,j,w,] <- 0
        x2[i,j,w,] <- 0
        x3[i,j,w,] <- 0
        x4[i,j,w,] <- 0
        x5[i,j,w,] <- 0
        x6[i,j,w,] <- 0
    }
}
x1[,1,,] <- 0
x2[,1,,] <- 0
x3[,1,,] <- 0
x4[,1,,] <- 0
x5[,1,,] <- 0
x6[,1,,] <- 0

#Sum up:
sum(x1)

#Plot out contrasts of spatial variations:
#1. NPP
file1='NPP_decomp.pdf'
pdf(file1,width=6, height=8)
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(2.5,3.5,2,2),
             mgp     = c(2.3,1,0),
             oma     = c(4,4,0,0),
             pch     = 16,
             cex.axis= 1.4,
             cex.lab = 1.4,
             mfrow   = c(1,1))
   NPPmax   <- 1.6E3
   NPP_ann0 <- integann(biofiles0, 'oPPt')
   NPP_ann1 <- integann(biofiles1, 'oPPt')
   NPP_ann0[NPP_ann0 > NPPmax] <- NPPmax
   NPP_ann1[NPP_ann1 > NPPmax] <- NPPmax

   image2D(NPP_ann0, Lon, Lat,    
              col = jet.colors(18),   zlim = c(0,8), 
             xaxt = 'n',frame = F,
             xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)

   image2D(NPP_ann1, Lon, Lat,    
              col = jet.colors(18),   zlim = c(0,8), 
             xaxt = 'n',frame = F,
             xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)


#2. PHY biomass
    cff <- apply(cff, c(1,2), mean)
   image2D(cff, Lon, Lat,    #Modeled size variance
             col = jet.colors(18),   zlim = c(0,8), 
            xaxt = 'n',frame = F,
            xlab = "", ylab = "")
   axis(side=1, at = lon1, labels=lon2)

#2. Growth rate at mean size

#3. QN at mean size

#4. Variance

#5. d2mu_QNdL2


