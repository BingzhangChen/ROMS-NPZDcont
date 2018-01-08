#Load data from the 11th year with daily resolution
TD    <- c(0,0.1)
pfix2 <- paste0('/',nameit,'_dbio_avgD.nc')
pfix1 <- paste0('/',nameit,'_avgD.nc')
prefix<- '~/Roms_tools/Run/NPacS1_'
avgFs <- paste0(prefix,TD,pfix1)
bioFs <- paste0(prefix,TD,pfix2)
source('~/Roms_tools/Rscripts/growth.R')
ND    <- 365
NO3_  <- array(NA, dim = c(ND,2,2))
TEMP  <- NO3_
PHY_  <- NO3_
Fer_  <- NO3_
LNV_  <- NO3_
VAR_  <- NO3_
PAR_  <- NO3_
 mu_  <- NO3_
d2mu_ <- NO3_
#Calculate optimal size:
Optsize <- function(NO3, temp, Fe, PAR){
   ESD  <- seq(0.59,60,0.01)
   PMU_ <- log(pi/6*ESD**3)
   r1   <- sapply(1:length(PMU_), 
                  function(x)mu(PMU = PMU_[x], 
                                PAR_= PAR, NO3 = NO3, temp = temp, Fe = Fe))
   optsize <- ESD[which.max(r1)]
   return(optsize)  #Unit: µm
}
OptESD <- matrix(NA, nr = ND, nc = 2)

#Get different variables for two regions:
for (i in 1:2){
    for (j in 1:2){
        avgfile <- avgFs[i]
        biofile <- bioFs[i]
        if (j == 1){         #Oligotrophic
           dlon <- c(x1,x2)
           dlat <- c(y1,y2)
        }else if (j == 2){   #Mesotrophic
           dlon <- c(X1,X2)
           dlat <- c(Y1,Y2)
        }
        NO3_[,i,j] <- get_box('NO3',      avgfile, dlon, dlat)
        TEMP[,i,j] <- get_box('temp',     avgfile, dlon, dlat)
        PHY_[,i,j] <- get_box('PHYTO',    avgfile, dlon, dlat)
        Fer_[,i,j] <- get_box('DFE',      avgfile, dlon, dlat)
        LNV_[,i,j] <- get_box('LNV',      avgfile, dlon, dlat)
        VAR_[,i,j] <- get_box('VAR',      avgfile, dlon, dlat)
         mu_[,i,j] <- get_box('omuNet',   biofile, dlon, dlat)
        PAR_[,i,j] <- get_box('oPAR',     biofile, dlon, dlat)
       d2mu_[,i,j] <- get_box('od2mudl2', biofile, dlon, dlat)

      if (i == 1){
         OptESD[,  j] <- sapply(1:ND,
                             function(x)Optsize(NO3=NO3_[x,i,j],
                                               temp=TEMP[x,i,j],
                                               Fe  =Fer_[x,i,j],
                                               PAR =PAR_[x,i,j]))
      }
    }
}

LNV_   <- LNV_/PHY_
VAR_   <- VAR_/PHY_ - LNV_^2
LNV_   <- exp(LNV_ - log(10))
LNV_   <- (LNV_*6/pi)**0.33333 
muAvg_ <- mu_ + VAR_*d2mu_/2 



