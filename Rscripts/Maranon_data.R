#Data from Maranon about size fractionated Chl and NPP:
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
SFchl       = dat  #Save size fractionated Chl data

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

#Add global obs. NPP data
#Read data:
file <- '~/Working/Global_PP/size_NPP.csv'
dat  <- read.csv(file)

#Remove NA data:
dat  <- dat[,c('chl0.2','chl2','chl20','NPP')]
dat  <- na.omit(dat)

#Calculate size diversity based on size-fractionated Chl:
dat$PMU <- NA   #Mean log size
dat$VAR <- NA   #Log size variance

X       <- c(2, 20)
X       <- log(X)
MUsig   <- function(Y, X){   #Y: the cumulative probability of size fraction
                             #X: Log size
    stopifnot(length(Y)==length(X))
    Yq   <- qnorm(Y)
    LM   <- lm(X ~ Yq)
    LM   <- as.numeric(coef(LM))
    PMU  <- LM[1]  # logESD
    VAR  <- LM[2]
    return(list(PMU = PMU, VAR = VAR))
}

#Calculate cumulative probability function:
dat$Tchl <- apply(dat[ , 1:3], 1, sum)
dat$p2   <- dat$chl0.2/dat$Tchl
dat$p20  <- (dat$chl0.2 + dat$chl2)/dat$Tchl
dat$p20  <- sapply(1:nrow(dat), function(i)min(dat$p20[i],.9999))

dat$PMU  <- sapply(1:nrow(dat), function(i)MUsig(Y = c(dat$p2[i], dat$p20[i]), X)$PMU)
dat$VAR  <- sapply(1:nrow(dat), function(i)MUsig(Y = c(dat$p2[i], dat$p20[i]), X)$VAR)


