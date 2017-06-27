source('~/Working/Global_PP/Interpolate_WOA.R')
#Get all the data of 4D NPP (final year):
get_sur <- function(Var,file, nameit){
    NPP = readnc(Var, sourcefile = file, nameit = nameit)
    
    #Number of vertical layers:
    N   = dim(NPP$depth)[3]
    
    #Plot surface NPP:
    NPP_s = NPP$data[,,N,]

    return(list(Lon=NPP$lon, Lat=NPP$lat, days=NPP$time,data=NPP_s))
}

get_ann_sur = function(Var, file, nameit){

    NPP = readnc(Var, sourcefile = file, nameit = nameit)
    
    if (length(dim(NPP$data)) >= 4){
       #Number of vertical layers:
       N   = dim(NPP$depth)[3]
       
       #Plot surface NPP:
       NPP_s = NPP$data[,,N,]
      
    } else{
       NPP_s = NPP$data
    }
    #Annual mean NPP:
    NPP_s_ann = apply(NPP_s, c(1,2), mean) 
    
    NPP_s_ann[NPP_s_ann==0] = NA
    return(list(Lon=NPP$lon, Lat=NPP$lat, data=NPP_s_ann))
}

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))  #A good color

DrawRect <- function(x0, y0, x1, y1){
   segments(x0, y0, x1, y0)
   segments(x0, y0, x0, y1)
   segments(x1, y0, x1, y1)
   segments(x0, y1, x1, y1)
}

#Draw the dynamics of mean trait values within these two boxes:
get_box <- function(Var, file, dLon, dLat){
    #dLon and dLat must be two end points defining the range
    NPP = readnc(Var, sourcefile = file)
    
    #Number of vertical layers:
    N   = dim(NPP$depth)[3]
    
    #Plot surface NPP:
    dat = NPP$data[,,N,]
    
    #Obtain the data:
    Klon= which(NPP$lon >= dLon[1] & NPP$lon <= dLon[2]) 
    Klat= which(NPP$lat >= dLat[1] & NPP$lat <= dLat[2])
     dat= dat[Klon,Klat,]
     dat[dat==0] = NA
     dat= apply(dat, 3, function(x)mean(x,na.rm=T))
     days = sort(NPP$time)
     dat  = dat[order(NPP$time)]
    return(list(days=days,dat=dat))
}
#Draw the dynamics of mean trait values within these two boxes (for HIGH RESOLUTION):
get_boxH <- function(Var, Hdiv=T, dLon, dLat){
    if (Hdiv){
      pwd1='npacific'
    } else{
      pwd1='npacific_0.001'
    }
    pwd2=getwd()

    #Load latlon data:
    load('~/Roms_tools/Run/npacific/latlon.Rdata')
    load('~/Roms_tools/Run/npacific/Days.Rdata')

    setwd(paste0('~/Roms_tools/Run/',pwd1))
    if (Var == 'NO3'){
      load('NO30m.Rdata') 
      dat = NO30m$Dat
    } else if (Var == 'LNV'){
      load('LNV0m.Rdata')
      dat = LNV0m$Dat
    } else if (Var == 'VAR'){
      load('VAR0m.Rdata')
      dat = VAR0m$Dat
    } else if (Var == 'muNet'){
      load('mu0m.Rdata')
      dat = mu0m$Dat
    } else if (Var == 'd2mudl2'){
      load('d2mu0m.Rdata')
      dat = d2mu0m$Dat
    }
    #dLon and dLat must be two end points defining the range
    
    #Obtain the data:
    Klon= which(latlon$Lon >= dLon[1] & latlon$Lon <= dLon[2]) 
    Klat= which(latlon$Lat >= dLat[1] & latlon$Lat <= dLat[2])
     dat= dat[Klon,Klat,]
     dat[dat==0] = NA
     dat  = apply(dat, 3, function(x)mean(x,na.rm=T))
    setwd(pwd2)
     days = sort(Days)
     dat  = dat[order(days)]

    return(list(days=days, dat=dat))
}

