pdf('FigS5_Timeseries_2box.pdf',width=6, height=10)
op <- par( font.lab  = 1,
             family  ="serif",
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             pch     = 16,
             oma     = c(4,4,0,0),
             mfcol   = c(7,2))
zones=matrix(1:14, ncol=2, byrow=FALSE)
layout(zones, widths=c(.55,.45))
#PAR_[PAR_ > 1000] <- NA
#Plot Seasonal NO3:
ff <- 0
mo <- 1:nrow(NO3_)
for (i in 1:2){
   if (i == 1){
      par(mar=c(2,6,1.5,1))
   }else{
      par(mar=c(2,2,1.5,1))
   }
   plot(mo, NO3_[,1,i], 
        ylim = range(NO3_[,,i]),
        xlab = '', 
        ylab = 'DIN (µmol/L)', type = 'l')
   points(mo, NO3_[,2,i], type = 'l', lty = 3, col=2)
   ff <- ff+1
   mtext(paste0(letters[ff],') DIN'), adj=0, cex=.8)
   if (i == 1){
      txt <- c('Low', 'High')
      txt <- paste(txt,'diversity')
      legend('topright',txt,col=c(1,2),lty=c(1,3))
      mtext('A',adj=1)
   }else{
      mtext('B',adj=1)
   }
   plot(mo, Fer_[,1,i], ylim = range(Fer_[,,i]),
        xlab = '', 
        ylab = 'Dissolved iron (nM)', type = 'l')
   points(mo, Fer_[,2,i], type = 'l', lty = 3, col=2)
   ff <- ff+1
   mtext(paste0(letters[ff],') Dissolved iron'), adj=0,cex=.8)

   #plot(mo, PAR_[,1,i], ylim = range(PAR_[,,i],na.rm = T),
   #     xlab = '', 
   #     ylab = bquote('Light (W '~m^-2~')'), type = 'l')
   #points(mo, PAR_[,2,i], type = 'l', lty = 3, col=2)
   #ff <- ff+1
   #mtext(paste0(letters[ff],') PAR'), adj=0,cex=.8)

   if (i == 1){
       LNVrange <- c(0.5, 1)
   }else{
       LNVrange <- c(1.5, 3)
   }
   plot(mo, LNV_[,1,i], ylim = range(LNV_[,,i]),
        xlab = '', 
        ylab = 'Mean size (µm)', type = 'l')
   points(mo, LNV_[,2,i], type = 'l', lty = 3, col=2)
   #points(mo, OptESD[,i], type = 'l', lty = 1, col=3)
   ff <- ff+1
   mtext(paste0(letters[ff],') Mean size'), adj=0,cex=.8)


   Varname  <- expression(paste("Size variance "*'(ln '*µm^3*')'^2))
   plot(mo, VAR_[,1,i], ylim = range(VAR_[,,i]),
        xlab = '', 
        ylab = Varname, 
        type = 'l')
   points(mo, VAR_[,2,i], type = 'l', lty = 3, col=2)
   ff <- ff+1
   mtext(paste0(letters[ff],') Size diversity'), adj=0,cex=.8)

   Varname  <- expression(paste('Growth rate ('*d^-1*')'))
   plot(mo, mu_[,1,i], ylim = range(mu_[,,i]),
        xlab = '', 
        ylab = Varname, 
        type = 'l')
   points(mo, mu_[,2,i], type = 'l', lty = 3, col=2)
   ff <- ff + 1
   mtext(paste0(letters[ff], ') Growth rate at mean size'),
         adj=0, cex=.8)

   plot(mo, d2mu_[,1,i],ylim = range(d2mu_[,,i]), 
        xlab = '', 
        ylab = expression(paste(frac(d^2*µ,d*L^2))), 
        type = 'l')
   points(mo, d2mu_[,2,i], type='l', lty = 3, col=2)
   ff <- ff + 1
   mtext(paste0(letters[ff],') Second derivative of growth'),
         adj=0,cex=.8)
   if (i == 1){
      par(mar=c(4,6,1,1))
   }else{
      par(mar=c(4,2,1,1))
   }
   plot(mo,muAvg_[,1,i], ylim = range(muAvg_[,,i]), 
        xlab = 'Days', 
        ylab = expression(paste('Productivity ('*d^-1*')')),
        type = 'l')
   points(mo,muAvg_[,2,i], type = 'l', lty = 3, col=2)
   ff <- ff+1
   mtext(paste0(letters[ff],
                ') Community productivity'), 
         adj=0,cex=.8)
}
dev.off()
