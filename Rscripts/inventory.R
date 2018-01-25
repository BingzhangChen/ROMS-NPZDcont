#Total N
TN      <- NO3+PHY+ZOO+DET  #Unit: mmol m-3
TLen    <- dim(NO3)[4]  #Total number of time indices

#Total Fe
Fe2N    <- 0.0265 #Unit: nM:µM 
TFe     <- DFe+DETFe+(PHY+ZOO)*Fe2N

for (tind in 1:TLen){
    for (i in 1:L){
        for (j in 1:M){
            cff <- Hz[i,j,]/(pm[i,j]*pn[i,j])
            TN[i,j,,tind] <- TN[i,j,,tind]*cff
           TFe[i,j,,tind] <- TFe[i,j,,tind]*cff
        }
    }
}

#Calculate total N and Fe for the whole domain:
TNt  <- apply( TN, 4, sum)/1E3/1E15   #Unit: Emol (Exa moles)
TFet <- apply(TFe, 4, sum)/1E6/1E12   #Unit: Emol

pdf('Total_inventory.pdf', width = 4, height = 7, paper = 'a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(2,4,1,.1),
             mgp     = c(2.3,1,0),
             oma     = c(4,4,0,0),
             cex.lab = 1,
             pch     = 16,
             mfrow   = c(2,1))

plot(year, TNt, type = 'l',
     xlab = '', ylab = 'Total nitrogen (Emol)')
plot(year, TFet, type = 'l',
     xlab = '', ylab = 'Total iron (Pmol)')
mtext('Year', side=1, outer=T, line=0)
dev.off()


