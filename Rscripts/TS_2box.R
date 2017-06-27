pdf('Timeseries_2box.pdf',width=6, height=10,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mgp     = c(2.3,1,0),
             cex.lab = 1,
             pch     = 16,
             oma     = c(4,4,0,0),
             mfrow   = c(6,2))
#Plot Seasonal NO3:

par(mar=c(2,6,1,1))
plot(NO3_A1$days,NO3_A1$dat, ylim=c(0,1.5),xlab = '', ylab = 'DIN (µmol/L)', type = 'l')
points(NO3_A2$days,NO3_A2$dat, type = 'l', lty = 3, col=2)
#points(NO3_A3$days,NO3_A3$dat, type = 'l', lty = 4, col=2)
mtext('a) DIN, Hotspot A', adj=0)
txt <- c('High','Low')
txt <- paste(txt,'diversity')
legend('topright',txt,col=c(1,2),lty=c(1,3))

par(mar=c(2,2,1,1))
plot(NO3_B1$days,NO3_B1$dat, xlab = '', ylab = '', type = 'l')
points(NO3_B2$days,NO3_B2$dat, type = 'l', lty = 3, col=2)
#points(NO3_B3$days,NO3_B3$dat, type = 'l', lty = 4, col=2)
mtext('b) DIN, Hotspot B', adj=0)

par(mar=c(2,6,1,1))
plot(L_A1$days, L_A1$dat, xlab = '', ylab = 'Mean size (µm)', type = 'l')
points(L_A2$days, L_A2$dat, type = 'l', lty = 3, col=2)
mtext('c) Mean size', adj=0)

par(mar=c(2,2,1,1))
plot(L_B1$days,L_B1$dat, xlab = '', ylab = '', type = 'l')
points(L_B2$days, L_B2$dat, type = 'l', lty = 3, col=2)
mtext('d) Mean size', adj=0)

par(mar=c(2,6,1,1))
plot(Var_A1$days, Var_A1$dat, xlab = '', ylab = 'Size diversity',
     type = 'l',ylim=c(0,max(Var_A1$dat)))
points(Var_A2$days,Var_A2$dat, type = 'l', lty = 3, col=2)
mtext('e) Size diversity', adj=0)

par(mar=c(2,2,1,1))
plot(Var_B1$days,Var_B1$dat, 
    xlab = '', ylab = '', type = 'l',ylim=c(0,max(Var_B1$dat)))
points(Var_B2$days,Var_B2$dat, type = 'l', lty = 3, col=2)
mtext('f) Size diversity', adj=0)

#Plot growth rate at the mean size:
par(mar=c(2,6,1,1))
plot(mu_A1$days, log(mu_A1$dat), xlab = '', 
         ylab = expression(paste('Log growth rate ('*d^-1*')')), type = 'l')
points(mu_A2$days, log(mu_A2$dat), type = 'l', lty = 3, col=2)
#points(muAvg_A3$days, muAvg_A3$dat, type = 'l', lty = 4, col=2)
mtext('g) Growth rate at mean size', adj=0)

par(mar=c(2,2,1,1))
plot(mu_B1$days,  log(mu_B1$dat), xlab = '', ylab = '', type = 'l')
points(mu_B2$days,log(mu_B2$dat), type = 'l', lty = 3, col=2)
#points(muAvg_B3$days, muAvg_B3$dat, type = 'l', lty = 4, col=2)
mtext('h) Growth rate at mean size', adj=0)

#Plot d2mudl2 at the mean size:
par(mar=c(2,6,1,1))
plot(d2mu_A1$days,d2mu_A1$dat, 
     xlab = '', ylab = expression(paste(frac(d^2*µ,d*L^2))),ylim=c(-0.1,0), type = 'l')
points(d2mu_A2$days, d2mu_A2$dat, type='l', lty = 3, col=2)
mtext('i) Second derivative of growth', adj=0)

par(mar=c(2,2,1,1))
plot(d2mu_B1$days,d2mu_B1$dat, xlab = '', ylab = '', ylim=c(-0.1,0),type = 'l')
points(d2mu_B2$days,d2mu_B2$dat, type = 'l', lty = 3,col=2)
mtext('j) Second derivative of growth rate', adj=0)

#Plot the overall community productivity:
par(mar=c(4,6,1,1))
plot(mu_A1$days,log(muAvg_A1), xlab = 'Days', 
    ylab = expression(paste('Log productivity ('*d^-1*')')), type = 'l')
points(mu_A2$days,log(muAvg_A2), type = 'l', lty = 3, col=2)
#points(muAvg_A3$days, muAvg_A3$dat, type = 'l', lty = 4, col=2)
mtext('k) Community productivity', adj=0)

par(mar=c(4,2,1,1))
plot(mu_B1$days,  log(muAvg_B1), xlab = 'Days', ylab = '', type = 'l')
points(mu_B2$days,log(muAvg_B2), type = 'l', lty = 3, col=2)
#points(muAvg_B3$days, muAvg_B3$dat, type = 'l', lty = 4, col=2)
mtext('l) Community productivity', adj=0)

dev.off()


