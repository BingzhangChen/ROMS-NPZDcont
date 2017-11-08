pdf('Opt_size.pdf', width = 6, height = 6)
op <- par(font.lab = 1,
            family ="serif", cex.axis=1.2, cex.lab=1.2,
            mar    = c(2,3.8,2.5,1.5),
            mgp    = c(2.3,1,0),
            mfcol  = c(2,2),
            oma    = c(4,4,1,0)) 

ESD  <- seq(0.5,60,0.01)
PMU_ <- log(pi/6*ESD**3)
ESD1 <- c(0.5, 2, 5, 10, 40)
PMU1 <- log(pi/6*ESD1**3)
r1   <- sapply(1:length(PMU_), 
               function(x)mu(PMU = PMU_[x], 
                             PAR_ = 50, NO3 = 0.03, temp = 25, Fe = 0.6))

r2   <- sapply(1:length(PMU_), 
               function(x)mu(PMU = PMU_[x], 
                             PAR_= 50, NO3 = 3.5, temp = 25, Fe = 0.06))

pdf('size_diagram.pdf',width=5,height=5,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             pch     = 16,
             oma     = c(4,4,0,0),
             mfcol   = c(1,1))

plot(PMU_, r1,  type = 'l', ylim=c(0,1.2),
     xaxt='n',
     xlab='ESD (µm)',
     ylab=expression(paste(µ*' ('*d^-1*')')) )
points(PMU_,r2,type='l',col=2)
axis(1, at=PMU1, labels=as.character(ESD1))
abline(v=c(0.6,.8), lty=2)
abline(v=c(1.9,2.1),lty=2, col=2)
legend('topright', c('A','B'), col=1:2,lty=2)
dev.off()
