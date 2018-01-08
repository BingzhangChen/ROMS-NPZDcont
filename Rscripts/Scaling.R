#Plot parameter scaling relationships:
setwd('~/Working/FlexEFT1D/Rscripts')
PMU <- seq(-3, 15, 0.01)

mumax <- function(x, mu0, a, b)mu0*exp(a*x+b*x**2)
Kn    <- function(x, K0,  a)    K0*exp(a*x)

pdf('Scaling.pdf', width = 6, height = 6)
op <- par(font.lab = 1,
            family ="serif", cex.axis=1.2, cex.lab=1.2,
            mar    = c(2,3.8,2.5,1.5),
            mgp    = c(2.3,1,0),
            mfcol  = c(2,2),
            oma    = c(4,4,1,0)) 

ESD <- c(0.5, 2,  5,  20, 100)
PMU_<- log(pi/6*ESD**3)
plot(PMU, mumax(PMU, 2, 0.2, -0.017),  type = 'l',  main = 'Maximal growth rate',
     xaxt='n',
     xlab='',
     ylab=expression(paste(µ[max]*' ('*d^-1*')')) )
axis(1, at=PMU_, labels=as.character(ESD))

plot(PMU, Kn(PMU, 0.2, 0.27), type = 'l', main = 'Half-saturation constant for N',
     xaxt='n',
     xlab='',
     ylab=expression(paste(K[n]*' (µmol N '*L^-1*')')) )
axis(1, at=PMU_, labels=as.character(ESD))

plot(PMU, Kn(PMU, 0.06, 0.27), type = 'l',
     main='Half-saturation constant for Fe',
     xaxt='n',
     xlab='',
     ylab=expression(paste(K[Fe]*' (nmol Fe '*L^-1*')')) )
axis(1, at=PMU_, labels=as.character(ESD))

plot(PMU, Kn(PMU, 0.02, 0.1), type = 'l', main = 'Slope of P-I curve',
     xaxt='n',
     xlab='',
     ylab=expression(paste(alpha[I]*' ('*W^-1*' '* m^2 *' '*d^-1*')'))) 
axis(1, at=PMU_, labels=as.character(ESD))
XLAB <- 'ESD (µm)'
mtext(XLAB, side=1, outer=T)
dev.off()
