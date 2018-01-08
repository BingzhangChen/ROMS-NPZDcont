#Generate a diagram of discrete phytoplankton size classes
PMU <- 10
SD  <- 5

x <- c(2, 5, 10, 15, 18) 
y <- dnorm(x, mean = PMU, sd = SD)
X <- seq(0,20,0.01)
Y <- dnorm(X, mean = PMU, sd = SD)
pdf('discrete_example.pdf', width=3,height=3)
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(4,4,0.5,.2),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             mfcol   = c(1,1),
             cex.axis=1) 

plot(X, Y, type = 'n', 
     xlab = bquote('Ln Volume ('~µm^3 ~ ')'),
     ylab = 'Biomass')
rect(x, 0, x+.3, y)
dev.off()

