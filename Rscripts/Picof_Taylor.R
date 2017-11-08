#Eq. 2 in Ward 2015:
Picofrac <- function(Ctot, Temp, Csm = 0.2414, Ds=1, As=0.668, bs=.302) {
    Cs <- 10**(log10(Csm) + 
                log10(1-exp(-Ds/Csm*Ctot)) -
                As*exp(-bs*Temp))
    return(Cs/Ctot)
}

Nanofrac <- function(Ctot, Temp){
    Pico <- Picofrac(Ctot, Temp)
    Nano <- Picofrac(Ctot, Temp, Csm = 0.921, Ds = 1, As=0.129, bs=0.173)
    return(Nano-Pico)
}
