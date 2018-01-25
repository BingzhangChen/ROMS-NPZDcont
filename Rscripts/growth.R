library(Rcpp)
sourceCpp("~/Working/FlexEFT1D/Rscripts/lambert.cpp")
mu  <- function(PMU,PAR_,NO3,temp,Fe,Ep=0.41,mu0a=0.8,K0Fe=.02,alphaFe=0.27,
                  alphamu=0.2,betamu=-0.01,aI0_C=.03,alphaI=-0.13,K0N=0.2,alphaK=0.27){
      tf     <- TEMPBOL(Ep,temp)
      mu0hat <- tf*mu0a*exp(alphamu*PMU + betamu*PMU**2)
      aI=ScaleTrait(PMU, aI0_C, alphaI)
      mu0hat1 = tf*mu0a  #To simplify the size dependence of aI
      cff=exp(-aI*PAR_/mu0hat1)
      SI   =1-cff
      mu0hatSI    = mu0hat*SI
      Kn  = ScaleTrait(PMU, K0N, alphaK)
      fN  = NO3/(NO3 + Kn)  #Nitrogen limitation index

      #Add iron limitation:
      KFe = ScaleTrait(PMU, K0Fe, alphaFe)
      fFe = Fe/(Fe + KFe)

      #Evaluate whether N or Fe is limiting (Liebig's law):
      #Only needs to call MM once:
      if (fFe < fN) {  #Fe is limiting
         Fe_lim = 1.
         fN     = fFe
      }else{
         Fe_lim = 0.
      }

      # Phytoplankton growth rate at the mean size:
       muNet = mu0hat*SI*fN
      return(muNet)
}


