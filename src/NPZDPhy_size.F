      subroutine NPZDPhy_size(PMU, NO3,Fe,tC,par_, muNet, SI, Lno3,
     &                        QN, theta)
      use bio_MOD
      implicit none
      
      !INPUT PARAMETERS:
      real, intent(in)    :: PMU, NO3, tC, par_, Fe
      
      !Output parameters:
      real, intent(out)   :: muNet, SI, Lno3, QN, theta
      
      !Local variable:
      real :: Kn, mu0hat,cff,cff1,fFe,KFe
      
      tf_p   = TEMPBOL(Ep,tC)
      mu0hat = tf_p*mu0a*exp(alphamu*PMU+betamu*PMU*PMU)
      Kn     = K0N*exp(alphaK * PMU)
      
      ! Effect of light limitation
      aI   = ScaleTrait(PMU, aI0_C, alphaI)
      SI   = 1.- max(exp(-aI*par_/tf_p/mu0a),0d0)
      Lno3 = NO3/(NO3+Kn)
       KFe = ScaleTrait(PMU, K0Fe, alphaFe)
       fFe = Fe/(Fe + KFe)
      Lno3 = min(fFe, Lno3)
      muNet= mu0hat*Lno3*SI

      !N:C ratio at avg. size
      cff1 = 1.-Qmin/Qmax
      cff  = 1.-cff1*Lno3
      QN   = Qmin/cff
      cff  = (thetamax - thetamin)/par_
      theta= thetamin+muNet/aI*cff   !Unit: gChl/molC
      return
      end subroutine NPZDPhy_size
