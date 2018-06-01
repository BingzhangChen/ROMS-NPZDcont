! THis is a module based on an Optimality-based NPZD model that can couple with
! other physical (0D,1D,3D) hosts.
! The optimal phytoplankton growth part based on Pahlow et al. MEPS 2013.
! The optimal zooplankton grazing part based on Pahlow & Prowe MEPS 2010.
! Original structure taken from GOTM and adapted for FABM by Jorn Bruggeman. 
! Add the component of optimal current feeding on 20150309 
! Modify the unimodal size functional relationship of growth rate in July, 2015 
! Add 5 submodels using cpp predirectives: 1) NPZDSIMPLE (the simplest NPZD)
! 2) NPZDCHL (NPZD with dynamic Chl); 3) EFTNOSIZE (adaptive model but no size)
! 4) EFTSIZE (adaptive model with size); 5)NPZDSIZE (non-adaptive model with size)
!include "cppdefs.h"
!if defined BIOLOGY && defined PZEFTFE
      MODULE NPZDCONTFE
      use lambert
      use bio_MOD
      implicit none
      private pi
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: dg = selected_real_kind(33, 4931)
      
      !Iron:Nitrogen ratio, Aumont et al. (2003) set Fe/C = 4D-6 (mol:mol). 
      !assume redfield ratio of C/N. Times 1000 to convert umol N to nmol Fe
      real, parameter :: Fe_N = 0.0265 ! Fe:Nitrogen molar ratio
      real, parameter :: lFe  = 0.6    ! Iron ligand concentration (nM). (TOM10 P. 19)
      real, parameter :: eps  = 1d-20, 
     &                   PI   = 3.141592653529D0
      contains
!The iron subroutine must be called after the ODE subroutine
      subroutine IRONCYCLE(Temp, DET, PP_NZ,PP_ND, PP_PN,PP_DZ, Fe_scav,
     &                     DETFe,DFe)
      implicit none
      real, intent(in)    :: Temp,DET, PP_ND, PP_NZ, PP_PN, PP_DZ
      real, intent(inout) :: DFe, DETFe
      real, intent(out)   :: Fe_scav
      real                :: keq, cff
      real, parameter     :: Kscm= 5D-3,  !Minimal scavenging rate (d-1)
     &                        Ksc= 3D-2   !Particle dependent scavenging rate (umolN-1 d-1)

      DFe  = max(DFe,Femin)

!dFedt = -phytouptake + Zoo excretion - scavenging + remineralization + dust deposition
! (From TOM10 and PISCES)

      !lFe: total ligand conc. (Nikelsen et al. Geosci. Model. Dev. 2015)
!     When iron concentration is above 0.6 nM, it is scavenged by DET
      !Following TOM10Appendix, Eq. 46
      !The equilibrium constant between free iron and ligands and organic complexes. 
      keq  = 10**(17.27-1565.7/(273.15 + Temp))
      !lFe: 0.6 nM
      !Following TOM10Appendix, Eq. 45
      !Iron scavenging rate = (Basal scavenging rate + particle asorbtion)*FEprime
      cff     = 1D0+(lFe-DFe)*keq 

      Fe_scav = (Kscm + Ksc*DET*TEMPBOL(Ez,Temp)) 
     &        * (-cff + sqrt(cff**2 + 4D0*DFe*keq))/2D0/keq  

      cff = Fe_N*PP_ND !The regeneration flux from DETFe ==> DFe
      
      DFe = DFe + (cff + (PP_NZ-PP_PN) * Fe_N - Fe_scav) * dtdays
      
      ! dDETFe/dt = Zooplankton defecation and mortality + scavenging - regeneration
      DETFe = DETFe + (PP_DZ*Fe_N + Fe_scav - cff)*dtdays

      end subroutine IRONCYCLE
!----------------------------------------------------------------
! The -1 branch of the Lambert-W function
      PURE REAL FUNCTION w1(x)
         IMPLICIT NONE
         REAL, PARAMETER  :: m1=0.3361, m2=-0.0042, m3=-0.0201
         REAL, INTENT(IN) :: x
         REAL             :: s
         s  = -1.0 - LOG(-x)
         w1 = -1.0 - s - 2.0/m1*(1.0 - 1.0/(1.0 + (m1*SQRT(0.5*s)
     &                         /(1.0 + m2*s*EXP(m3*SQRT(s))))))
      END FUNCTION w1
!-----------------------------------------------------------------------
      subroutine Opt_feed(gmax,Cf,Emax,Ca,beta,ghat,RMzoo,
     &                    INGES,RES,EGES)
        implicit none
        real, intent(in) :: gmax,Cf,Emax,Ca,beta,ghat,RMzoo
        real, intent(out):: INGES,RES,EGES
        real :: g_th,aAt,aAf,AE
        ! Calculate feeding threshold of ghat (dimensionless)
          g_th = Cf/Emax/(1. - Ca)
   
        ! Calculate total foraging activity (d-1)
          aAt = gmax/beta*(-1. - w1(-(1.-g_th)/EXP(1.+beta)))
   
        ! Compute the optimal foraging activity to maximize net growth of zooplankton
        ! ghat is Ihat in Pahlow et al. (2010)
          IF (ghat .GT. g_th) THEN
      ! Optimal foraging, aAf: Specific foraging activity (d-1) = 1/h (h: biomass specific handling time)
           !With increasing food conc., optimal aAf slightly increases with food conc.
              aAf = aAt*beta/(-1. - w1(-(1.0-g_th/ghat)/EXP(1.0 +beta)))
         ! Assimilation efficiency relates with two factors
         ! 1) The difference between total activity (aAt) and foraging activity (aAf)
         ! 2) Gut passage time
         ! AE is negatively related with foraging activity
         ! because digestion enzymes decrease with increasing food
         ! Another factor, gut passange time linearly increases with handling time
            AE  = Emax*(1. - exp(-beta*(aAt/aAf - 1.)))  
         ! If not optimal grazing, set AE = 0.45 (aAf/aAt = 0.6)    
            ELSE
            aAf = 0.0
            AE  = Emax
          ENDIF
        
        !Ingestion rate (d-1) of zooplankton (gmax*Q(PT)), now should be a function of size 
        INGES = aAf * ghat 
        ! Zooplankton respiration (-> NO3)
        RES = RMzoo + Ca*AE*INGES + Cf*aAf    
        ! ZOOPLANKTON EGESTION (-> Detritus)
        EGES= INGES*(1.-AE)
        return
      end subroutine Opt_feed
!--------------------------------------------------------------------------------------
      pure real function pnorm(mean,var,y)
        implicit none
        real, intent(in) :: mean,var,y
        real             :: sumx,cff,x
        integer          :: i
        
        cff =1.
        sumx=cff
        x   =(y-mean)/sqrt(var)
        do i = 1,100
           cff = cff*x*x/float(2*i+1)
           sumx= cff+sumx
        enddo
        pnorm=0.5+1./sqrt(2.*pi)*exp(-x*x/2.)*x*sumx
        return
      end function pnorm
!--------------------------------------------------------------------------------------
      pure real function normal(mean,var,l)
        implicit none
        real, intent(in) :: mean,var,l
        normal = 1D0/sqrt(2D0*var*pi)*exp(-(l-mean)**2/var/2D0) 
        return
      end function normal
!--------------------------------------------------------------------------------------
      subroutine MONOD(Temp, PAR, NO3, mu0_,aI0_,bI0_,KN, muNet, QN,
     &                                                  theta, SI, Lno3)
      ! This subroutine calculate phytoplankton growth rates, QN and Theta based on the simplest MONOD equation
      real, intent(in)  :: Temp, PAR, NO3, mu0_, aI0_,bI0_, KN
      real, intent(out) :: QN, theta, SI, Lno3,muNet
      real              :: rmax_T
      
      ! The maximal growth rate (rmax_T) under temperature tC 
       rmax_T = mu0_*TEMPBOL(Ep,Temp)
      
      !The light limitation index (SI)
      ! Unit of aI0_C: (W m-2)-1, different from other models (aI0_Chl)
      ! Include photoinhibition (Platt et al. J Phycol 1980)
      SI = (1d0 - exp(-aI0_*PAR/rmax_T))*exp(-bI0_*PAR/rmax_T)
      
       ! Growth rate (d-1) is the product function of temperature, light, and nutrient   
       Lno3  = NO3/(NO3 + KN)
       muNet = rmax_T*Lno3*SI
       QN    = Qmin/(1d0-(1d0-Qmin/Qmax)*NO3/(NO3+Kn))
      theta  = thetamin+muNet/PAR/aI0_*(thetamax-thetamin)   !Unit: gChl/molC
      return
      end subroutine

      ! The subroutine for phytoplankton with three traits (size L,
      ! light log Iopt, Temp Topt
      subroutine citrate_PHY3(NO3,PAR_,Temp_, Fe, L, x, Z, mu,
     &     dmudL,d2mudL2, dMudx, d2Mudx2, dMudZ, d2MudZ2, theta, QN,
     &     dthetadl, d2thetadl2, dQNdL, d2QNdL2, dthetadZ, d2thetadZ2,
     &     dthetadx, d2thetadx2)
      implicit none
      real(dp), intent(in)  :: L, x, Z, NO3, PAR_,Temp_, Fe 
      real(dp), intent(out) :: mu, dmudL, d2mudL2
      real(dp), intent(out) :: dMudx, d2Mudx2, dMudZ, d2MudZ2 
      real(dp), intent(out) :: theta, QN, dthetadl, d2thetadl2
      real(dp), intent(out) :: dQNdL, d2QNdL2,dthetadZ,d2thetadZ2
      real(dp), intent(out) :: dthetadx, d2thetadx2
      real(dp) :: umax, alpha,H, G, Iopt, hat, dhatdl, d2hatdl2,cff 
      real(dp) :: fN, Kn,dfNdl,d2fNdl2,d3fNdl3, d4fNdl4
      real(dp) :: KFe,fFe
      real(dp) :: dGdL, d2GdL2, Y, dYdx, d2Ydx2, dHdx, d2Hdx2
      real(dp) :: P, ht, f, dPdZ, d2PdZ2, dfdZ, d2fdZ2
      real(dp) :: cff1, dcffdl, cff2 
      real(dg) :: cff3
      real(dp), parameter   ::  b  = .21
      real(dp), parameter   ::  Eb = .0633
      real(dp), parameter   ::  a0 = .34, k = -.47 
      real(dp), parameter   ::  w  = 14.,  T0= 15d0

      !x: ln optimal light (umol photons m-2 s-1)
      !Z: optimal temperature (ºC)
      !L: log cell volume (um^3)
      !Light term
      umax = mu0b*exp(b*x)
      alpha=   a0*exp(k*x)
      Iopt =      exp(x)

      H  = 1./(PAR_/(alpha*Iopt*Iopt) + 1./umax - 2./(alpha*Iopt) +
     &           1./(PAR_*alpha))
      !Note here suboptimal is optimal (even if I < Iopt, because umax is higher, mu is also higher)

      !Size (nutrient) term (L: lnvolume)
      hat      = exp(alphamu*L + betamu*L**2)
      dhatdl   = hat*(alphamu + 2.*betamu * L)
      d2hatdl2 = 2.*hat*betamu+hat*(alphamu+2.*betamu*L)**2

      ! Nitrogen limitation:
      Kn  = ScaleTrait(L, K0N, alphaK)
      fN  = NO3/(NO3 + Kn)  !Nitrogen limitation index

      ! Iron limitation:
      KFe = ScaleTrait(L, K0Fe, alphaFe)
      fFe = Fe/(Fe + KFe)
      

      if (fN .ge. fFe) then !Iron limitation
         call MM_2der(Fe, K0Fe,alphaFe, L, Kn,fN,dfNdl,d2fNdl2)
      else                  !Nitrogen limitation
         call MM_2der(NO3,K0N,  alphaK, L, Kn,fN,dfNdl,d2fNdl2)
      endif

      G      = hat * fN
      !Calculate the derivatives of G
      dGdL   = hat*dfNdl + dhatdl*fN
      d2GdL2 = hat*d2fNdl2 + 2.*dfNdl*dhatdl + fN*d2hatdl2

      !Calculate the derivatives of H
      Y      = 1./H
      cff    = PAR_/(alpha*Iopt*Iopt)
      dYdx   = -(k+2.)*cff
     &   - b/umax +2.*(k+1.)/(alpha*Iopt) - k/(PAR_*alpha)

      d2Ydx2 = (k+2.)**2*cff + b*b/umax 
     &   - 2.*(k+1.)**2/(alpha*Iopt) + k**2/(PAR_*alpha)

      dHdx   = -H*H*dYdx
      d2Hdx2 = (2.*H*H*H)*dYdx*dYdx - d2Ydx2*H*H

      !Temperature term
      ht      = exp(Eb*(Temp_ - T0))
      P       = 1.-((Temp_ - Z)/w)**2
      f       = ht * P
      dPdZ    = 2.*(Temp_  - Z)/(w*w)
      d2PdZ2  = -2./(w*w)
      dfdZ    = ht*dPdZ 
      d2fdZ2  = ht*d2PdZ2 

      f       = max(f, 0d0)
      Mu      = H*G*f
      dMudx   = G*f*dHdx
      d2Mudx2 = G*f*d2Hdx2
      dMudL   = f*H*dGdL
      d2MudL2 = f*H*d2GdL2
      dMudZ   = H*G*dfdZ
      d2MudZ2 = H*G*d2fdZ2

     !N:C ratio at avg. size (assumed independent of light and temp.)
      cff1    = 1.-Qmin/Qmax
      cff     = 1.-cff1*fN
      QN      = Qmin/cff

      dcffdl  = -cff1*dfNdl
      dQNdL   = cff1/cff**2 * dfNdl * Qmin
      cff2    = 1./cff**2
      cff3    = cff2/cff
      !d2QNdL2 = Qmin*cff1*(d2fNdl2/cff**2 - 2.D0/cff**3*dfNdl*dcffdl)
      d2QNdL2 = Qmin*cff1*(d2fNdl2*cff2 - 2.*cff3*dfNdl*dcffdl)  !Correct

      cff     = (thetamax - thetamin)/PAR_
      !Chl:C ratio at avg. traits
      theta      = thetamin+Mu/alpha*cff   !Unit: gChl/molC
       
      !alpha only depends on Iopt
      !Derivatives at mean size
      dthetadl   = cff/alpha*dMudL
      d2thetadl2 = cff/alpha*d2MudL2
 
      !Derivatives at mean temp.
      dthetadZ   = cff/alpha*dMudZ
      d2thetadZ2 = cff/alpha*d2MudZ2

      !Derivatives at mean Iopt.
      dthetadx   = cff*dY_Xdl(mu, alpha, dMudx, k*alpha) 
      d2thetadx2 = cff*d2Y_Xdl2(mu,alpha,dMudx, k*alpha, d2Mudx2,
     &                                             k**2*alpha)

      end subroutine 
      ! The subroutine only for phytoplankton (size as the master trait)
      subroutine PHY_NPZDCONT(NO3,PAR_,Temp_,Fe, PMU, muNet,dmudl,
     &                       d2mudl2,d3mudl3,d4mudl4,
     &                       theta,QN,dthetadL,dQNdL,d2thetadl2,d2QNdL2,
     &                       Fe_lim, fN, SI)
      implicit none
      real, intent(in)  :: PMU, NO3, PAR_,Temp_, Fe 
      real, intent(out) :: muNet, dmudl, d2mudl2
      real, intent(out) :: d3mudl3, d4mudl4
      real, intent(out) :: theta, QN, dthetadl, dQNdL
      real, intent(out) :: d2thetadl2, d2QNdL2
      real, intent(out) :: Fe_lim
      real, intent(out) :: fN, SI
      real :: mu0hat, mu0hat1
      real :: dmu0hatdl, d2mu0hatdl2, aI, cff,daI_mu0hatdl
      real :: d3mu0hatdl3, d4mu0hatdl4
      real :: dSIdl
      real :: mu0hatSI 
      real :: dmu0hatSIdl, dmu0hat_aIdl,d2mu0hat_aIdl2
      real :: d2SIdl2,d2mu0hatSIdl2
      real :: d3SIdl3,d4SIdl4
      real :: d3muIhatdl3, d4muIhatdl4
      real :: Kn,dfNdl,d2fNdl2,d3fNdl3, d4fNdl4,tf
      real :: fFe 
      real :: KFe, dcffdl, cff1
        
      tf         = TEMPBOL(Ep,Temp_)
      cff1       =alphamu + 2.* betamu * PMU
      mu0hat     =tf*mu0a*exp(alphamu*PMU + betamu*PMU**2)
      dmu0hatdl  =mu0hat*cff1
      d2mu0hatdl2=mu0hat*2.*betamu+mu0hat*cff1**2
      d3mu0hatdl3=(2.*betamu+cff1**2)*dmu0hatdl+4.*betamu*mu0hat*cff1
      d4mu0hatdl4=dmu0hatdl*8.*betamu*cff1
     &           + (2.*betamu+cff1**2)*d2mu0hatdl2+8.*betamu**2*mu0hat

      ! Initial slope of P-I curve
      aI=ScaleTrait(PMU, aI0_C, alphaI)
      
      !The light limitation index (SI)
#   ifdef AIUNIMODAL
      cff=exp(-aI*PAR_/mu0hat)
      daI_mu0hatdl=aI*(alphaI-alphamu-2D0*betamu*PMU)/mu0hat
#   else
      mu0hat1 = tf*mu0a  !To simplify the size dependence of aI
      cff=exp(-aI*PAR_/mu0hat1)
      daI_mu0hatdl=aI*alphaI/mu0hat1
#   endif
      SI   =1D0-cff
      dSIdl=cff*PAR_*daI_mu0hatdl
      
      mu0hatSI    = mu0hat*SI
      dmu0hatSIdl = mu0hatSI*cff1 + mu0hat*dSIdl

#   ifdef AIUNIMODAL
       dmu0hat_aIdl  = (dmu0hatdl - alphaI*mu0hat)/aI
      d2mu0hat_aIdl2 = mu0hat/aI*2D0*betamu
     &  + (alphamu-alphaI+2D0*betamu*PMU)*dmu0hat_aIdl
      
      daI_mu0hat2dl  = daI_mu0hatdl/mu0hat - aI/mu0hat**3*dmu0hatdl !Correct
      d2aI_mu0hatdl2 = alphaI*daI_mu0hatdl - aI/mu0hat**2*d2mu0hatdl2 
     &  - daI_mu0hat2dl*dmu0hatdl
      
      d3aI_mu0hatdl3 = d2aI_mu0hatdl2*(alphaI-alphamu-2d0*betamu*PMU)
     &  - 4d0*betamu*daI_mu0hatdl
      
      d2SIdl2 = -PAR_*dSIdl*daI_mu0hatdl+PAR_*cff*d2aI_mu0hatdl2
      d3SIdl3 = PAR_*(-2.*dSIdl*d2aI_mu0hatdl2 - d2SIdl2*daI_mu0hatdl 
     & + (1.-SI)*d3aI_mu0hatdl3)  !Correct

      daI_mu0hat3dl  = daI_mu0hat2dl/mu0hat - aI*dmu0hatdl/mu0hat**4 !Correct
      d2aI_mu0hat2dl2= d2aI_mu0hatdl2/mu0hat 
     & - daI_mu0hatdl/mu0hat**2*dmu0hatdl - daI_mu0hat3dl * dmu0hatdl 
     & - aI/mu0hat**3*d2mu0hatdl2 !Correct 

      d3aI_mu0hatdl3 = alphaI*d2aI_mu0hatdl2 - aI/mu0hat**2*d3mu0hatdl3
     & - 2d0*daI_mu0hat2dl*d2mu0hatdl2 -d2aI_mu0hat2dl2*dmu0hatdl !Correct

      daI_mu0hat4dl  = daI_mu0hat3dl/mu0hat - aI/mu0hat**5*dmu0hatdl  !Correct
      d2aI_mu0hat3dl2= d2aI_mu0hat2dl2/mu0hat 
     & - daI_mu0hat2dl/mu0hat**2*dmu0hatdl - aI/mu0hat**4*d2mu0hatdl2
     & - daI_mu0hat4dl*dmu0hatdl   !Correct
      
      d3aI_mu0hat2dl3= (d3aI_mu0hatdl3*mu0hat 
     & - daI_mu0hatdl*d2mu0hatdl2)/mu0hat**2 
     & - 2d0/mu0hat**3*(d2aI_mu0hatdl2*mu0hat - daI_mu0hatdl*dmu0hatdl)
     & * dmu0hatdl - (aI/mu0hat**3*d3mu0hatdl3
     & + 2d0*daI_mu0hat3dl*d2mu0hatdl2+d2aI_mu0hat3dl2*dmu0hatdl)  !Correct
      
      d4aI_mu0hatdl4 = alphaI*d3aI_mu0hatdl3 
     & - (aI/mu0hat**2*d4mu0hatdl4 + 3.* daI_mu0hat2dl*d3mu0hatdl3
     & + 3.*d2aI_mu0hat2dl2 * d2mu0hatdl2)- d3aI_mu0hat2dl3*dmu0hatdl  !Correct
      
      d4SIdl4 = PAR_*((1.-SI)*d4aI_mu0hatdl4 - 3.*dSIdl*d3aI_mu0hatdl3
     & - 3.*d2aI_mu0hatdl2*d2SIdl2-daI_mu0hatdl*d3SIdl3)   !Correct

#   else     
      d2SIdl2 = PAR_*alphaI/mu0hat1*aI*(cff*alphaI-dSIdl)
      d3SIdl3 = PAR_/mu0hat1*alphaI**3*aI*cff*((1.-PAR_/mu0hat1*aI)**2
     & - aI*PAR_/mu0hat1)
      d4SIdl4 = PAR_/mu0hat1*alphaI**4*aI*cff*((1.-PAR_/mu0hat1*aI)
     & * ((1.-PAR_/mu0hat1*aI)**2 - aI*PAR_/mu0hat1)
     & + aI*PAR_/mu0hat1*(2.*PAR_/mu0hat1*aI-3.))
#   endif      

      d2mu0hatSIdl2 = cff1*dmu0hatSIdl + 2.*betamu*mu0hatSI
     & + mu0hat*cff1*dSIdl + mu0hat*d2SIdl2  !Correct
      
      d3muIhatdl3 = d2mu0hatSIdl2*cff1 + dmu0hatSIdl*4.*betamu 
     & + d2mu0hatdl2*dSIdl + 2.*dmu0hatdl*d2SIdl2 + mu0hat*d3SIdl3  !Correct
      
      d4muIhatdl4 = mu0hat*d4SIdl4 + 4.*dmu0hatdl*d3SIdl3 
     & + 6.*d2mu0hatdl2*d2SIdl2 + 4.*dSIdl*d3mu0hatdl3+SI*d4mu0hatdl4  !Correct
      
      Kn  = ScaleTrait(PMU, K0N, alphaK)
      fN  = NO3/(NO3 + Kn)  !Nitrogen limitation index

      !Add iron limitation:
      if (DO_IRON) then
         KFe = ScaleTrait(PMU, K0Fe, alphaFe)
         fFe = Fe/(Fe + KFe)
      endif

      !Evaluate whether N or Fe is limiting (Liebig's law):
      !Only needs to call MM once:
      if (DO_IRON .and. (fFe < fN)) then  !Fe is limiting
         Fe_lim = 1.
         call MM_4der(Fe,  K0Fe, alphaFe, PMU, fN, dfNdl, d2fNdl2, 
     &                                             d3fNdl3, d4fNdl4)
      else
         Fe_lim = 0.
         call MM_4der(NO3, K0N,  alphaK,  PMU, fN, dfNdl, d2fNdl2,
     &                                             d3fNdl3, d4fNdl4)
      endif 

      ! Phytoplankton growth rate at the mean size:
       muNet = mu0hat*SI*fN
       dmudl = dmu0hatSIdl*fN + mu0hatSI*dfNdl
      d2mudl2=2.*dmu0hatSIdl*dfNdl+d2mu0hatSIdl2*fN+mu0hatSI*d2fNdl2 
      d3mudl3=3.*(d2mu0hatSIdl2*dfNdl+dmu0hatSIdl*d2fNdl2) 
     &       + fN*d3muIhatdl3 + mu0hatSI*d3fNdl3 !Correct
      d4mudl4=4.*d3muIhatdl3*dfNdl+6.*d2mu0hatSIdl2*d2fNdl2
     &       + 4.*dmu0hatSIdl*d3fNdl3 + fN*d4muIhatdl4 
     &       + mu0hatSI*d4fNdl4  !Correct

     !N:C ratio at avg. size
      cff1    = 1.-Qmin/Qmax
      cff     = 1.-cff1*fN
      QN      = Qmin/cff

      dcffdl  = -cff1*dfNdl
      dQNdL   = cff1/cff**2 * dfNdl * Qmin
      d2QNdL2 = Qmin*cff1*(d2fNdl2/cff**2 - 2./cff**3*dfNdl*dcffdl)  !Correct

      cff     = (thetamax - thetamin)/PAR_

      !Chl:C ratio at avg. size
      theta   = thetamin+muNet/aI*cff   !Unit: gChl/molC
      
      dthetadl = cff*dY_Xdl(muNet, aI, dmudl, aI*alphaI)  !Correct
      
      d2thetadl2 = cff*d2Y_Xdl2(muNet,aI,dmudl,aI*alphaI, d2mudl2,
     &                                               aI*alphaI**2) !Correct
 
      return
      end subroutine PHY_NPZDCONT

      pure real function dY_Xdl(Y, X, dYdl, dXdl) 
      implicit none
      real, intent(in)    :: Y, X, dYdl, dXdl
         dY_Xdl = dYdl/X - Y/X**2*dXdl
      end function

      pure real function d2Y_Xdl2(Y, X, dYdl, dXdl, d2Ydl2, d2Xdl2) 
      implicit none
      real, intent(in)    :: Y, X, dYdl, dXdl, d2Ydl2, d2Xdl2
      real(dg)            :: cff
       cff = 1./(X*X*X)
       d2Y_Xdl2 = d2Ydl2/X - 2.*dYdl*dXdl/(X**2) - Y/(X**2)*d2Xdl2 +
     &                                            2.*Y*cff*(dXdl**2)
      end function

      !Mechaelis-Mention functions and derivatives
      subroutine MM_4der(N,K0,alphaK,L,fN, dfNdl, d2fNdl2, d3fNdl3,
     &                                                      d4fNdl4)
      implicit none
      real,   intent(in)  :: N, K0, alphaK, L
      real,   intent(out) :: fN, dfNdl, d2fNdl2, d3fNdl3, d4fNdl4
      real :: Kn
         ! Half saturation constant for growth at avg. size
         Kn = ScaleTrait(L, K0, alphaK)
         fN = N/(N + Kn) ! Nutrient limitation index at avg. size
         dfNdl = -alphaK*Kn*N/(N+Kn)**2
       d2fNdl2 = -alphaK**2*N*Kn*(1./(N+Kn)**2 - 2.*Kn/(N+Kn)**3)
       d3fNdl3 = alphaK**3*N*Kn*(2.*N*Kn-(Kn-N)**2)/(Kn+N)**4
       d4fNdl4 = alphaK**4*N*Kn*(11.*Kn*N*(N-Kn)+Kn**3-N**3)/(N+Kn)**5  !Correct
      end subroutine

      !Mechaelis-Mention functions and derivatives
      subroutine MM_2der(N,K0,alphaK,L,Kn,fN, dfNdl, d2fNdl2)
      implicit none
      real(dp),   intent(in)  :: N, K0, alphaK, L
      real(dp),   intent(out) :: fN, dfNdl, d2fNdl2
      real(dp),   intent(out) :: Kn
      real(dp) :: cff
         ! Half saturation constant for growth at avg. size
         Kn = ScaleTrait(L, K0, alphaK)
         cff= N+Kn
         fN = N/(N + Kn) ! Nutrient limitation index at avg. size
         dfNdl = -alphaK*Kn*N/(N+Kn)**2
       d2fNdl2 = -alphaK**2*N*Kn*(1./(N+Kn)**2 - 2.*Kn/(cff*cff*cff))
      end subroutine

      END MODULE
!endif
