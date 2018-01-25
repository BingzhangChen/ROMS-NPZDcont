! $Id: bio_PZEFT.F 1458 2015-04-30 15:01:25Z bzchen $
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://www.romsagrif.org
!======================================================================
!
#include "cppdefs.h"
#if defined BIOLOGY && defined BIO_PZEFT

      subroutine biology_tile (Istr,Iend,Jstr,Jend)
      use lambert
!
! Compute biological forcing functions
!
! In this particular implementation there are 6 state variables:
! NO3, PHYTOplankton, ZOOplanknton, DETritus, Mean size, and variance.
!
      implicit none
      integer Istr,Iend,Jstr,Jend
#include "param.h"
#include "grid.h"
#include "ocean3d.h"
#include "ocean2d.h"
#include "diagnostics.h"
#include "scalars.h"
#include "forces.h"
#include "mixing.h"
      real :: mu0star,alphamu,aIstar,alphaI,A0star,alphaA,          
     &        K0star,alphaK
      real :: V0star,alphaV,Q0star,alphaQ,fV_initial
      real :: w_p,w_d,kc
      real :: gmax,kp,rdn,rpdu,rpdl,Ea,zetaN,RMchl,zetaChl,mz,Tr,alphaG
      real :: beta,Ca,Cf,Emax,RMzoo
      integer :: nutrient_uptake,grazing_formulation
      logical :: fV_opt,CF_opt
#ifdef OXYGEN
      real   CN_Z
#endif
      integer ITERMAX
      integer nsink
      real, external :: tf,w1,ScaleTrait,PenAff
#ifdef DIAGNOSTICS_BIO

      real trend_no3,trend_phy,trend_zoo,trend_det,somme
      real bilan_no3,bilan_phy,bilan_zoo,bilan_det, sinking_loss, 
     &     trend_total
      integer l, iflux
      real ThisVSinkFlux(N, NumVSinkTerms),  ! [mmol m-2 s-1], upward flux is positive
     &     ThisFlux(N, NumFluxTerms)
# ifdef OXYGEN
      real ThisGasExcFlux(NumGasExcTerms), trend_o2, bilan_o2
# endif
#endif
      parameter (
     &  ITERMAX = 5,      ! number of small implicit time step
     &  nsink   = NumVSinkTerms, ! no Chlorophyll
     &  fV_initial = 0.25,    !Initial value for the fraction of energy devoted to nutrient uptake
     &  wPstar   = 0.1,    !Sinking rate of phytoplankton (m d-1, has to be converted to seconds!!)
     &  wDET     = 1.0,    !Sinking rate of detritus (m d-1, has to be converted to seconds!!)
     &  kw      = 0.04,
     &  kc      = 0.025,   !Light attenuation of phytoplankton (mg N m-2)
     & mu0star    = 5.0
     & alphamu    = 0.1
     & aIstar     = 0.125      !Initial slope of the P-I curve (mol C gChl-1 (W m-2)-1 d-1).
     & alphaI     = 0.0
     & A0star     = 1.0        
     & alphaA     = -0.3
     & V0star     = 5.0
     & alphaV     = 0.2
     & K0star     = 0.17  !Normalization constant for half-saturation constant for nitrate (mmolN m-3), from Litchman et al. 2007
     & alphaK     = 0.27  ! Litchman et al. 2007
     & alphaG     = 2.0
     & Q0star     = 0.034
     & alphaQ     = -0.18
     &  gmax    = 2.0,     !Maximal biomass-specific ingestion rate of zooplankton at 15 ºC (d-1)
     &  rdn     = 0.03,    !The conversion rate from detritus to nutrient at reference temperature (d-1,CONVERT TO SECONDS!)
     &  rpdu    = 0.10,    !Phytoplankton mortality rate at reference temperature in the euphotic zone (d-1,CONVERT TO SECONDS!)
     &  rpdl    = 0.02,    !Phytoplankton mortality rate at reference temperature below the euphotic zone(d-1,CONVERT TO SECONDS!) 
     &  Ea        = 0.65,  !Activation energy (eV)
     &  mz        = 0.025, !Coefficient for zooplankton mortality at 15 C (d-1,CONVERT TO SECONDS!)
     &  zetaN     = 0.6,   !Cost of DIN uptake (mol C/mol N)   
     &  RMchl     = 0.1,   !Cost of chl maintenance (d-1,CONVERT TO SECONDS!) 
     &  zetaChl   = 0.4,   !Cost of photosynthesis coefficient (mol C/g Chl)
     &  beta      = 0.2,   !Assimilation coefficient of zooplankton 
     &  fV_opt    = .true.,!logic variable indicating whether phytoplankton N and carbon uptake should be optimality-based
     &  nutrient_uptake = 1,     !The model used for nutrient uptake
     &  grazing_formulation = 1, !The variable indicating the grazing formulation used by zooplankton  
     &  V0        = 5.0,  !Maximal nutrient uptake rate (mol N/mol C/d, convert to seconds!)
     &  A0        = 1.0,  !Maximal nutrient uptake affinity (m3/mmol C/d, convert to seconds!)
     &  Tr        = 15.0, !Reference temperature (ºC)
     &  CF_opt    = .true.,!logic variable indicating whether zooplankton feeding should be optimality-based
     &  Ca        = 0.3,  !Cost of assimilation coefficient
     &  Cf        = 0.3,  !Effective cost of foraging coefficient
     &  Emax      = 0.99, !Maximal assimilation efficiency
     &  RMZoo     = 0.16 )!Specific maintenance respiration (d-1)
        

      integer i,j,k, ITER, iB
      real NO3(N), PHY(N), ZOO(N), DET(N),PMU(N),VAR(N),TEM(N),FC(0:N),
     &     PARsup,par,I_0, theta, rmax_T, thetaHat, dtdays,attn,
     &     Q, larg, I_min, Ith, rpd,fpar,rIhat,vNhat,fA,ZINT,fV,
     &     Rchl,mu,dn,
     &     dp,dd,dz,ghat,g_th,aAt,aAf,INGES,AE,RES,NGR,EGES,
     &     SB(N,nsink),dSB(0:N,nsink),wSB(nsink),cff,cff1,cff2,cff6,
     &     TwPHYT, TwDET, Trmax,Tgmax,Talpha,Trdn,Trpdu,Trpdl,Tmz,
     &     aR,aL,cu
      real :: I_0,tC, par,I_min,I_zero   !Environmental variables
      real :: V0hat,Kn,Lmin,A0hat,dA0hatdl,d2A0hatdl2,fA                          
      real :: VNhat,dVNhatdl,d2VNhatdl2 ! Nutrient uptake variables
      real :: mu0hat,muIhat,mu0hatSI,dmu0hatSIdl,d2mu0hatSIdl2
      real :: dmuIhatdl,d2muIhatdl2,muNet,dmuNetdl,d2muNetdl2 ! Growth rate variables
      real :: fV,dfVdl,d2fVdl2  !fV
      real :: ZINT,dZINdl,d2ZINdl2 !ZINT
      real :: aI,SI,dSIdl         !Light dependent component
      real :: RMchl,Theta,ThetaAvg,ThetaHat,dThetaHatdl,d2ThetaHatdl2 !Chl related variables
      real :: dPMUdt,dVARdt,dPHYdt,dNO3dt,dZOOdt,dDETdt        ! For ODE
      real :: Qopt,Qs,QNavg,dQoptdl,d2Qoptdl2  ! cell quota related variables
      real :: rpd,Cpp,Chltot
      real :: INGES,AE,RES,NGR,EGES,gbar,dgdlbar,d2gdl2bar,gtot  ! Grazing related variables
      real            :: larg
      real, parameter :: secs_pr_day = 86400., errr = 1E-10, Penfac=3.0,pi=3.1415926535

!# ifdef OXYGEN
!     &      , O2(N), tem(N), sal(N), den(N)
!     &      , O2satu_loc, Kv_O2_loc
!     &      , eos80, u10_loc, Sc
!     
!#define OCMIP_OXYGENSAT
!
!#  ifdef OCMIP_OXYGENSAT
!     &       , o2sato   ! OCMIP function, calculates O2 saturation
!#  else /* OCMIP_OXYGENSAT */
!     &       , satpc    ! oxygen saturation in % (calculated, but unused)
!     &       , AOU      ! Apparent oxygen utilization (calc., but unused)
!#  endif /* OCMIP_OXYGENSAT */
!# endif /* OXYGEN */
!# if defined OXYGEN || defined DIAGNOSTICS_BIO
!     &      , dtsec     ! length of time step in seconds (for gas exchange)
!# endif
!# if defined DIAGNOSTICS_BIO
!     &      , LastVSinkFlux,ColumnMassOld(NumVSinkTerms)
!     &      , ColumnMassNew(NumVSinkTerms)
!# endif /* DIAGNOSTICS_BIO */
!
# include "compute_auxiliary_bounds.h"
!
      dtdays= dt/(24.*3600.*float(ITERMAX))  ! time step as fraction of day.
# if defined DIAGNOSTICS_BIO || defined OXYGEN
      dtsec = dt / float(ITERMAX)           ! time step in seconds
# endif /* DIAGNOSTICS_BIO || OXYGEN */
      
! Convert the rates into the time scale of dtdays.
      TwPHYT = wPHYT*dtdays
      TwDET = wDET*dtdays
      Trmax= rmax*dtdays
      Tgmax= gmax*dtdays
      Talpha=alpha*dtdays
      Trdn = rdn*dtdays
      Trpdu= rpdu*dtdays
      Trpdl= rpdl*dtdays
      Tmz  = mz*dtdays
      TRMchl=RMchl*dtdays
      TRMzoo=RMzoo*dtdays
      TV0  = V0*dtdays
      TA0  = A0*dtdays
!
! Since the following solver is iterative to achieve implicit
! discretization of the biological interaction, two time slices are
! required, BIO where BIO is understood as vector of
! biological state variables: BIO=[NO3,Phyt,Zoo,Det]. Assume
! that the iterations converge, the newly obtained state variables
! satisfy equations
!
!           BIO = BIO + dtdays * rhs(BIO)
! 
! where rhs(BIO) is the vector of biological r.h.s. computed at
! the new time step. During the iterative procedure a series of
! fractional time steps is performed in a chained mode (splitting
! by different biological conversion processes) in sequence NO3 -- 
! Phyt -- Zoo -- Det, that is the main food chain. In all 
! stages the concentration of the component being consumed is
! treated in fully implicit manner, so that the algorithm guarantees
! non-negative values, no matter how strong is the concentration of
! active consuming component (Phyto or Zoo).
!
! The overall algorithm, as well as any stage of it is formulated
! in conservative form (except explicit sinking) in sense that the
! sum of concentration of all five components is conserved.
!

/*
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
*/
#  define I_RANGE Istr,Iend
#  define J_RANGE Jstr,Jend

      

      do j=J_RANGE
        do i=I_RANGE
#ifdef DIAGNOSTICS_BIO
! Reset the biogeochemical fluxes. This is necessary because the
! biological routine uses multiple. time steps for each physical time
! step.
          do k=1,N
            do l=1,NumFluxTerms
              bioFlux(i,j,k,l) = 0.0
            enddo
          end do
          do k=0,N
            do l=1,NumVSinkTerms
              bioVSink(i,j,k,l) = 0.0
            enddo
          enddo
# ifdef OXYGEN
          do l=1,NumGasExcTerms
            GasExcFlux(i,j,l) = 0.0
          enddo
# endif
#endif /* DIAGNOSTICS_BIO */
!
! Extract biological variables from tracer arrays; place them into
! scratch variables; restrict their values to be positive definite.
! t: tracer array!

      do k=1,N                                  ! N: the number of vertical layers
        NO3(k) =max(t(i,j,k,nnew,ino3_),0.)   ! Nitrate
        PHY(k) =max(t(i,j,k,nnew,iphy1),0.)   ! Phytoplankton
        ZOO(k) =max(t(i,j,k,nnew,izoo1),0.)   ! Zooplankton
        DET(k) =max(t(i,j,k,nnew,idet1),0.)   ! Detritus
        PMU(k) =    t(i,j,k,nnew,ipmu1)       ! Log cell volume [µm**3]
        VAR(k) =max(t(i,j,k,nnew,ivar1),0.)   ! Variance
        CHL(k) =max(t(i,j,k,nnew,ichl1),0.)
        TEM(k) =max(t(i,j,k,nnew,itemp),0.)   ! temperature; [deg. C]
# ifdef OXYGEN
        sal(k) =max(t(i,j,k,nnew,isalt),0.)      ! salinity; [PSU]
#  ifndef OCMIP_OXYGENSAT
        den(k) =1000.+ rho1(i,j,k)   		 ! potential density; [kg m-3]
#  endif
        O2(k)  =max(t(i,j,k,nnew,iO2),0.)        ! Oxygen;  [mmol O2 m-3]
# endif 
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO ITER=1,ITERMAX      !--> Start internal iterations to achieve
                             !    nonlinear backward-implicit solution.

        I_0 = srflx(i,j)*rho0*Cp*0.43              ! Surface par
! srflx obtained from the Routine set_forces.F
! rho0: average density of seawater
! Cp  : specific heat of seawater
        I_0 = min(I_0,errr)  ! Force I_0 to be positive
        I_min = 0.01*I_0

!        if (I_0 .gt. 0.) then
!
!   *** SUN IS UP ***
!
! Calulate mu: Set Photosynthetically Available Radiation (PAR) at
! surface from solar radiation x 0.43. Then, within each grid box
! compute attenuation coefficient based on the concentration of
! Phytoplankton inside the grid box, and attenuate PAR from surface
! down (thus, PAR at certain depth depends on the whole distribution
! of Phytoplankton above). To compute aJ, one needs PAR somewhat in
! the middle of the gridbox, so that attenuation "attn" corresponds
! to half of the grid box height, while PAR is multiplied by it
! twice: once to get it in the middle of grid-box and once the
! compute on the lower grid-box interface;
!
        par = I_0
        do k=N,1,-1     !<-- irreversible

          attn=exp(-0.5*(kw + kc*CHL(k))*     ! Attenuation in this grid
     &            (z_w(i,j,k)-z_w(i,j,k-1)))   ! z_w: depth of the grid (negative values!).

          PARsup = par * attn       ! Light level at the middle of the grid

         ! Loss rate of phytoplankton to detritus depends on local light intensity.
          if (PARsup .ge. I_min) then
              rpd = rpdu*tf(Ea,TEM(k))*dtdays
           else
              rpd = rpdl*tf(Ea,TEM(k))*dtdays
          endif
      mu0hat = dtdays*ScaleTrait(PMU(k),tf(Ea,TEM(k))*mu0star,alphamu)
      V0hat  = dtdays*ScaleTrait(PMU(k),tf(Ea,TEM(k))*V0star, alphaV)
      aI     = dtdays*ScaleTrait(PMU(k),            aIstar,alphaI)   ! Initial slope of P-I curve
      Qs     =        ScaleTrait(PMU(k),            Q0star,alphaQ)/2.
      RMChl  = dtdays*tf(self,tC)*RMChl

      ! Define I_min as 1% of I_0 (bottom of euphotic layer)
      I_min = 0.01*I_0

      ! Threshold irradiance and RMchl is set temperature dependent
      Ith = zetaChl*RMchl/aI  

   ! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)
      IF( fV_opt .EQ. .true. ) THEN 
           if( PARsup .gt. Ith ) then
            larg = (1.0 + RMchl/mu0hat) 
     &           * exp(1.0 + aI*PARsup/(mu0hat*zetaChl))

            ThetaHat=1.0/zetaChl+(1.0-WAPR(larg,0,0))*mu0hat/(aI*PARsup)
           ! First derivative of ThetaHat
           dThetaHatdl = (alphamu - alphaI)                 
     &  *         (ThetaHat - (1/zetaChl)/(1 + WAPR(larg,0,0)))   
     &  +        alphamu*( WAPR(larg,0,0)/(1 + WAPR(larg,0,0)))  
     &  *         (mu0hat/(aI*PARsup))*RMchl/(mu0hat + RMchl)   
           ! Second derivative of ThetaHat 
         d2ThetaHatdl2 = (alphamu-alphaI) * dThetaHatdl                  
     &  -          ( (aI*PARsup/mu0hat)                                    
     &  *            ( (alphaI-alphamu) /zetaChl)**2 )    
     &  *             WAPR(larg,0,0)/(1 + WAPR(larg,0,0))**3                         
     &  -             2*alphamu                                                  
     &  *            ( (alphaI-alphamu)/zetaChl)         
     &  *            ( RMchl/(RMchl + mu0hat) )                     
     &  *            ( mu0hat/(aI*PARsup) )                               
     &  *              WAPR(larg, 0, 0)/(1 + WAPR(larg, 0, 0))**3                        
     &  +              alphamu                                                    
     &  *            ( (alphamu - alphaI))                                  
     &  -               alphamu       
     &  *            ( mu0hat/(RMchl + mu0hat) )  
     &  *            (  RMchl/(RMchl + mu0hat) )                           
     &  *            ( mu0hat / (aI*PARsup) )                                                    
     &  *            WAPR(larg,0,0)/(1.+WAPR(larg, 0, 0))
            
           else
              ThetaHat      = 1e-6  !  a small positive value 
              dThetaHatdl   = errr
              d2ThetaHatdl2 = errr
           endif
        ELSE
           ThetaHat = 0.6
           dThetaHatdl   = errr
           d2ThetaHatdl2 = errr

        ENDIF
   
       !The light limitation index (SI)
       SI = 1.-exp(-aI*PARsup*ThetaHat/mu0hat)
       dSIdl = ( (alphaI-alphamu)*ThetaHat + dThetaHatdl)
    &      *(1. - SI)*aI*PARsup/mu0hat    !confirmed
    
       ! Light dependent growth rate 
       ! (needs to take into account the cost of dark and light-dependent chl maintenance)
       mu0hatSI = mu0hat*SI  ! Gross specific carbon uptake (photosynthesis)
       muIhat   = mu0hatSI - (mu0hat*SI+RMchl)*zetaChl*ThetaHat ! Net specific carbon uptake
    
       dmu0hatSIdl = alphamu*mu0hatSI
     &      +      aI*PARsup*(1.0 - SI)
     &      *      ((alphaI - alphamu)*ThetaHat 
     &      +     dThetaHatdl)
    !
         d2mu0hatSIdl2 = alphamu*dmu0hatSIdl
     & + PARsup*(alphaI-alphamu)
     & * (aI*(1.-SI)*dThetaHatdl-aI*ThetaHat*dSIdl
     & + (1.-SI)*ThetaHat*alphaI*aI) + d2ThetaHatdl2
       
        dmuIhatdl =  dmu0hatSIdl                      
     &    * (1.-zetaChl*ThetaHat*(2.*mu0hatSI+RMchl))
     &    - mu0hatSI*zetaChl*(RMchl+mu0hatSI)*dThetaHatdl
     !
           d2muIhatdl2 = -dmu0hatSIdl*zetaChl
     &      * (2.*ThetaHat*dmu0hatSIdl+(2.*mu0hatSI+RMchl)*dThetaHatdl)
     &      + d2mu0hatSIdl2*(1.-zetaChl*ThetaHat*(2.*mu0hatSI+RMchl)) 
     &      - zetaChl*(RMchl*mu0hatSI*d2ThetaHatdl2                    
     &      +               RMchl*dThetaHatdl*dmu0hatSIdl
     &      +               mu0hatSI**2*d2ThetaHatdl2
     &      +               2.*mu0hatSI*dThetaHatdl*dmu0hatSIdl)
     
       ! Define vNhat: the saturation function of ambient nutrient concentration
       select case(nutrient_uptake) 
       ! case 1: Classic Michaelis Menton 
         case(1)
       ! Potential maximal nutrient-limited uptake
       Kn = ScaleTrait(PMU(k),K0star, alphaK)  ! Half-saturation constant of nitrate uptake
       VNhat = V0hat*NO3(k)/(NO3(k) + Kn)
       dVNhatdl = -VNhat*alphaK*Kn/(NO3(k)+Kn) + alphaV*VNhat
       d2VNhatdl2 = -alphaK*(Kn/(NO3(k)+Kn)*dVNhatdl
     &  - Kn*VNhat*alphaK*NO3(k)/(NO3(k)+alphaK)**2)
     &  + alphaV*dVNhatdl
       ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
         case(2)
       Lmin = log(0.2**3/6*pi) + log(1+Penfac)/( Penfac*alphaA) 
       A0hat  = PenAff(PMU(k), alphaA, Penfac, Lmin)
     &        * ScaleTrait(PMU(k),dtdays*tf(self,tC)*A0star,alphaA)
       dA0hatdl = alphaA*A0hat
     &        - A0star*exp(PMU(k)*alphaA)*Penfac*alphaA
     &        * exp(Penfac*alphaA*(PMU(k)-Lmin))
       d2A0hatdl2 = alphaA*dA0hatdl
     &        - Penfac*alphaA*exp(alphaA*((1.+Penfac)*PMU(k)-Penfac*Lmin))
     &        * (dA0hatdl + A0hat*alphaA*(1.+Penfac))  
    
            !Define fA
            fA = 1.0/( 1.0 + sqrt(A0hat * NO3(k)/V0hat) ) 
         VNhat = (1.0-fA)*V0hat*fA*A0hat*NO3(k)
     &           /((1.0-fA)*V0hat + fA*A0hat*NO3(k)) 
     
        dVNhatdl= alphaV*VNhat
     &   - VNhat**2/NO3(k)/A0hat*(1.+sqrt(NO3(k)*A0hat/V0hat))
     &   * (alphaV-dA0hatdl/A0hat)  !Should be OK
        
        d2VNhatdl2 = alphaV*dVNhatdl-
     &   (VNhat**2/NO3(k)/A0hat*(1.+sqrt(NO3(k)*A0hat/V0hat))
     &   * (dA0hatdl**2/A0hat**2-d2A0hatdl2/A0hat)
     &   + VNhat**2/NO3(k)/A0hat*(alphaV-dA0hatdl/A0hat)
     &   * 0.5*sqrt(NO3(k)/A0hat/V0hat)*(dA0hatdl-A0hat*alphaV)
     &   + (1.+sqrt(NO3(k)*A0hat/V0hat))*(alphaV-dA0hatdl/A0hat)
     &   * VNhat**2/NO3(k)/A0hat*(2.*alphaV-dA0hatdl/A0hat))  !Should be OK
        end select  
       
        ! Define optimal fV
      IF( fV_opt .EQ. .TRUE. ) THEN
          ! Define ZINT
          ZINT = Qs*(muIhat/vNhat+zetaN)
          Qopt = ( 1.0 + sqrt(1.0 + 1.0/ZINT) )* Qs
        dZINdl = Qs*(dmuIhatdl/VNhat - muIhat*dVNhatdl/VNhat**2)
     &           +  alphaQ*ZINT    
     
      dQoptdl  = alphaQ*Qopt - dZINdl*Qs/(2.*ZINT*sqrt(ZINT*(1.+ZINT))) !confirmed  
     
      d2ZINdl2  = Qs/VNhat*((alphaQ-dVNhatdl/VNhat)*dmuIhatdl
     &  + d2muIhatdl2)
     &  - Qs/VNhat**2*(muIhat*d2VNhatdl2
     &  + dVNhatdl*(dmuIhatdl+alphaQ*muIhat-2.*muIhat/VNhat*dVNhatdl))
     &  + alphaQ*dZINdl
     
        d2Qoptdl2 =  alphaQ*dQoptdl
     & -  Qs/(2.*ZINT*sqrt(ZINT*(ZINT+1.)))
     & *  (d2ZINdl2+alphaQ*dZINdl
     & -  (2.*ZINT+1.5)/(ZINT*(ZINT+1.))*dZINdl**2)      ! Confirmed
      
           if (PARsup .gt. I_zero) then   
              fV = (-1.0 + sqrt(1.0 + 1.0/ZINT))*Qs*muIhat/vNhat
     !
           dfVdl = alphaQ*Qs*(1/Qopt + 2*zetaN)
     &             -  (zetaN + Qs/Qopt**2) * dQoptdl  !Confirmed
     !
           d2fVdl2 = (alphaQ**2)*Qs*(1/Qopt + 2*zetaN )
     &      -    2.*alphaQ*Qs* dQoptdl/Qopt**2
     &      +    2.*(dQoptdl**2)*Qs/Qopt**3
     &      -   (zetaN + Qs/Qopt**2.)*d2Qoptdl2  ! Confirmed
     ! 
             else if (PARsup .le. I_zero) then
              fV = 0.01
            dfVdl= errr
           d2fVdl2 = errr
           endif
     
           ELSE
             fV   = fV_initial 
             Qopt = 16.0/106.0 
             dfVdl= errr
            d2fVdl2=errr
            dQoptdl=errr
            d2Qoptdl2=errr
         ENDIF
     
        ! Net growth rate (d-1) of phytoplankton at the average size as a function of temperature, light, and nutrient   
        muNet = muIhat*(1.0-fV-Qs/Qopt) - zetaN*fV*vNhat
        write(*,*) 'muNet = ', muNet
       ! Here the derivative of muNet includes respiratory costs of both N Assim and Chl maintenance       
       dmuNetdl = muIhat*(Qs/Qopt**2*dQoptdl-alphaQ*Qs/Qopt-dfVdl)
     & +          (1.-fV-Qs/Qopt)*dmuIhatdl
     & -          zetaN*(fV*dVNhatdl+VNhat*dfVdl)

            d2muNetdl2 = Qs/Qopt**2*dQoptdl*dmuIhatdl
     &  + muIhat*dQoptdl*(alphaQ*Qs/Qopt**2-2.*Qs/Qopt**3*dQoptdl)
     &  + muIhat*Qs/Qopt**2*d2Qoptdl2
     &  - (alphaQ*(Qs/Qopt*(dmuIhatdl+alphaQ*muIhat)
     &  - muIhat*Qs/Qopt**2*dQoptdl))
     &  - (muIhat*d2fVdl2+dfVdl*dmuIhatdl)
     &  + dmuIhatdl*(Qs/Qopt**2*dQoptdl-alphaQ*Qs/Qopt-dfVdl)
     &  + (1.-fV-Qs/Qopt)*d2muIhatdl2
     &  - zetaN*(fV*d2VNhatdl2+2.*dfVdl*dVNhatdl+VNhat*d2fVdl2)  !dC
    !
           !  chl:C ratio [ g chl / mol C ] of the whole cell, at the mean cell size 
           Theta = ThetaHat*(1.0 - fV - Qs/Qopt)
           
           ! Calculate the mean chl:C ratio of phytoplankton (averaged over all sizes) 
           ! using the second derivative of the cell quota w.r.t. log size. 
           ThetaAvg = Theta + (VAR/2.)
     &    *  ( d2ThetaHatdl2*(1.0 - fV - Qs/Qopt)
     &    -      ThetaHat   *( d2fVdl2 + Qs*((dQoptdl**2)*(2./Qopt) - d2Qoptdl2)/Qopt**2 ))
    
           ! Calculate the mean N:C ratio of phytoplankton (averaged over all sizes) 
           ! using the second derivative of the cell quota w.r.t. log size. (How to derive?)
           QNavg = Qopt/(1.+((2/Qopt)*dQoptdl**2 - d2Qoptdl2 )*VAR/(2*Qopt))  
           
           ! Phytoplankton primary production
           Cpp  = muNet * PHY / QNavg
           
           ! Total Chl a (unit: mg Chl m-3)
           Chltot = PHY/QNavg*Theta


! END OF COMPUTING PHYTOPLANKTON GROWTH RATE.  

             
#ifdef DIAGNOSTICS_BIO
            ThisFlux(k, NFlux_NewProd) = cff*NO3(k)
!
# ifdef OXYGEN
!           production of O2 by phyto growth
            ThisFlux(k, OGain_NewProd) =
     &           ThisFlux(k, NFlux_NewProd) * (CN_Phyt + 2.)
# endif /* OXYGEN */
#endif /* DIAGNOSTICS_BIO */
# ifdef OXYGEN
            O2(k) = O2(k) + cff*NO3(k)*(CN_Phyt + 2.)
# endif
            ! Calculate new state variables
                par = PARsup*attn
!      Calcul of the euphotic depth
            if (PARsup .ge. I_min) then
              if (par  .ge. I_min) then
                hel(i,j)=-z_w(i,j,k-1)
              else
                hel(i,j)=-z_r(i,j,k)
              endif
            endif
!
! Start zooplankton component. Zooplankton grazing was modeled
! following Pahlow et al. MEPS (2010) Optimal current feeding model

                ! The total amount of phytoplankton grazed by zooplankton
                ! (molN;gmax is the maximal specific ingestion rate!)
                ! ghat: the ingestion function describing the relationship between 
                ! phytoplankton concentration and zooplankton ingestion rate. 
             SELECT CASE(grazing_formulation)
             
             ! Holling Type I
             case (1)
             ghat = min(PHYT(k)/2.0/kp,1.0)
             
             ! Holling Type II
             case (2)
             ghat = PHYT(k)/(kp + PHYT(k))  
             
             ! Holling Type III
             case (3) 
             ghat = PHYT(k)**2.0/(kp**2.0 + PHYT(k)**2.0)
             
             ! Ivlev
             case (4)
             ghat = 1.0 - exp(-log(2.0)*PHYT(k)/kp)  
             
             END SELECT
             ! Calculate feeding threshold of ghat
             g_th = Cf/Emax/(1.0 - Ca)  
             
             ! Calculate total foraging activity
             aAt = tf(Ea,TEMP(k))*Tgmax/beta
     &           * (-1.0 - w1(-(1.0 - g_th)/exp(1.0 + beta)))
             
             ! Compute the optimal foraging activity to maximize net growth of zooplankton
             IF (ghat .gt. g_th) THEN
              if(CF_opt .eq. .true.) then
              ! Optimal foraging
              aAf = aAt*beta/(-1.0-w1(-(1.0-g_th/ghat)/exp(1.0+beta)))
              ! Control model
              else
              aAf = aAt * 0.3
              endif
 
              AE  = Emax*(1.0 - exp(-beta*(aAt/aAf - 1.0)))   ! Assimilation efficiency
               
             ELSE
               aAf = 0.0
               AE  = Emax
             ENDIF
                
               ! Ingestion rate of zooplankton 
               INGES = aAf * ghat
               
               ! Zooplankton respiration (-> Nut)
               RES = TRMzoo*tf(Ea,TEMP(k)) + Ca*AE*INGES + Cf*aAf    
             
               ! Zooplankton net growth
               NGR = AE*INGES - RES
             
               ! ZOOPLANKTON EGESTION (-> Detritus)
               EGES = INGES - NGR - RES
             
               ! Nutrient changes. 
               dn = -PHYT(k)*mu+ZOO(k)*RES+Trdn*DET(k)*tf(Ea,TEMP(k))
             
               ! dp = phytoplankton production - zooplankton grazing - phytoplankton->detritus (molN d-1)
               dp = PHYT(k)*mu - ZOO(k)*INGES -
     &              rpd*PHYT(k)*tf(Ea,TEMP(k))
             
               ! dz = zooplankton growth(excluding excretion and egestion) - zooplankton mortality (molN d-1)
               dz = ZOO(k)*NGR - ZOO(k)**2.*Tmz*tf(Ea,TEMP(k))
             
               ! dD = phytoplankton mortality + zooplankton mortality + zooplankton egestion - conversion-to-nutrient (molN d-1)
               dd = rpd*PHYT(k)*tf(Ea,TEMP(k))   + 
     &              ZOO(k)**2.*Tmz*tf(Ea,TEMP(k)) +
     &              ZOO(k)*EGES - Trdn*DET(k)*tf(Ea,TEMP(k)) 
 
!
!        else
!#ifdef DIAGNOSTICS_BIO
!          do k = N, 1, -1
!            ThisFlux(k, NFlux_NewProd) = 0.0
!# ifdef OXYGEN
!            ThisFlux(k, OGain_NewProd) = 0.0
!# endif /* OXYGEN */
!          enddo
!#endif /* DIAGNOSTICS_BIO */
!#ifdef AVERAGES
!! if the sun is down, set the logical variable "sun_down"
!! to true for not taking into account this time step in the averaging
!          if ((ZEROTH_TILE).and.(srflx(Istr,Jstr).eq.0.)) then
!            sun_down=.true.
!          endif
!#endif
!          hel(i,j)=0.0
!        endif
!


!          cff1=dtdays*gmax*Zoo(k)/(K_Phyt+Phyt(k))
!          cff2=dtdays*mu_P_D
!          Phyt(k)=Phyt(k)/(1.+cff1+cff2)
!          Zoo(k)=Zoo(k)+Phyt(k)*cff1*beta
!#ifdef DIAGNOSTICS_BIO
!          ThisFlux(k, NFlux_Grazing)=Phyt(k)*cff1*beta
!          ThisFlux(k, NFlux_SlopFeed) = Phyt(k) * cff1 * (1.-beta)
!          ThisFlux(k, NFlux_Pmort) = Phyt(k) * cff2
!#endif /* DIAGNOSTICS_BIO */
!          Det(k)=Det(k)+Phyt(k)*(cff1*(1.-beta)+cff2)
!
! (1) Zoo excretion to NO3  (rate mu_Z_A)
! (2) Zoo mortality to Det (rate mu_Z_D)
!
!          cff1=dtdays*mu_Z_A
!          cff2=dtdays*mu_Z_D*Zoo(k)
!          Zoo(k)=Zoo(k)/(1.+cff1+cff2)
!#ifdef DIAGNOSTICS_BIO
!          ThisFlux(k, NFlux_Zmetab)=cff1*Zoo(k)
!          ThisFlux(k, NFlux_Zmort)=cff2*Zoo(k)
!# ifdef OXYGEN
!!         Zoo uptake of O2 (rate t_Zbmet + R_C)
!!         there is no control yet for assuring non-negative Oxygen
!!         values!
!          ThisFlux(k, OLoss_Zmetab) =
!     &         ThisFlux(k, NFlux_Zmetab) * CN_Z
!# endif /* OXYGEN */
!#endif /* DIAGNOSTICS_BIO */
!#ifdef OXYGEN 
!          O2(k)  = O2(k) - cff1 * Zoo(k) * CN_Z
!#endif
!          NO3(k)=NO3(k)+Zoo(k)*cff1
!          Det(k)=Det(k)+Zoo(k)*cff2
!
! (1) Det remineralization to N03
!
!          cff1=dtdays*mu_D_N
!          Det(k)=Det(k)/(1.+cff1)
!#ifdef DIAGNOSTICS_BIO
!          ThisFlux(k, NFlux_ReminD)=Det(k)*cff1
!# ifdef OXYGEN
!!         Loss of O2 in Det/Det remineralization
!          ThisFlux(k, OLoss_ReminD)=
!     &             ThisFlux(k,NFlux_ReminD) * (CN_Phyt+2.)
!# endif /* OXYGEN */
!#endif /* DIAGNOSTICS_BIO */
!#ifdef OXYGEN
!          O2(k) = O2(k) - Det(k) * cff1 * (CN_Phyt+2.)
!#endif
!          NO3(k)=NO3(k)+Det(k)*cff1
             PHY(k) = PHY(k) + dp
             NO3(k) = NO3(k) + dn
             ZOO(k) = ZOO(k) + dz 
             DET(k) = DET(k) + dd
        enddo
!
!#ifdef OXYGEN
!#  ifdef OCMIP_OXYGEN_SC
!!*********************************************************************
!!  alternative formulation (Sc will be slightly smaller up to about 35
!!  C)
!!  Computes the Schmidt number of oxygen in seawater using the
!!  formulation proposed by Keeling et al. (1998, Global Biogeochem.
!!  Cycles, 12, 141-163).  Input is temperature in deg C.
!!
!        Sc = 1638.0 - 81.83*tem(N) + 1.483*(tem(N)**2) -
!     &       0.008004*(tem(N)**3)
!!*********************************************************************
!#  else /* OCMIP_OXYGEN_SC */
!!       calculate the Schmidt number for O2 in sea water [Wanninkhof,
!!       1992]
!        Sc=1953.4 - 128.0*tem(N) + 3.9918*(tem(N)**2) -
!     &     0.050091*(tem(N)**3)
!#  endif /* OCMIP_OXYGEN_SC */
!!
!!       calculate the wind speed from the surface stress values
!        u10_loc = sqrt(sqrt( (0.5*(sustr(i,j)+sustr(i+1,j)))**2
!     &                         +(0.5*(svstr(i,j)+svstr(i,j+1)))**2)
!     &       * rho0 * 550.)   ! 550 = 1 / (1.3 * 0.0014) (=rho_air * CD)
!!       calculate the gas transfer coef for O2
!        Kv_O2_loc=0.31*u10_loc*u10_loc*sqrt(660./Sc)/(100.*3600.) 
!!  denominator: convert Kv from [cm/h] to [m/s]
!!       calculate the saturation oxygen level
!#  ifdef OCMIP_OXYGENSAT
!        O2satu_loc = o2sato(tem(N), sal(N))
!#  else /* OCMIP_OXYGENSAT */
!        call O2sato(O2(N),tem(N),sal(N),den(N),O2satu_loc,satpc,AOU)
!#  endif /* OCMIP_OXYGENSAT */
!!       air-sea flux of O2
!!       abs(z_w(i,j,N-1))==> volume of upper layer
!#  ifdef DIAGNOSTICS_BIO
!        ThisGasExcFlux(OFlux_GasExc) = Kv_O2_loc * (O2satu_loc - O2(N))
!     &                        * dtsec / ( z_w(i,j,N) - z_w(i,j,N-1) )
!!       ThisGasExcFlux is positive if ocean takes up O2 from the
!!       atmosphere
!#  endif
!        O2(N) = O2(N) + Kv_O2_loc * (O2satu_loc - O2(N))
!     &                        * dtsec / ( z_w(i,j,N) - z_w(i,j,N-1) )  
!#endif /* OXYGEN */
!
! ---------------------------------------------------------------------------
! ---------------------------------------------------------------------------
! Vertical sinking: Vertical advection algorithm based on monotonic,
! continuous conservative parabolic splines.
! 1: Put biological tracers into SB array.
        do k=1,N
!          SB(k,1)=theta(i,j,k)*Phyt(k)*CN_Phyt*12.
          SB(k,1)=PHYT(k)
          SB(k,2)=DET(k)
        enddo
        wSB(1)=TwPHYT
        wSB(2)=TwDET
!        wSB(3)=wDet

        do iB=1,nsink

! Part (i): Construct parabolic splines: compute vertical derivatives
! of the fields SB. The derivatives are located at W-points;
! Neumann boundary conditions are assumed on top and bottom.

! Neumann boundary condition: the 1st derivatives at the boundaries are constants.
! SB: biological tracers
! Hz: vertical level thickness
! FC: 
! dSB: Vertical derivatives of fields SB at the grid box interfaces

          dSB(0,iB)=0.
          FC(0)=0.
          cff6=6.
          do k=1,N-1
            cff=1./(2.*Hz(i,j,k+1)+Hz(i,j,k)*(2.-FC(k-1)))
            FC(k)=cff*Hz(i,j,k+1)
            dSB(k,iB)=cff*(cff6*(SB(k+1,iB)-SB(k,iB))
     &                             -Hz(i,j,k)*dSB(k-1,iB))
          enddo
          dSB(N,iB)=0.
          do k=N-1,1,-1     !<-- irreversible
            dSB(k,iB)=dSB(k,iB)-FC(k)*dSB(k+1,iB)
          enddo
!
! Part (ii): Convert dSB [which are now vertical derivatives
! of fields SB at the grid box interfaces] into field values
! at these interfaces, assuming parabolic profiles within each grid
! box. Restrict these values to lie between bounds determined from
! box-averaged values of grid boxes adjscent from above and below.
! (This restriction is part of PPM-like monotonization procedure.)
!
          cff=1./3.
          dSB(0,iB)=SB(1,iB) !-cff*Hz(1)*(dSB(0,iB)+0.5*dSB(1,iB))
          dSB(N,iB)=SB(N,iB) !+cff*Hz(N)*(dSB(N,iB)+0.5*dSB(N-1,iB))
          do k=2,N          !<-- irreversible
            dSB(k-1,iB)=SB(k,iB)
     &                 -cff*Hz(i,j,k)*(0.5*dSB(k,iB)+dSB(k-1,iB))
            dSB(k-1,iB)=max(dSB(k-1,iB),min(SB(k-1,iB),SB(k,iB)))
            dSB(k-1,iB)=min(dSB(k-1,iB),max(SB(k-1,iB),SB(k,iB)))
          enddo
!  
! Part (iii): Convert dSB into flux-integrated values,
! complete PPM flux limiting. This procedure starts from assigning
! Left and Right (aR,aL) values of the interpolating parabolae, then
! monotonicity conditions are checked and aL,aR are modified to fit.
! Overall, from this moment and further on it follows Colella--
! --Woodward, 1984 bombmaking code almost exactly.
!
          do k=1,N           !<-- irreversible
!            FC(k)=dtdays/Hz(i,j,k)
            FC(k)=1/Hz(i,j,k)
            aR=dSB(k,iB)
            aL=dSB(k-1,iB)
            cff1=(aR-aL)*6.*(SB(k,iB)-.5*(aR+aL))
            cff2=(aR-aL)**2
            if ((aR-SB(k,iB))*(SB(k,iB)-aL).lt.0.) then
              aL=SB(k,iB)
              aR=SB(k,iB)
            elseif (cff1.gt.cff2) then
              aL=3.*SB(k,iB)-2.*aR
            elseif (cff1.lt.-cff2) then
              aR=3.*SB(k,iB)-2.*aL
            endif
            cu=wSB(iB)*FC(k)
            dSB(k-1,iB)=SB(k,iB)-(1.-cu)*(.5*(aR-aL)-(.5*(aR+aL)
     &                                   -SB(k,iB) )*(1.-2.*cu))
          enddo
          dSB(N,iB)=0.   ! Set no-flux boundary conditions at top.
!
! Apply fluxes:
!
          do k=1,N
            SB(k,iB)=SB(k,iB)+wSB(iB)*FC(k)*(dSB(k,iB)-dSB(k-1,iB))
          enddo
        enddo  ! <-- iB

#ifdef DIAGNOSTICS_BIO
        do iflux = 1, NumVSinkTerms
          ColumnMassOld(iflux) = 0.0
          ColumnMassNew(iflux) = 0.0
        end do
#endif /* DIAGNOSTICS_BIO */

        do k=1,N
!          theta(i,j,k)= SB(k,1)/(SB(k,2)*CN_Phyt*12.+1.E-20)
!          if (theta(i,j,k).gt.theta_m) theta(i,j,k)=theta_m
#ifdef DIAGNOSTICS_BIO
! ColumnMassOld and ColumnMassNew are needed to compute the sinking flux
! into the sediment
          ColumnMassOld(1)=ColumnMassOld(1)
     &                    +PHYT(k)
          ThisVSinkFlux(k, NFlux_VSinkP1)=PHYT(k)-SB(k,2)
#endif /* DIAGNOSTICS_BIO */

          Phy(k) = SB(k,2)

#ifdef DIAGNOSTICS_BIO
          ColumnMassNew(1)=ColumnMassNew(1)
     &                    +Phyt(k)
#endif /* DIAGNOSTICS_BIO */
          ! detritus
#ifdef DIAGNOSTICS_BIO
          ColumnMassOld(2)=ColumnMassOld(2)
     &                    +Det(k)
          ThisVSinkFlux(k, NFlux_VSinkD1)=Det(k)-SB(k,3)
#endif /* DIAGNOSTICS_BIO */

          Det(k) = SB(k,3)

#ifdef DIAGNOSTICS_BIO
          ColumnMassNew(2)=ColumnMassNew(2)
     &                    +Det(k)
#endif /* DIAGNOSTICS_BIO */
        enddo
!
#ifdef DIAGNOSTICS_BIO
! Transfer fluxes to global arrays at the end of each biological time step
! for computational efficiency, divide now by dtsec to get the correct units
        do iflux = 1, NumFluxTerms
          do k = 1, N
            bioFlux(i,j,k,iflux) = ( bioFlux(i,j,k,iflux) +
     &                               ThisFlux(k, iflux) / dt )
# ifdef MASKING
     &             * rmask(i,j)
# endif /* MASKING */
          end do
        end do
        do iflux = 1, NumVSinkTerms
! Compute the vertical sinking flux into the sediment by comparing
! previous and current mass in this (i,j) column
! The flux is positive if upward, so usually it will be
! negative, i.e. into the sediment.
          LastVSinkFlux = ( ColumnMassNew(iflux) -
     &                      ColumnMassOld(iflux) )
          bioVSink(i,j,0,iflux) = (  bioVSink(i,j,0,iflux) +
     &                               LastVSinkFlux / dt )
# ifdef MASKING
     &          * rmask(i,j)
# endif /* MASKING */
          do k = 1, N
            LastVSinkFlux = LastVSinkFlux +
     &                      ThisVSinkFlux(k,iflux)
            bioVSink(i,j,k,iflux) =(  bioVSink(i,j,k,iflux) +
     &                                LastVSinkFlux / dt )
# ifdef MASKING
     &          * rmask(i,j)
# endif /* MASKING */
          end do
        end do
# ifdef OXYGEN
! ThisGasExcFlux is already in the correct units [mmol s^-1]
        do iflux = 1, NumGasExcTerms
           GasExcFlux(i,j,iflux) = ( GasExcFlux(i,j,iflux) +
     &                               ThisGasExcFlux(iflux) / dt )
#  ifdef MASKING
     &          * rmask(i,j)
#  endif /* MASKING */
        end do
# endif
#endif /* DIAGNOSTICS_BIO */
!
      ENDDO  ! <-- ITER
!
! Write back
!
!          print*,'N=',N
#undef DEBUG_BIO
# if defined DIAGNOSTICS_BIO && defined DEBUG_BIO
!
      k=N
      if ((i.eq.13).and.(j.eq.15)) then
        bilan_no3 = bioFlux(i,j,k,NFlux_Zmetab)
     &           + bioFlux(i,j,k,NFlux_ReminD)
     &           - bioFlux(i,j,k,NFlux_NewProd)
!
        bilan_phy = bioFlux(i,j,k,NFlux_NewProd)
     &           - bioFlux(i,j,k,NFlux_Pmort)
     &           - bioFlux(i,j,k,NFlux_Grazing)
     &           - bioFlux(i,j,k,NFlux_SlopFeed)
     &           - ( bioVSink(i,j,k,NFlux_VSinkP1) 
     &               - bioVSink(i,j,k-1,NFlux_VSinkP1) )
!
        bilan_zoo = bioFlux(i,j,k,NFlux_Grazing) 
     &            - bioFlux(i,j,k,NFlux_Zmetab)
     &            - bioFlux(i,j,k,NFlux_Zmort)
!
        bilan_det = bioFlux(i,j,k,NFlux_SlopFeed)
     &            + bioFlux(i,j,k,NFlux_Zmort)
     &            + bioFlux(i,j,k,NFlux_Pmort)
     &            - bioFlux(i,j,k,NFlux_ReminD)
     &            - ( bioVSink(i,j,k,NFlux_VSinkD1) 
     &                - bioVSink(i,j,k-1,NFlux_VSinkD1) )
#   ifdef OXYGEN
        bilan_o2 = bioFlux(i,j,k, OGain_NewProd)
     &           - bioFlux(i,j,k, OLoss_Zmetab)
     &           - bioFlux(i,j,k, OLoss_ReminD)
        if (k.eq.N) then
           bilan_o2 =  bilan_o2 + GasExcFlux(i,j, OFlux_GasExc)
        endif                     
        trend_o2 = ( (min(t(i,j,k,nnew,iO2),0.) +O2(k))
     &               -   t(i,j,k,nnew,iO2) ) / dt
#   endif
        somme = bilan_no3+bilan_phy+bilan_zoo+bilan_det
        trend_no3 = ( (min(t(i,j,k,nnew,iNO3_),0.) +NO3(k))
     &               - t(i,j,k,nnew,iNO3_) ) / dt
        trend_phy = ( (min(t(i,j,k,nnew,iPhy1),0.) +Phyt(k))
     &               - t(i,j,k,nnew,iPhy1) )/ dt
        trend_zoo = ( (min(t(i,j,k,nnew,izoo1),0.) +Zoo(k))
     &               - t(i,j,k,nnew,izoo1) )/ dt
        trend_det = ( (min(t(i,j,k,nnew,idet1),0.) +Det(k))
     &       - t(i,j,k,nnew,idet1) )/ dt
!
        trend_total =  trend_no3 + trend_phy + trend_zoo + trend_det
!
        sinking_loss = - ( bioVSink(i,j,k,NFlux_VSinkD1) 
     &                     - bioVSink(i,j,k-1,NFlux_VSinkD1) )
     &                 - ( bioVSink(i,j,k,NFlux_VSinkP1) 
     &                     - bioVSink(i,j,k-1,NFlux_VSinkP1) )
! 
        print*, '=================='
        print*, 'i=',i,' j=',j,' k=',k
        print*, 'somme SMS(of N)=         ',somme
        print*, 'trend_total of N=        ',trend_total
        print*, 'Sinking_loss of N=       ',sinking_loss
        print*, 'These three values have to be the same!'
        print*, 'Error for N = trend_total-somme=',trend_total-somme
        print*, "-----------------"
        print*, 'bilan_no3 - trend_no3=', bilan_no3-trend_no3
        print*, 'bilan_phy-trend_phy=  ', bilan_phy-trend_phy
        print*, 'bilan_zoo-trend_zoo=  ', bilan_zoo-trend_zoo
        print*, 'bilan_det-trend_det=  ', bilan_det-trend_det
        print*, '-----------------'
        print*, 'bioFlux(i,j,k,NFlux_NewProd)= ', 
     &                               bioFlux(i,j,k,NFlux_NewProd)
        print*, 'bioFlux(i,j,k,NFlux_Grazing)= ', 
     &                               bioFlux(i,j,k,NFlux_Grazing)
        print*, 'bioFlux(i,j,k,NFlux_SlopFeed)=', 
     &                               bioFlux(i,j,k,NFlux_SlopFeed)
        print*, 'bioFlux(i,j,k,NFlux_Pmort)=   ', 
     &                               bioFlux(i,j,k,NFlux_Pmort)
        print*, 'bioFlux(i,j,k,NFlux_Zmetab)=  ', 
     &                               bioFlux(i,j,k,NFlux_Zmetab)
        print*, 'bioFlux(i,j,k,NFlux_Zmort)=   ', 
     &                               bioFlux(i,j,k,NFlux_Zmort)
        print*, 'bioFlux(i,j,k,NFlux_ReminD)=  ', 
     &                               bioFlux(i,j,k,NFlux_ReminD)
        print*, 'bioVSink(i,j,k,NFlux_VSinkP1)=', 
     &                                bioVSink(i,j,k,NFlux_VSinkP1)
        print*, 'bioVSink(i,j,k,NFlux_VSinkD1)=', 
     &                                bioVSink(i,j,k,NFlux_VSinkD1)
#    ifdef OXYGEN
        print*, '-----------------'
        print*, 'Error for O2 = trend_O2-SMS(O2)=',trend_o2-bilan_o2
        print*, 'bioFlux(i,j,k,OGain_NewProd)=',
     &                               bioFlux(i,j,k,OGain_NewProd)
        print*, 'bioFlux(i,j,k,OLoss_Zmetab)= ',
     &                               bioFlux(i,j,k,OLoss_Zmetab)
        print*, 'bioFlux(i,j,k,OLoss_ReminD)= ',
     &                               bioFlux(i,j,k,OLoss_ReminD)
        if (k.eq.N) then
        print*, 'GasExcFlux(i,j,OFlux_GasExc)= ',
     &                                GasExcFlux(i,j,OFlux_GasExc)
        endif 
#    endif
        print*, '=================='
      endif
# endif

          do k=1,N
            t(i,j,k,nnew,iNO3_)=min(t(i,j,k,nnew,ino3_),0.) +NO3(k)
            t(i,j,k,nnew,iPhy1)=min(t(i,j,k,nnew,iphy1),0.) +PHYT(k)
            t(i,j,k,nnew,iZoo1)=min(t(i,j,k,nnew,izoo1),0.) +ZOO(k)
            t(i,j,k,nnew,iDet1)=min(t(i,j,k,nnew,idet1),0.) +DET(k)
!            t(i,j,k,nnew,iChla)=min(t(i,j,k,nnew,iChla),0.) +
!     &                              CN_Phyt*12.*Phyt(k)*theta(i,j,k)
# ifdef OXYGEN
            t(i,j,k,nnew,iO2)  =min(t(i,j,k,nnew,iO2),0.)   +O2(k)
# endif
          enddo
# ifdef OXYGEN
          O2satu(i,j) = O2satu_loc
          Kv_O2(i,j) = Kv_O2_loc
          u10(i,j) = u10_loc
!         if (i==10.and.j==10) then 
!            write (*,*) 'O2satu(10,10)=',O2satu(10,10)
!            write (*,*) 'Kv_O2(10,10)=',KV_O2(10,10)
!            write (*,*) 'u10(10,10)=',u10(10,10)
!          endif
# endif /* OXYGEN */
        enddo
      enddo


#else
      subroutine biology_empty ()
#endif
      return
      end

!---------------------------------------------------------
      pure real function tf(Ea, tC)
      implicit none
! !IROUTINE: Temperature dependence of plankton rates. Added by Bingzhang Chen on 20150213
! !DESCRIPTION:
! Here, the temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
! !INPUT PARAMETERS:
      real, intent(in)         :: tC, Ea
      real, parameter          :: kb = 8.62E-5 ! boltzman constant constant [ eV /K ]
      real, parameter          :: Tr = 15.0    ! reference temperature [ ºC ]     
      tf=exp(-(Ea/kb)*(1.0/(273.15+tC) - 1.0/(273.15+Tr)))

      end function tf
!-----------------------------------------------------------------------
      ELEMENTAL real FUNCTION w1(x)
      IMPLICIT NONE
      REAL, PARAMETER  :: m1=0.3361, m2=-0.0042, m3=-0.0201
      REAL, INTENT(IN) :: x
      REAL :: s
      s = -1.0 - LOG(-x)
      w1 = -1.0 - s - 2.0/m1*(1.0 - 1.0/(1.0 + (m1*SQRT(0.5*s)
     &    /(1.0 + m2*s*EXP(m3*SQRT(s))))))
      END FUNCTION w1
!-----------------------------------------------------------------------    
       real function PenAff( logsize, alpha, Pfac, lmin ) 
         implicit none
         real, intent(IN) :: logsize, alpha, Pfac, lmin 
     
         ! A 'penalty' function to reduce the value of affinity for nutrient at very small cell sizes
         ! in order to avoid modeling unrealistically small cell sizes.  This is needed because affnity
         ! increases with decreasing cell size, which means that under low-nutrient conditions, without
         ! such a penalty, unrealistically small cell sizes could be predicted.
         ! This penalty function becomes zero at logsize = lmin.   
      
         PenAff = 1.0 - exp(Pfac*alpha*(logsize - lmin))
     
         return
       end function PenAff
!-----------------------------------------------------------------------    
      subroutine OptCurrFeed(self,Ptot,tC,AE,INGES,RES,NGR,EGES)
! This subroutine calculates the ingestion, respiration, egestion, and net growth rate
! of zooplankton given the total amount of prey available 
       use lambert
       implicit none
       type (type_gotm_PZEFT), intent(in) :: self
       real(rk), intent(IN)  :: Ptot,tC
       real(rk), intent(OUT) :: AE,INGES,RES,NGR,EGES
       ! Declare local variables
       real(rk) :: ghat,g_th,aAt,aAf 
      
       if (Ptot .le. 0.0) then
        write(*,*) 'error: Total phytoplankton biomass being negative!'
        stop 
       endif
     ! Calculate the feeding probability function based on Ptot
     ! kp relates to the capture coefficient
       SELECT CASE(self%grazing_formulation)
       
       ! Holling Type I
       case (1)
       ghat = min(Ptot/2.0/self%kp,1.0)
       
       ! Holling Type II
       case (2)
       ghat = Ptot/(self%kp + Ptot)  

       ! Holling Type III
       case (3) 
       ghat = Ptot**2.0/(self%kp**2.0 + Ptot**2.0)

       ! Ivlev
       case (4)
       ghat = 1.0 - exp(-LOG(2.0)*Ptot/self%kp)  
      
       END SELECT
     ! Calculate feeding threshold of ghat (dimensionless)
       g_th = self%Cf/self%Emax/(1.0 - self%Ca)  

     ! Calculate total foraging activity (d-1)
       aAt = tf(self,tC)*self%gmax/self%beta*      &
     &      (-1.0 - w1(-(1.0-g_th)/EXP(1.0+self%beta)))

     ! Compute the optimal foraging activity to maximize net growth of zooplankton
       IF (ghat .GT. g_th) THEN
         if (self%CF_opt .eq. .true.) then
        ! Optimal foraging, aAf: Specific foraging activity (d-1) = 1/h (h: biomass specific handling time)
        !With increasing food conc., optimal aAf slightly increases with food conc.
           aAf = aAt*self%beta/(-1.0 - w1(-(1.0 - g_th/ghat) &
     &         /EXP(1.0 + self%beta)))
           ! Control model
           else
           aAf = aAt * 0.3
         endif 
      ! Assimilation efficiency relates with two factors
      ! 1) The difference between total activity (aAt) and foraging activity (aAf)
      ! 2) Gut passage time
      ! AE is negatively related with foraging activity
      ! because digestion enzymes decrease with increasing food
      ! Another factor, gut passange time linearly increases with handling time
         AE  = self%Emax*(1.0 - exp(-self%beta*(aAt/aAf - 1.0)))  
          
         ELSE

         aAf = 0.0
         AE  = self%Emax
       ENDIF
     
       ! Ingestion rate of zooplankton 
       INGES = aAf * ghat
       
       ! Zooplankton respiration (-> Nut)
       RES = self%RMzoo*tf(self,tC) + self%Ca*AE*INGES + self%Cf*aAf    

       ! Zooplankton net growth
       NGR = AE*INGES - RES

       ! ZOOPLANKTON EGESTION (-> Detritus)
       EGES = INGES - NGR - RES
 
  end subroutine OptCurrFeed

