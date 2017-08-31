! THis is a module based on an Optimality-based NPZD model that can couple with
! other physical (0D,1D,3D) hosts.
! The optimal phytoplankton growth part based on Pahlow et al. MEPS 2013.
! The optimal zooplankton grazing part based on Pahlow & Prowe MEPS 2010.
! Original structure taken from GOTM and adapted for FABM by Jorn Bruggeman. 
! !Add the component of optimal current feeding on 20150309 
! Modify the unimodal size functional relationship of growth rate in July, 2015 
! Add 5 submodels using cpp predirectives: 1) NPZDSIMPLE (the simplest NPZD)
! 2) NPZDCHL (NPZD with dynamic Chl); 3) EFTNOSIZE (adaptive model but no size)
! 4) EFTSIZE (adaptive model with size); 5)NPZDSIZE (non-adaptive model with size)
#include "cppdefs.h"
#if defined BIOLOGY && defined BIO_PZEFT
      MODULE PZEFT
      use lambert
      implicit none
      private
      integer,parameter :: NENV = 5,
#ifdef NPZDSIMPLE
     & NPAR=14, NVAR=4
#elif defined EFTNOSIZE
     & NPAR=18, NVAR=4
#elif defined NPZDSIZE
     & NPAR=24, NVAR=6
#elif defined EFTSIZE
     & NVAR=6,  NPAR=
#  ifdef OCF
     & 34
#  else
     & 29
#  endif 
#endif

! PUBLIC DERIVED TYPES:
      type,public :: type_PZEFT   ! Define a type
!     Variable identifiers
! Input variables: (These variables must be provided in the host model)
!-------------------------------------------------------------------------------
! Modify input data and parameters as vectors (to facilitate MCMC fitting)
      real    :: env(NENV),vars(NVAR)
!-------------------------------------------------------------------------------
! The following parameters are internal and should be parameterized in the 'initialize' subroutine
      real    :: pars(NPAR),PMU,VAR
      integer :: nutrient_uptake,grazing_formulation,zoo_mortality
#if defined EFTNOSIZE || defined EFTSIZE
      real :: fV_initial
#endif
! Output variables (to be used by the host model): 
!-------------------------------------------------------------------------------
!      real :: deriv(NDERIV,2) !1st and 2nd derivatives
      real :: ThetaAvg,QNavg, muNet,w_pAvg
#if defined EFTSIZE || defined EFTNOSIZE
      real :: fV
#endif
#if defined NPZDSIZE || defined EFTSIZE 
      real :: dmuNetdl,d2muNetdl2,gbar,dgdlbar,d2gdl2bar,dwdl,d2wdl2
#endif
      real :: AE=0.,INGES=0.,RES=0.,NGR=0.,EGES=0.
      real :: odeout(NVAR)
      real :: pp(4,4), dd(4,4),  cc(4)      ! For PPDD
!-------------------------------------------------------------------------------
      contains
         procedure :: initialize
         procedure :: tf=>temp
         procedure :: phyto
         procedure :: Zooplankton
         procedure :: ODE
         procedure :: modified_patankar
         procedure :: RK2
         procedure :: modified_patankar_2
      end type type_PZEFT
      integer, parameter, public :: iNO3 =1,  iPHY = 2, iZOO = 3, 
     &       iDET =4
#if defined NPZDSIZE || defined EFTSIZE
     &       ,iPMUPHY =5,!iPMUPHY is the product of PHY and PMU (in biovolume)
                       !that will enter diffusion and advection 
                       !in the 3D physical model
     &       iVARPHY =6!iVARPHY is the product of PHY and cell variance
#endif

      integer, parameter, public ::  ! For environmental variables
     &        iCHL =1,  ipar = 2,
     &        itC  =3,  iI_0 = 4, 
     &      idtdays=5

      integer, parameter, public :: imu0 =1, iaI0 = 2, iK0 = 3,
     & iV0   = 4, iw_p = 5, iw_d = 6,  ikc = 7,      
     & igmax = 8, ikp  = 9, irdn =10,
     & irpd  =11,
     & iEa   =12, imz  =13, ikzb =14
#if defined NPZDSIZE || defined EFTSIZE
     & ,ialphamu = 15,
     & ibetamu  = 16,
     & ialphaI  = 17,
     & ialphaW  = 18,
     &      iA0 = 22,
     &  ialphaA = 23,
     & ialphaK  = 19,
     & ialphaV  = 20,
     & ialphaG  = 21,
     &      iS  = 24    !Size based grazing selectivity
#endif
#ifdef EFTNOSIZE
     & ,  iQ0   = 15,
     & iRMchl   = 16,
     & izetaChl = 17,
     & izetaN   = 18
#endif 
#if defined EFTSIZE
     & ,iQ0     =27, 
     &  ialphaQ =28,
     &  iRMchl  =29,
     &  izetaChl=25,
     &  izetaN  =26
#  ifdef OCF
     &,iCa      =32, iCf     =33,
     & iEmax    =34, iRMzoo  =30,
     & ibeta    =31
#  endif /* OCF */
#endif
      integer, parameter :: stdout= 6
      real,    parameter :: eps  = 1d-5, pi=3.141592653529D0,
     &                      PMUmin=-2.,VARmin=0.01 
!-----------------------------------------------------------------------
   
      contains
   
!-----------------------------------------------------------------------
!
! !IROUTINE: Initialise the PZEFT model parameters
! All the values were obtained from adaptive MH fitting to global Chl:C and growth rate data
! Currently, only the optimality based model also with optimal nutrient uptake were fitted. 
! !INTERFACE:
      subroutine initialize(self)
! !INPUT PARAMETERS:
      class (type_PZEFT), intent(inout), target :: self 
      real :: dtdays,temp
      real :: w31_p, w31_d
      integer :: i

      w31_p  = 0.1   !Stoke's law: V = 2/9*(rho_p - rho_f)*g*R**2/mu (mu: dynamic viscosity, which is function
                     !of temperature.) w_p0 = w(31)*(0.451 + 0.0178*T) (Gregg et al. 2008). 
                     !We also set this to be size dependent
                     !such that w_p = w_p0* Wstar * exp(logsize * alphaW)
      w31_d  = 5.0 
      dtdays = self%env(idtdays) 
      temp   = self%env(itC    )
      self%pars(ikc    )    = 0.025      !Light attenuation of chlorophyll (mg Chl m-2)
      self%pars(imu0   )    = 
#if defined NPZDSIMPLE || defined NPZDSIZE
     &         1.5D0*dtdays
#else
     &         5D0*dtdays      !Smith et al. 2015
#endif
! Note:
! Initial slope of the P-I curve (mol C gChl-1 (W m-2)-1 d-1). 
! 1 W m-2 = 4.6 uE m-2 s-1 = 1.84 E m-2 d-1
! Pahlow et al. (2013) estimated the aI0 0.5~3.7 m2 mol C gChl-1 E-1
! corresponding to 0.27~2.01 mol C gChl-1 (W m-2)-1 d-1.
! Ward et al. (L&O) used the size independent slope of 3.83E-7 mol C (gChl-1) (uE m-2)-1. 
! Edwards et al. L&O (2015) estimated a size dependent value of 10**(-1.36) d-1 (uE m-2 s-1)-1 (= 5.05E-7 (uE m-2)-1)
! So I choose a value of 10**(-1.36)*4.6/(50/12)= 0.05 
      self%pars(iaI0   )    = 0.25*dtdays   ! Fennel et al. (2006)
      self%pars(iV0    )    = 5d0*dtdays      ! Smith et al. 2015
      self%pars(iK0    )    = 
#if defined NPZDSIMPLE || defined NPZDSIZE
     &           0.02    !Half saturation constant for growth
#else
     &           0.17    !Half saturation constant for nitrogen uptake (Litchman et al. 2007)
#endif
!      self%pars(iw_p   )    = w31_p*dtdays*(0.451+0.0178*temp) !Temperature dependent sinking rate of phytoplankton, POSITIVE!!
      self%pars(iw_p   )    = w31_p*dtdays
      self%pars(iw_d   )    = w31_d*dtdays
!      self%pars(iw_d   )    = w31_d*dtdays*(0.451+0.0178*temp)!Sinking rate of detritus
      self%pars(igmax  )    = 
#ifdef OCF
     &  2.5D0*dtdays !Maximal biomass-specific ingestion rate of zooplankton at 15 ºC
#else 
     &  1D0*dtdays   !Fennel et al. (2006)
#endif
      self%pars(ikp    )    = 1D0   !Fennel et al. (2006),half-saturation constant for zooplankton grazing
      self%pars(irdn   )    = 0.05*dtdays!The conversion rate from detritus to nutrient at reference temperature (Fennel et al. 2006)
      self%pars(irpd   )    = 0.02*dtdays!Phytoplankton mortality rate at reference temperature in the euphotic zone (Ward et al. 2012)
      self%pars(iEa    )    = 0.65             !Activation energy (eV)
      self%pars(imz    )    = 0.025*dtdays   !Fennel et al. 2006. Coefficient for zooplankton mortality at 15 C
      self%pars(ikzb   )    = 0.5            !The mortality term b for zooplankton 
#if defined NPZDSIZE || defined EFTSIZE
      self%pars(ialphamu)   = 1.2d-1
      self%pars(ibetamu )   = -0.02
      self%pars(ialphaI)    = 0.0     ! Smith et al. 2015,but -0.13 Edwards et al. L&O 2015
      self%pars(iA0    )    = 2D1*dtdays     ! From Pahlow et al. MEPS 2013, not use it now.    
      self%pars(ialphaA)    = -0.289
      self%pars(ialphaV)    = -0.27

!Normalization constant for half-saturation constant for nitrate (mmolN m-3), from Litchman et al. 2007
      self%pars(ialphaK)    = 0.27    
      self%pars(ialphaG)    = 1.1       !Not estimated
      self%pars(ialphaW)    = 0.39       ! Laws 1975
      self%pars(    iS )    = -0.25     ! Acevedo-Trejos et al. 2015
#endif
#if defined EFTSIZE || defined EFTNOSIZE
      self%fV_initial = 0.25 !Initial value for the fraction of energy devoted to nutrient uptake
      self%pars(iQ0    )    = 0.07    ! 0.07 in Litchman et al. 2007; Ward et al. L&O 2012 
#  ifdef EFTSIZE
      self%pars(ialphaQ)    = -0.17   ! -0.17 in Litchman et al. 2007; Ward et al. L&O 2012 
#  endif
      self%pars(izetaN  )   = 0.6     ! Smith et al. 2015;Cost of DIN uptake (mol C/mol N)   
      self%pars(iRMchl  )   = 0.1*dtdays    !Cost of chl maintenance (d-1) 
      self%pars(izetaChl)   = 0.8         !Cost of photosynthesis coefficient (mol C/g Chl)
#endif
#  ifdef OCF 
      self%pars(ibeta   )   = 0.2             !Assimilation coefficient 
      self%pars(iCa     )   = 0.3             !Cost of assimilation coefficient
      self%pars(iCf     )   = 0.3             !Effective cost of foraging coefficient  
      self%pars(iEmax   )   = 0.99            !Maximal assimilation efficiency
      self%pars(iRMzoo  )   = 0.16*dtdays     !Specific maintenance respiration (d-1)
#  endif /* OCF */
! End of pars definition.

! Choose Different models:
      self%nutrient_uptake=1      !The model used for nutrient uptake
      self%grazing_formulation=3
      self%zoo_mortality=2


      do i = 1,NVAR
         self%vars(i) = max(self%vars(i),eps)   !To avoid negative values
      enddo

! Calculate real PMU and VAR:
#if defined NPZDSIZE || defined EFTSIZE
      if (self%vars(iPHY) .le. eps) then
        self%PMU = PMUmin
        self%VAR = VARmin   !To avoid very large PMU and VAR
      else
        self%PMU = log(self%vars(iPMUPHY)/self%vars(iPHY)) ! Convert the unit to log(mu**3)
        self%VAR =     self%vars(iVARPHY)/self%vars(iPHY)
      endif
#endif 

      end subroutine initialize
   
!---------------------------------------------------------
      function temp(self) result(tf)
!IROUTINE: Temperature dependence of plankton rates. Added by Bingzhang Chen on 20150213
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
        class (type_PZEFT), intent(in), target :: self
        real  :: tf,Ea,tC
        real, parameter:: kb = 8.62d-5, ! boltzman constant constant [ eV /K ]
     &                    Tr = 15d0
       
       Ea=self%pars(iEa)
       tC=self%env(itC)
       tf=exp(-(Ea/kb)*(1D0/(273.15 + tC)-1D0/(273.15 + Tr)))
     
       end function temp
!EOC
 
!-----------------------------------------------------------------------
! !IROUTINE:Computing phytoplankton growth components  
!
! !INTERFACE:
      subroutine phyto(self)
!INPUT PARAMETERS:
      class (type_PZEFT),intent(inout) :: self
!LOCAL VARIABLES:
      real :: I_zero, tf, PMU,VAR,par,NO3,zetaChl,alphaI,alphamu,
     &        betamu, PHY
      real :: dmu0hatdl, d2mu0hatdl2
      real :: V0hat,Kn,Lmin,A0hat,dA0hatdl,d2A0hatdl2,fA
      real :: VNhat,dVNhatdl,d2VNhatdl2, fN, d2fNdl2 ! Nutrient uptake variables
      real :: mu0hat,muIhat,mu0hatSI,dmu0hatSIdl,d2mu0hatSIdl2
      real :: dmuIhatdl,d2muIhatdl2! Growth rate variables
      real :: fV,dfVdl,d2fVdl2  !fV
      real :: ZINT,dZINdl,d2ZINdl2 !ZINT
      real :: aI,SI,dSIdl,d2SIdl2         !Light dependent component
      real :: RMchl,Theta,ThetaHat,dThetaHatdl,d2ThetaHatdl2 !Chl related variables
      real :: Qopt,Qs,dQoptdl,d2Qoptdl2  ! cell quota related variables
      real :: larg1,larg,w_p,dmu0hat_aIdl,d2mu0hat_aIdl2,dlargdl,
     &        d2largdl2,W_Y,dWYYdl,daI_mu0hatdl,d2aI_mu0hatdl2,
     &        alphaK,alphaV,alphaA,alphaQ,zetaN,d2wpdl2
      real, parameter :: Penfac=3D0
   
!-----------------------------------------------------------------------
! Phytoplankton codes below mostly come from S. Lan Smith
! calculate the trait values at the mean size, lmean 
! multiplying RATES only by tC (Temperature dependence)

! Only for pot. max. affinity, apply a penalty function (to avoid unrealistic increase at v. small cell size). 
! To get the max. affinity at ESD = 0.2 micrometers: 
! The maximal growth rate (mu0hat) under temperature tC 

      tf     = self%tf()
      par    = self%env(ipar)
      NO3    = self%vars(iNO3)
      PHY    = self%vars(iPHY)
#ifdef NPZDSIMPLE
      mu0hat = tf*self%pars(imu0)     
      ThetaHat = 0.6
      aI     = self%pars(iaI0)
      !The light limitation index (SI)
      SI    = 1D0-exp(-aI*par*ThetaHat/mu0hat)
      ! Half-saturation constant of nitrate uptake
      Kn = self%pars(iK0)
      self%muNet = mu0hat*SI*NO3/(NO3 + Kn)
      self%ThetaAvg = 0.24
      self%QNavg    = 16D0/106D0
! Total Chl a (unit: mg Chl m-3)
      self%env(iCHL) = self%vars(iPHY)/self%QNavg*self%ThetaAvg
      self%w_pAvg    = self%pars(iw_p) ! Sinking rate of phytoplankton
#endif

#ifdef EFTNOSIZE
      mu0hat = tf*self%pars(imu0)     
      aI     = self%pars(iaI0)
      zetaChl= self%pars(izetaChl)
      V0hat  = tf*self%pars(iV0)
      zetaChl= self%pars(izetaChl)
      zetaN  = self%pars(izetaN )
      aI     = self%pars(iaI0)
      Qs     = self%pars(iQ0)/2D0
      RMchl  = tf*self%pars(iRMchl)
      ! Threshold irradiance and RMchl is set temperature dependent
      I_zero = zetaChl*RMchl/aI  
      ! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)

      if( par .gt. I_zero ) then
      
         larg1 = exp(1.0 + aI*par/(mu0hat*zetaChl))
         larg  = (1.0 + RMchl/mu0hat)*larg1 
         W_Y      = WAPR(larg,0,0)
         ThetaHat = 1.0/zetaChl + (1.0- W_Y) * mu0hat/(aI * par)
      else
         ThetaHat = 0.01  !  a small positive value 
      endif
      SI       = 1D0-exp(-aI*par*ThetaHat/mu0hat)
      mu0hatSI = mu0hat*SI  ! Gross specific carbon uptake (photosynthesis)
      muIhat   = max((mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat),
     &            errr) ! Net specific carbon uptake
      Kn    = self%pars(iK0)
      VNhat = V0hat*NO3/(NO3 + Kn)
      ZINT = Qs*(muIhat/VNhat + zetaN)
      Qopt = (1D0 + sqrt(1D0 + 1D0/ZINT))*Qs
      if (par .gt. I_zero) then   
         fV = (-1.0 + sqrt(1.0 + 1.0/ZINT))*Qs*muIhat/VNhat
      else
         fV = 1D-2
      endif
      self%muNet = muIhat*(1D0-fV-Qs/Qopt)- zetaN*fV*VNhat
      self%ThetaAvg = ThetaHat*(1D0 - fV - Qs/Qopt)
      self%QNavg    = Qopt
! Total Chl a (unit: mg Chl m-3)
      self%env(iCHL) = self%vars(iPHY)/self%QNavg*self%ThetaAvg
      self%w_pAvg    = self%pars(iw_p) ! Sinking rate of phytoplankton
#endif

#ifdef NPZDSIZE
      PMU    = self%PMU
      VAR    = self%VAR
      ThetaHat = 0.6
      alphaK = self%pars(ialphaK) 
      alphamu= self%pars(ialphamu)
      betamu = self%pars(ibetamu )
      alphaI = self%pars(ialphaI )
      mu0hat = tf*self%pars(imu0)
     &   *exp(alphamu * PMU + betamu*PMU**2)

      dmu0hatdl = mu0hat*(alphamu + 2D0 * betamu * PMU)

      d2mu0hatdl2=mu0hat * 2D0*betamu
     & + mu0hat*(alphamu + 2D0*betamu*PMU)**2

      ! Initial slope of P-I curve
      aI     = ScaleTrait(PMU, self%pars(iaI0), alphaI)

          !The light limitation index (SI)
      SI    = 1D0-exp(-aI*par*ThetaHat/mu0hat)

      daI_mu0hatdl = aI*(alphaI - alphamu - 2D0*betamu*PMU)/mu0hat

      dSIdl = exp(-aI*par*ThetaHat/mu0hat)*ThetaHat*par*daI_mu0hatdl
 
      mu0hatSI = mu0hat*SI

      dmu0hatSIdl = mu0hatSI*(alphamu+2.*betamu*PMU) + mu0hat*dSIdl

      dmu0hat_aIdl   = (dmu0hatdl - alphaI*mu0hat)/aI

        d2mu0hat_aIdl2 = mu0hat/aI*2.*betamu 
     &  + (alphamu-alphaI+2.*betamu*PMU)* dmu0hat_aIdl

        d2aI_mu0hatdl2 = -2.*betamu*aI/mu0hat
     &  - (2.*betamu*PMU + (alphamu-alphaI))*daI_mu0hatdl
  
       d2SIdl2 = ThetaHat*par*(-dSIdl*daI_mu0hatdl
     &  + (1D0-SI)*d2aI_mu0hatdl2)

       d2mu0hatSIdl2 = (alphamu+2D0*betamu*PMU)*dmu0hatSIdl
     &  + 2D0*betamu*mu0hatSI + mu0hat*(alphamu+2D0*betamu*PMU)*dSIdl
     &  + mu0hat*d2SIdl2


      ! Half saturation constant for growth:
         Kn  = ScaleTrait(PMU, self%pars(iK0),alphaK)
      ! Phytoplankton growth rate at the mean size:
      self%muNet = mu0hat*SI*NO3/(NO3 + Kn)

      self%dmuNetdl = dmu0hatSIdl*NO3/(NO3+Kn) 
     &  - mu0hatSI*alphaK*Kn*NO3/(NO3 + Kn)**2

      fN = NO3/(NO3 + Kn)

      d2fNdl2 = -alphaK**2*NO3*Kn
     &  * (1D0/(NO3 + Kn)**2 - 2D0*Kn/(NO3+Kn)**3)

      self%d2muNetdl2 = 2D0*dmu0hatSIdl*dfNdl + d2mu0hatSIdl2*fN
     &  + mu0hatSI*d2fNdl2 

      self%ThetaAvg = 0.24
      self%QNavg    = 16D0/106D0
! Total Chl a (unit: mg Chl m-3)
      self%env(iCHL) = self%vars(iPHY)/self%QNavg*self%ThetaAvg
! Calculate the community sinking rate of phytoplankton
      ! Phytoplankton sinking rate at the average size
      w_p = ScaleTrait(PMU,self%pars(iw_p),self%pars(ialphaW))
      d2wpdl2 = self%pars(ialphaW)**2*w_p
      self%w_pAvg = w_p+0.5*VAR*d2wpdl2 ! Sinking rate of phytoplankton

#endif
#ifdef EFTSIZE
      alphaQ = self%pars(ialphaQ)
      alphaI = self%pars(ialphaI )
      alphamu= self%pars(ialphamu)
      betamu = self%pars(ibetamu )
      
      PMU=self%PMU
      VAR=self%VAR

      Qs     = ScaleTrait(PMU, self%pars(iQ0), alphaQ)/2D0
      mu0hat = tf*self%pars(imu0)*exp(alphamu*PMU + betamu*PMU**2)

      dmu0hatdl = mu0hat*(alphamu + 2D0*betamu*PMU)

      d2mu0hatdl2=mu0hat * 2D0*betamu
     & + mu0hat*(alphamu + 2D0*betamu*PMU)**2

      zetaChl= self%pars(izetaChl)
      zetaN  = self%pars(izetaN )

      V0hat  = ScaleTrait(PMU, tf*self%pars(iV0), self%pars(ialphaV))
      ! Initial slope of P-I curve
      aI     = ScaleTrait(PMU, self%pars(iaI0), alphaI)
      RMchl  = tf*self%pars(iRMchl)
      ! Threshold irradiance and RMchl is set temperature dependent
      I_zero = zetaChl*RMchl/aI  
! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)

      if( par .gt. I_zero ) then
        
        larg1 = exp(1D0 + aI*par/(mu0hat*zetaChl))
        larg  = (1D0 + RMchl/mu0hat)*larg1   
        
        dmu0hat_aIdl   = (dmu0hatdl - alphaI*mu0hat)/aI
        d2mu0hat_aIdl2 = mu0hat/aI*2.*betamu 
     &  + (alphamu-alphaI+2.*betamu*PMU)* dmu0hat_aIdl

        daI_mu0hatdl = -(mu0hat/aI)**(-2)*dmu0hat_aIdl

        d2aI_mu0hatdl2 = -2.*betamu*aI/mu0hat
     &  - (2.*betamu*PMU + (alphamu-alphaI))*daI_mu0hatdl

        dlargdl = -RMchl*larg1/mu0hat**2 *dmu0hatdl
     &  + (1.+RMchl/mu0hat)*larg1 * par/zetaChl*daI_mu0hatdl
        
        d2largdl2 = -RMchl*(larg1*mu0hat**(-2)*d2mu0hatdl2
     &  + larg1*par/zetaChl*daI_mu0hatdl*mu0hat**(-2)*dmu0hatdl
     &  + larg1*dmu0hatdl*(-3.*mu0hat**(-3)*dmu0hatdl))
     &  + par/zetaChl*((1+RMchl/mu0hat)*larg1*d2aI_mu0hatdl2
     &  + (1.+RMchl/mu0hat)*larg1*par/zetaChl*daI_mu0hatdl*daI_mu0hatdl
     &  + RMchl*(-mu0hat**(-2)*dmu0hatdl)*larg1*daI_mu0hatdl)

       W_Y      = WAPR(larg,0,0)

       ThetaHat = 1.0/zetaChl + (1.0- W_Y) * mu0hat/(aI * par)
       
       dThetaHatdl = 1/par
     &  *  (-W_Y/larg/(1.+W_Y)*dlargdl*mu0hat/aI
     &  +  (1.-W_Y)*dmu0hat_aIdl)
       
       dWYYdl = dlargdl*(-W_Y**2/larg**2/(1.+W_Y)**3
     &  -  W_Y/larg**2/(1.+W_Y) + W_Y/(larg*(1.+W_Y))**2)

       d2ThetaHatdl2 = 1/par*(-(W_Y/larg/(1.+W_Y)*dlargdl*dmu0hat_aIdl
     &  +  W_Y/larg/(1.+W_Y)*d2largdl2*mu0hat/aI
     &  +  dWYYdl*dlargdl*mu0hat/aI)
     &  -  W_Y/larg/(1.+W_Y)*dlargdl * dmu0hat_aIdl
     &  +  (1.-W_Y)*d2mu0hat_aIdl2)

       else
          ThetaHat      = 0.01  !  a small positive value 
          dThetaHatdl   = 0.
          d2ThetaHatdl2 = 0.
       endif

        SI    = 1.-exp(-aI*par*ThetaHat/mu0hat)

        dSIdl = ( (alphaI-alphamu-2.*betamu*PMU)
     &   * ThetaHat + dThetaHatdl) * (1.-SI)*aI*par/mu0hat    !confirmed

       d2SIdl2 = par*(- dSIdl*aI/mu0hat*(ThetaHat*
     & (alphaI-alphamu-2.*betamu*PMU) + dThetaHatdl)
     & + (1.-SI)*(ThetaHat*(alphaI-alphamu-2.*betamu*PMU) + dThetaHatdl)
     & * (aI*alphaI/mu0hat-aI/mu0hat**2*dmu0hatdl)
     & + (1.-SI)*aI/mu0hat*(-2.*betamu*ThetaHat
     & + (alphaI-alphamu-2.*betamu*PMU)*dThetaHatdl+ d2ThetaHatdl2))

      ! Light dependent growth rate 
      ! (needs to take into account the cost of dark and light-dependent chl maintenance)
      mu0hatSI = mu0hat*SI  ! Gross specific carbon uptake (photosynthesis)
      muIhat   = mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat ! Net specific carbon uptake

       dmu0hatSIdl = mu0hatSI*(alphamu+2.*betamu*PMU) + mu0hat*dSIdl

       d2mu0hatSIdl2 = d2mu0hatdl2*SI+2.*dmu0hatdl*dSIdl+mu0hat*d2SIdl2

       dmuIhatdl = (1D0- zetaChl*ThetaHat)*dmu0hatSIdl  
     & - dThetaHatdl*zetaChl*(mu0hatSI+RMchl) 
   
       d2muIhatdl2=d2mu0hatSIdl2 - zetaChl 
     &  *(ThetaHat*d2mu0hatSIdl2+2.*dThetaHatdl*dmu0hatSIdl
     &  +mu0hatSI*d2ThetaHatdl2)
     &  -zetaChl*RMchl*d2ThetaHatdl2
    
      !Define VNhat: the saturation function of ambient nutrient concentration
       selectcase(self%nutrient_uptake)  
      ! case 1: Classic Michaelis Menton 
        case(1)
      ! Potential maximal nutrient-limited uptake
          alphaK = self%pars(ialphaK) 
          alphaV = self%pars(ialphaV)
          ! Half-saturation constant of nitrate uptake
          Kn    = ScaleTrait(PMU,self%pars(iK0), alphaK ) 
          VNhat = V0hat*NO3/(NO3 + Kn)
       dVNhatdl = -VNhat*alphaK*Kn/(NO3+Kn) + alphaV*VNhat

          d2VNhatdl2 = -alphaK*(Kn/(NO3+Kn)*dVNhatdl    
     &   - Kn*VNhat*alphaK* NO3/(NO3+alphaK)**2)
     &   + alphaV*dVNhatdl

      ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
        case(2)
         alphaA= self%pars(ialphaA)
         alphaV= self%pars(ialphaV)

          Lmin = log(0.5**3/6*pi) + log(1+Penfac)/( Penfac*alphaA) 

          A0hat  = PenAff(PMU, alphaA, Penfac, Lmin )
     &        * ScaleTrait(PMU, tf*self%pars(iA0), alphaA)

       dA0hatdl = alphaA*A0hat 
     &        - self%pars(iA0)*exp(PMU*alphaA)*Penfac*alphaA
     &        * exp(Penfac* alphaA *(PMU-Lmin))

       d2A0hatdl2 = alphaA*dA0hatdl 
     &        - Penfac*alphaA*exp(alphaA*((1.+Penfac)*PMU-Penfac*Lmin)) 
     &        * (dA0hatdl + A0hat*alphaA*(1.+Penfac))  
   
          !Define fA
          fA = 1.0/( 1.0 + sqrt(A0hat * NO3/V0hat) ) 
       VNhat = (1.0-fA)*V0hat*fA*A0hat*NO3/
     &        ((1.0-fA)*V0hat + fA*A0hat*NO3) 
   
      dVNhatdl= alphaV*VNhat
     &   - VNhat**2/NO3/A0hat*(1.+sqrt(NO3*A0hat/V0hat))
     &   * (alphaV-dA0hatdl/A0hat)  !Should be OK
        
        d2VNhatdl2 = alphaV*dVNhatdl-
     &   (VNhat**2/NO3/A0hat*(1.+sqrt(NO3*A0hat/V0hat))
     &   * (dA0hatdl**2/A0hat**2-d2A0hatdl2/A0hat) 
     &   + VNhat**2/NO3/A0hat*(alphaV-dA0hatdl/A0hat)
     &   * 0.5*sqrt(NO3/A0hat/V0hat)*(dA0hatdl-A0hat*alphaV) 
     &   + (1.+sqrt(NO3*A0hat/V0hat))*(alphaV-dA0hatdl/A0hat)
     &   * VNhat**2/NO3/A0hat*(2.*alphaV-dA0hatdl/A0hat))  !Should be OK
      endselect  

           ZINT = Qs*(muIhat/VNhat + zetaN)
           Qopt = (1.0 + sqrt(1.0 + 1.0/ZINT))*Qs
         dZINdl = Qs*(dmuIhatdl/VNhat - muIhat*dVNhatdl/VNhat**2)
     &           +  alphaQ*ZINT    
   
       dQoptdl  = alphaQ*Qopt-dZINdl*Qs/(2.*ZINT*sqrt(ZINT*(1.+ZINT))) !confirmed  
   
       d2ZINdl2 = Qs/VNhat*((alphaQ-dVNhatdl/VNhat)*dmuIhatdl 
     &  + d2muIhatdl2)
     &  - Qs/VNhat**2*(muIhat*d2VNhatdl2 
     &  + dVNhatdl*(dmuIhatdl+alphaQ*muIhat-2.*muIhat
     &  / VNhat*dVNhatdl)) + alphaQ*dZINdl
     
        d2Qoptdl2 = alphaQ*dQoptdl -  Qs/(2.*ZINT*sqrt(ZINT*(ZINT+1.)))
     & *  (d2ZINdl2 + alphaQ*dZINdl     
     & -  (2.*ZINT  + 1.5)/(ZINT*(ZINT+1.))*dZINdl**2)      ! Confirmed
    
         if (par .gt. I_zero) then   
            fV = (-1.0 + sqrt(1.0 + 1.0/ZINT))*Qs*muIhat/VNhat
!
         dfVdl = alphaQ*Qs*(1/Qopt + 2*zetaN)   
     &             -  (zetaN + Qs/Qopt**2) * dQoptdl  !Confirmed
!
           d2fVdl2 = (alphaQ**2)*Qs*(1/Qopt + 2*zetaN ) 
     &      -    2.*alphaQ*Qs* dQoptdl/Qopt**2    
     &      +    2.*(dQoptdl**2)*Qs/Qopt**3           
     &      -       (zetaN + Qs/Qopt**2.) * d2Qoptdl2  ! Confirmed
! 
           else 
            fV   = 0.01
          dfVdl  = 0.
         d2fVdl2 = 0.
         endif
      self%fV = fV
! Net growth rate (d-1) of phytoplankton at the average size as a function of temperature, light, and nutrient   
      self%muNet = muIhat*(1.0-fV-Qs/Qopt) - zetaN*fV*VNhat
! Here the derivative of muNet includes respiratory costs of both N Assim and Chl maintenance       
        self%dmuNetdl = muIhat*(Qs/(Qopt**2)*dQoptdl
     & -          alphaQ*Qs/Qopt-dfVdl) + (1.-fV-Qs/Qopt)*dmuIhatdl 
     & -          zetaN*(fV*dVNhatdl+VNhat*dfVdl)

         self%d2muNetdl2 = Qs/(Qopt**2)*dQoptdl*dmuIhatdl
     &  + muIhat*dQoptdl*(alphaQ*Qs/(Qopt**2)
     &  - 2.*Qs/(Qopt**3)*dQoptdl)
     &  + muIhat*Qs/(Qopt**2)*d2Qoptdl2   
     &  - (alphaQ*(Qs/Qopt*(dmuIhatdl+alphaQ*muIhat) 
     &  - muIhat*Qs/(Qopt**2)*dQoptdl))  
     &  - (muIhat*d2fVdl2+dfVdl*dmuIhatdl)  
     &  + dmuIhatdl*(Qs/(Qopt**2)*dQoptdl-alphaQ*Qs/Qopt-dfVdl) 
     &  + (1.-fV-Qs/Qopt)*d2muIhatdl2     
     &  - zetaN*(fV*d2VNhatdl2+2.*dfVdl*dVNhatdl+VNhat*d2fVdl2)  !dC

!  chl:C ratio [ g chl / mol C ] of the whole cell, at the mean cell size 
          Theta = ThetaHat*(1.0 - fV - Qs/Qopt)
! Calculate the mean chl:C ratio of phytoplankton (averaged over all sizes) 
! using the second derivative of the cell quota w.r.t. log size. (Why negative?)
            self%ThetaAvg = Theta + VAR*.5          
     &    *  ( d2ThetaHatdl2*(1.0 - fV - Qs/Qopt)         
     &    -      ThetaHat   *(d2fVdl2 + Qs*((dQoptdl**2)*(2./Qopt)
     &    - d2Qoptdl2)/(Qopt**2)))
            
            if (self%ThetaAvg .lt. 0) self%ThetaAvg = Theta
   
! Calculate the mean N:C ratio of phytoplankton (averaged over all sizes) 
! using the second derivative of the cell quota w.r.t. log size. (How to derive?)
       self%QNavg = Qopt/(1.+((2/Qopt)*dQoptdl**2 - d2Qoptdl2)
     &    * VAR/(2.*Qopt))  
      ! Calculate the community sinking rate of phytoplankton
      ! Phytoplankton sinking rate at the average size
      w_p         = ScaleTrait(PMU,self%pars(iw_p),self%pars(ialphaW))
      self%dwdl   = self%pars(ialphaW)*w_p
      self%d2wdl2 = self%pars(ialphaW)**2*w_p
      self%w_pAvg = w_p+0.5*VAR*self%d2wdl2 ! Sinking rate of phytoplankton

#endif /* EFTSIZE */
      end subroutine phyto
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
 
!---------------------------------------------------------
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
!
!-----------------------------------------------------------------------
! Adopted from the model of Lan. Smith
      PURE REAL FUNCTION ScaleTrait( logsize, star, alpha) 
          implicit none
          real, intent(IN) :: logsize, star, alpha
       
          ! Calculate the size-scaled value of a trait
          ! for the given log (natural, base e) of cell volume as pi/6*ESD**3 (micrometers). 
       
          ScaleTrait = star * exp( alpha * logsize )
          return
      END FUNCTION ScaleTrait
!-----------------------------------------------------------------------    
      subroutine Zooplankton(self)
! Calculate the feeding probability Q(Pt) based on the optimal feeding current model 
! (Pahlow MEPS 2010)
! The total amount of phytoplankton grazed by zooplankton (molN m-3)
! This does not depend on size classes and the parameter 'alphaG'
! that determines prey switching

! This subroutine calculates the ingestion, respiration, egestion, and net growth rate
! of zooplankton given the total amount of prey available 
! Change this subroutine to include prey size dependence of food palatability 
! Based on Smith and Vallina (in prep.)
! Also ignore the optimal feeding model by Pahlow (2010)
       implicit none
       class (type_PZEFT), intent(inout), target :: self
       ! Declare local variables
       real :: tf,ghat,g_th,aAt,aAf,AE,PHY,VAR,kp,PMU 
       integer :: i
       integer, parameter :: N=100    !discretize the continous normal distribution
       real :: x(N), Ptot, Ptot1,dx
       real, parameter :: lmin = log(pi/6D0*2**3)
       ! For control model:
       real, parameter :: GGE = 0.3, unass = 0.24

          tf  = self%tf()
          PHY = self%vars(iPHY)
          kp  = self%pars(ikp)
           
#if defined NPZDSIMPLE || defined EFTNOSIZE 
        ! Calculate the feeding probability function based on self%PHY
        ! kp relates to the capture coefficient
          SELECT CASE(self%grazing_formulation)
            ! Holling Type I
            case (1)
            ghat = min(PHY/2.0/kp,1.0)
            ! Holling Type II
            case (2)
            ghat = PHY/(kp + PHY)  
            ! Holling Type III
            case (3) 
            ghat = PHY**2/(kp**2 + PHY**2)
            ! Ivlev
            case (4)
            ghat = 1D0 - exp(-LOG(2D0)*PHY/kp)  
          END SELECT

          self%INGES = tf*self%pars(igmax)*ghat
          self%RES   = self%INGES*(1D0-GGE-unass)
          self%NGR   = self%INGES*GGE
          ! ZOOPLANKTON EGESTION (-> Detritus)
          self%EGES = self%INGES - self%NGR - self%RES

#else
          PMU = self%PMU
          VAR = self%VAR
         ! Calculate the integration  
          x(1) = PMU-3D0*sqrt(VAR)  ! Minimal size
          x(N) = PMU+3D0*sqrt(VAR)  ! Maximal size
          dx   = (x(N) - x(1))/float(N-1)
          do i=2,N-1
            x(i) = x(i-1) + dx 
          enddo
        ! Calculate total phytoplankton palatable prey (Ptot,Eq. 13):
        ! Ptot roughly equals to Phy
         Ptot = 0.
         Ptot1= 0.
         do i=1,N-1

           Ptot = Ptot 
     &   + normal(PMU,VAR,x(i))*dx*exp(self%pars(iS)*(x(i)-lmin))
           ! Considering killing-the-winner strategy
           Ptot1= Ptot1 + normal(PMU,VAR,x(i))**self%pars(ialphaG)
     &              *dx*exp(self%pars(iS)*(x(i)-lmin))
         enddo
         Ptot  =  Ptot*PHY*1.0027  ! Total palatable prey
         Ptot1 = Ptot1*PHY*1.0027  ! denominator of gbar
!         write(*,*) 'Ptot = ',Ptot,'Phy = ',PHY, 'Ptot1= ',Ptot1
        ! Calculate the feeding probability function based on Ptot
        ! kp relates to the capture coefficient
!         SELECT CASE(self%grazing_formulation)
!           ! Holling Type I
!           case (1)
!           ghat = min(Ptot/2.0/kp,1.0)
!           ! Holling Type II
!           case (2)
!           ghat = Ptot/(kp + Ptot)  
!           ! Holling Type III
!           case (3) 
           ghat = Ptot**2/(kp**2 + Ptot**2)
           ! Ivlev
!           case (4)
!           ghat = 1.0 - exp(-LOG(2D0)*Ptot/kp)  
!         END SELECT

#  if defined OCF && defined EFTSIZE
        ! Calculate feeding threshold of ghat (dimensionless)
          g_th = self%pars(iCf)/self%pars(iEmax)/(1.0 - self%pars(iCa))
   
        ! Calculate total foraging activity (d-1)
          aAt = tf*self%pars(igmax)/self%pars(ibeta)*      
     &      (-1.0 - w1(-(1.0-g_th)/EXP(1.0+self%pars(ibeta))))
   
        ! Compute the optimal foraging activity to maximize net growth of zooplankton
          IF (ghat .GT. g_th) THEN
           ! Optimal foraging, aAf: Specific foraging activity (d-1) = 1/h (h: biomass specific handling time)
           !With increasing food conc., optimal aAf slightly increases with food conc.
              aAf = aAt*self%pars(ibeta)/(-1.0 - w1(-(1.0-g_th/ghat) 
     &         /EXP(1.0 + self%pars(ibeta))))
         ! Assimilation efficiency relates with two factors
         ! 1) The difference between total activity (aAt) and foraging activity (aAf)
         ! 2) Gut passage time
         ! AE is negatively related with foraging activity
         ! because digestion enzymes decrease with increasing food
         ! Another factor, gut passange time linearly increases with handling time
            AE  = self%pars(iEmax)*
     &     (1.0 - exp(-self%pars(ibeta)*(aAt/aAf - 1.0)))  
         ! If not optimal grazing, set AE = 0.45 (aAf/aAt = 0.6)    
            ELSE
   
            aAf = 0.0
            AE  = self%pars(iEmax)
          ENDIF
        
        !Ingestion rate (d-1) of zooplankton (gmax*Q(PT)), now should be a function of size 
          self%INGES = aAf * ghat 
          
          ! Zooplankton respiration (-> NO3)
          self%RES = self%pars(iRMzoo)*tf
     & +             self%pars(iCa)*AE*self%INGES
     & +             self%pars(iCf)*aAf    
          ! Zooplankton net growth (To force it positive)
          self%NGR = AE*self%INGES - self%RES
          ! ZOOPLANKTON EGESTION (-> Detritus)
          self%EGES = self%INGES - self%NGR - self%RES
#  else
          !Zooplankton mass-specific ingestion rate (= m*Phy/Zoo in dilution exp., unit: d–1)
          self%INGES = tf*self%pars(igmax)*ghat
        

          self%RES   = self%INGES*(1.-GGE-unass)
          self%NGR   = self%INGES*GGE
         ! ZOOPLANKTON EGESTION (-> Detritus)
          self%EGES  = self%INGES - self%NGR - self%RES

!          
#  endif /* OCF  */
! Grazing rate on the mean size (specific to N-based Phy biomass, unit: d-1) (Eq. 12)
          self%gbar = self%INGES*self%vars(iZOO)
     &  *   exp(self%pars(iS)*(PMU - lmin))
     &  *   normal(PMU,VAR,PMU)**(self%pars(ialphaG)-1d0)/Ptot1 

          self%dgdlbar   = self%pars(iS)*self%gbar  ! for all cases (when evaluated at the mean size

          self%d2gdl2bar = -(self%pars(ialphaG)-1d0)*self%gbar/VAR 
     &  +     self%pars(iS)**2*self%gbar 
#endif
      end subroutine Zooplankton
!--------------------------------------------------------------------------------------
      real function normal(mean,var,l)
        implicit none
        real, intent(in) :: mean,var,l
        normal = 1D0/sqrt(2D0*var*pi)*exp(-(l-mean)**2/var/2D0) 
        return
      end function normal
!--------------------------------------------------------------------------------------
!Calculate the initial estimate of first and second derivatives of ghat based on Richardson extrapolation
!      subroutine INI_DERIV(self,h)
!        implicit none
!!INPUT PARAMETERS:
!        real, intent(IN) :: h
!        class (type_PZEFT),intent(inout) :: self
!!The number of variables that need to calculate 1st and 2nd derivatives
!        !Array storing f(x+h),f(x),f(x-h) at h,h/2,h/4
!        real,allocatable :: f(:,:)
!        real :: l,PMU
!        integer :: i
!
!       allocate(f(NDERIV,3))
!!First, Get an initial estimate of first and second derivatives:
!!A0(h) = (f(x+h) - f(x-h))/(2*h) , and then:
!       !Get Average size:
!       PMU = log(self%vars(iPMUPHY)/self%vars(iPHY))
!       l   = PMU - h
!!1:f(x-h);2:f(x);3:f(x+h)
!       do i = 1,3
!          l = l + h
!          self%vars(iPMUPHY) = self%vars(iPHY)*exp(l)
!          call self%Zooplankton
!          f(1,i)   = self%gbar
!       enddo
!
!!Estimate the initial values of the first derivatives:
!       self%deriv(:,1) = (f(:,3) - f(:,1))/(2d0*h)
!!For second derivatives:
!!A0(h) = (f(x+h) - 2*f(x) + f(x-h))/h**2       
!       self%deriv(:,2) = (f(:,3)-2d0*f(:,2)+f(:,1))/h**2
!!restore PMU:
!       self%vars(iPMUPHY) = self%vars(iPHY)*exp(PMU)
!       deallocate(f)
!      end subroutine INI_DERIV
!-------------------------------------------------------------------------------------------
!      subroutine RICHARDDERIV(self)
!Calculate approximate first and second derivatives based on Richardson extrapolation
!The initial estimate of first and second derivatives are called from ini_deriv
!A1(h) = (4**n*A0(h/2) - A0(h))/(4**n-1), then:
!The code copied from wikipedia:
!https://en.wikipedia.org/wiki/Richardson_extrapolation
!       implicit none
!       class (type_PZEFT),intent(inout) :: self
!!       integer, parameter :: M = 20  !Maximal Nrows
!       real :: h = 0.05
!       !Matrix storing Richardson extrapolations
!!       real, allocatable :: A(:,:,:,:)
!       real :: A(2,NDERIV,2)
!!       real, parameter :: tol = 1d-4  !tolerance
!!       integer :: i,j
!
!!       allocate(A(2,NDERIV,2))
!       call self%INI_DERIV(h)
!       A(1,:,:)  = self%deriv(:,:)
!!       A(1,1,:,:) = self%deriv(:,:)
!!       iloop: do i=1,M-1
!       h = h/2d0
!       call self%INI_DERIV(h)
!!       A(i+1,1,:,:) = self%deriv(:,:)
!       A(2,:,:) = self%deriv(:,:)
!       self%deriv(:,:) = 2D0*A(2,:,:) - A(1,:,:)
!         jloop: do j = 1,i !Go across the row until the diagonal is reached
!Use the last value computed (i.e. A(i + 1, j)) and 
!the element from the row above it (i.e. A(i, j)) 
!to compute the next Richardson extrapolate
!         A(i+1,j+1,:,:)=(2**j*A(i+1,j,:,:) 
!     &                      - A(i,  j,:,:))/(2**j-1d0)
!
!         enddo jloop

!         if(all(abs(A(i+1,i+1,:,:)-A(i,i,:,:)) .lt. tol )) then
          !Write final derivatives:
!              self%deriv(:,:) = A(i+1,i+1,:,:) 
!              Findsolution = .true.
!              exit iloop
!         endif
!        enddo iloop

!        if (Findsolution .eq. .false.) then
!          write(stdout,*) 'Failed to find solution.'
!          self%deriv(:,:) = A(10,10,:,:)
!        endif
!       deallocate(A)
!
!       self%dgdlbar   = self%deriv(1,1)
!       self%d2gdl2bar = self%deriv(1,2)
!
!      end subroutine RICHARDDERIV
!======================================================================
      subroutine ODE(self)
       implicit none
       class(type_PZEFT),intent(inout),target :: self
       real :: rpd,tf,INGES
       integer, parameter :: numc=4
       integer :: i,j
       real :: PP_nz,PP_nd,PP_dp,PP_dz,PP_pn,PP_zp,PP_np,
     &         NO3,PHY,ZOO,DET,VAR,PMU,dPMUdt,dVARdt,mort,dPMUPHYdt,
     &         dVARPHYdt

      tf   = self%tf()
      NO3  = self%vars(iNO3)
      DET  = self%vars(iDET)
      ZOO  = self%vars(iZOO)
      PHY  = self%vars(iPHY)
      PMU  = self%PMU
      VAR  = self%VAR
      rpd  = tf*self%pars(irpd)
! dz = zooplankton growth(excluding excretion and egestion) - zooplankton mortality (molN d-1)
!      selectcase(self%zoo_mortality)
!        case(1)
!        ! Linear mortality term
!        mort = ZOO*self%pars(imz)*tf
!        case (2)
        ! Quadratic term
        mort = ZOO**2*self%pars(imz)*tf
        
!        case (3)
!        ! MM Type II term 
!        mort = tf*ZOO*self%pars(imz)*ZOO/(ZOO+self%pars(ikzb)) 
!
!        case (4)
!        ! Type III term
!        mort = tf*ZOO*self%pars(imz)*ZOO**2
!     &        /(ZOO**2 + (self%pars(ikzb))**2)
!      endselect

     
   
      call self%phyto
      call self%Zooplankton

      INGES = self%INGES
#if defined NPZDSIMPLE || defined EFTNOSIZE
      self%odeout(iPHY) = PHY*(self%muNet -rpd) - ZOO*INGES
      ! Nutrient changes (Eq. 23)
      self%odeout(iNO3) = self%RES*ZOO + self%pars(irdn)*DET*tf
     &  - PHY*self%muNet
#else
      self%odeout(iPHY) = 0.5*PHY*(2D0*(self%muNet - rpd) 
     &  + VAR*self%d2muNetdl2) - ZOO*INGES  ! The total amount of phytoplankton grazed by zoo is independent of average size    

      ! Nutrient changes (Eq. 23)
      self%odeout(iNO3) = self%RES*ZOO + self%pars(irdn)*DET*tf
     &  - PHY*(self%muNet + 0.5*VAR*self%d2muNetdl2)
!
!! Since we treat the product of (PMU*PHY) as a single unit in the differential eq.
!! d(PMUPHY)dt = PMU*dPHYdt + PHY*dPMUdt (here the PMU is exp(PMU)), the same applies to d(VARPHY)dt
!!      call self%RICHARDDERIV

      !Eq. 18,sinking considered.
      dPMUdt = VAR   *(self%dmuNetdl   - self%dgdlbar   - self%dwdl)  
      dVARdt = VAR**2*(self%d2muNetdl2 - self%d2gdl2bar - self%d2wdl2)  !Eq. 19 

      dPMUPHYdt=exp(PMU)*(self%odeout(iPHY) + PHY*dPMUdt)
      dVARPHYdt=    VAR * self%odeout(iPHY) + PHY*dVARdt
#endif

      self%odeout(iZOO) = ZOO*self%NGR - mort
! dD = phytoplankton mortality + zooplankton mortality + zooplankton egestion - conversion-to-nutrient (molN d-1)
      self%odeout(iDET) = -(self%odeout(iNO3)+self%odeout(iPHY)
     &        + self%odeout(iZOO))
! For PPDD matrix:
      PP_nz = self%RES*ZOO
      PP_nd = self%pars(irdn)*DET*tf
      PP_dp = PHY*rpd
      PP_dz = mort + ZOO*self%EGES
      PP_zp = ZOO*self%INGES
      PP_pn = 
#if defined NPZDSIMPLE || defined EFTNOSIZE
     &   PHY*self%muNet      
#else
     &   PHY*(self%muNet + 0.5*VAR*self%d2muNetdl2)

      PP_pn = 0.5*( PP_pn + abs(PP_pn))
      PP_np = 0.5*(-PP_pn + abs(PP_pn))
!      if (PP_pn .lt. 0) then
!        PP_pn = 0. 
!        PP_np = -PP_pn
!      endif
#endif
! Write out production/destruction matrix
      ! 1: NO3; 2: PHY; 3: ZOO; 4: DET
      self%pp(:,:) = 0.
      self%pp(1,2) = max(PP_np,0.)
      self%pp(1,3) = max(PP_nz,0.)
      self%pp(1,4) = max(PP_nd,0.)
      self%pp(2,1) = max(PP_pn,0.)
      self%pp(3,2) = max(PP_zp,0.)
      self%pp(4,2) = max(PP_dp,0.)
      self%pp(4,3) = max(PP_dz,0.)
#ifndef EULER      
      self%dd(:,:) = 0.
      do i = 1,numc
        do j = 1,numc
          self%dd(j,i)=self%pp(i,j) 
        end do
      end do 
      self%odeout(iPMUPHY)=dPMUPHYdt
      self%odeout(iVARPHY)=dVARPHYdt
#endif
#ifdef EULER      
      self%vars(iNO3) = NO3*(1.+(PP_np+PP_nz+PP_nd)/NO3)/(1.+PP_pn/NO3)
      self%vars(iPHY) = PHY*(1.+PP_pn/PHY)/(1.+(PP_zp+PP_dp+PP_np)/PHY)
      self%vars(iZOO) = ZOO*(1.+PP_zp/ZOO)/(1.+(PP_nz+PP_dz)/ZOO)
      self%vars(iDET) = DET*(1.+(PP_dp+PP_dz)/DET)/(1.+PP_nd/DET)
      if (dPMUPHYdt .ge. 0.) then
          self%vars(iPMUPHY) = self%vars(iPMUPHY) + dPMUPHYdt
      else
          self%vars(iPMUPHY) = self%vars(iPMUPHY)
     &     /(1.-dPMUPHYdt/self%vars(iPMUPHY)) 
      endif

      if (dVARPHYdt .ge. 0.) then
          self%vars(iVARPHY) = self%vars(iVARPHY) + dVARPHYdt
      else
          self%vars(iVARPHY) = self%vars(iVARPHY)
     &     /(1.-dVARPHYdt/self%vars(iVARPHY)) 
      endif
#endif
      end subroutine ODE
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
      subroutine modified_patankar(self)

        IMPLICIT NONE
!
! !INPUT PARAMETERS:
       integer, parameter  :: numc=4,stdout=6
       class(type_PZEFT),intent(inout),target :: self
!
! !INPUT/OUTPUT PARAMETER:
        REAL :: cc(1:numc)
!     
!! !REVISION HISTORY:
!!  Original author(s): Hans Burchard, Karsten Bolding
!!
!! !LOCAL VARIABLES:
       REAL :: pp1(1:numc,1:numc),dd1(1:numc,1:numc)
       REAL :: a(1:numc,1:numc),r(1:numc)
       REAL :: cc1(1:numc)
       integer  :: i,j
!!EOP
!!-----------------------------------------------------------------------
!!BOC
      do i = 1,numc
         cc(i) = self%vars(i)
      enddo
!  absolutely essential since not all elements are calculated
        pp1(:,:)=0.
        dd1(:,:)=0.
     
        call self%ODE
         
        pp1 = self%pp
        dd1 = self%dd

      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd1(i,j)
            if (i.ne.j) a(i,j)=-pp1(i,j)/cc(j)
         end do
         a(i,i)=a(i,i)/cc(i)
         a(i,i)=1.+a(i,i)
         r(i)  =cc(i)+pp1(i,i)
      end do

      call matrix(numc,a,r,cc1(:))
      
      self%vars(iNO3)   =cc1(iNO3)
      self%vars(iPHY)   =cc1(iPHY)
      self%vars(iZOO)   =cc1(iZOO)
      self%vars(iDET)   =cc1(iDET) 
      self%vars(iPMUPHY)=self%vars(iPMUPHY) +self%odeout(iPMUPHY)
      self%vars(iVARPHY)=self%vars(iVARPHY) +self%odeout(iVARPHY)

      end subroutine modified_patankar
!--------------------------------------------------------------------------------------     
!--------------------------------------------------------------------------------------     
      subroutine RK2(self)
        implicit none
        class(type_PZEFT),intent(inout),target :: self
        integer :: i
        real :: newode(NVAR) ! temporary array holding intermediate values
        real :: vars(NVAR)

      do i = 1, NVAR
        vars(i)=self%vars(i)
      enddo
     
      call self%ODE

      do i = 1,NVAR
        newode(i)    = self%odeout(i)
        self%vars(i) = self%vars(i) + self%odeout(i)
      enddo 

      call self%ODE

      do i = 1,NVAR
        self%vars(i) = vars(i)+0.5*(self%odeout(i)+newode(i))
      enddo 
      
      end subroutine
!--------------------------------------------------------------------------------------------
      subroutine modified_patankar_2(self)

        IMPLICIT NONE
!
! !INPUT PARAMETERS:
       integer, parameter  :: numc=4,stdout=6
       class(type_PZEFT),intent(inout),target :: self
!
! !INPUT/OUTPUT PARAMETER:
        REAL :: cc(1:numc)
!     
!! !REVISION HISTORY:
!!  Original author(s): Hans Burchard, Karsten Bolding
!!
!! !LOCAL VARIABLES:
       REAL :: pp2(1:numc,1:numc),dd2(1:numc,1:numc)
       REAL :: pp1(1:numc,1:numc),dd1(1:numc,1:numc)
       REAL :: a(1:numc,1:numc),r(1:numc)
       REAL :: cc1(1:numc)
       real :: dPMUPHYdt1,dVARPHYdt1,PMUPHY,VARPHY
       integer  :: i,j
!!EOP
!!-----------------------------------------------------------------------
!!BOC
      do i = 1,NVAR
         self%vars(i) = max(self%vars(i),eps)   !To avoid negative values
      enddo
      do i = 1,numc
         cc(i) = self%vars(i)
      enddo
        PMUPHY = self%vars(iPMUPHY)
        VARPHY = self%vars(iVARPHY)
!  absolutely essential since not all elements are calculated
        pp2(:,:)=0. ; pp1(:,:)=0.
        dd2(:,:)=0. ; dd1(:,:)=0.
     
        call self%ODE
         
        pp1        = self%pp
        dd1        = self%dd
        dPMUPHYdt1 = self%odeout(iPMUPHY)
        dVARPHYdt1 = self%odeout(iVARPHY)

      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            a(i,i)=a(i,i)+dd1(i,j)
            if (i.ne.j) a(i,j)=-pp1(i,j)/cc(j)
         end do
         a(i,i)=a(i,i)/cc(i)
         a(i,i)=1.+a(i,i)
         r(i)  =cc(i)+pp1(i,i)
      end do

      call matrix(numc,a,r,cc1(:))
      
      self%vars(iNO3)   =cc1(iNO3)
      self%vars(iPHY)   =cc1(iPHY)
      self%vars(iZOO)   =cc1(iZOO)
      self%vars(iDET)   =cc1(iDET) 
      self%vars(iPMUPHY)=self%vars(iPMUPHY) +self%odeout(iPMUPHY)
      self%vars(iVARPHY)=self%vars(iVARPHY) +self%odeout(iVARPHY)

      call self%ODE
     
      pp2=0.5*(self%pp+pp1)
      dd2=0.5*(self%dd+dd1)
     
      do i=1,numc
         a(i,i)=0.
         do j=1,numc
            
            a(i,i)=a(i,i)+dd2(i,j)
            if (i.ne.j) a(i,j)=-pp2(i,j)/cc1(j)
         enddo
         a(i,i)=a(i,i)/cc1(i)
         a(i,i)=1.+a(i,i)
         r(i)  =cc(i)+pp2(i,i)
      enddo

      call matrix(numc,a,r,cc(:))
     
      self%vars(iNO3)=cc(1)
      self%vars(iPHY)=cc(2)
      self%vars(iZOO)=cc(3)
      self%vars(iDET)=cc(4) 
      self%vars(iPMUPHY)=PMUPHY 
     & + 0.5*(self%odeout(iPMUPHY) + dPMUPHYdt1)
      self%vars(iVARPHY)=VARPHY
     & + 0.5*(self%odeout(iVARPHY) + dVARPHYdt1)

      end subroutine modified_patankar_2
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Matrix solver
!
! !INTERFACE:
        subroutine matrix(n,a,r,c)
!!
!! !DESCRIPTION:
!! This is a Gaussian solver for multi-dimensional linear equations.
!!
!! !USES:
        IMPLICIT NONE
!!
!! !INPUT PARAMETERS:
        integer, intent(in)             :: n
!!
!! INPUT/OUTPUT PARAMETERS:
       REAL                             :: a(1:n,1:n),r(1:n)
!!
!! OUTPUT PARAMETERS:
       REAL, intent(out)                :: c(1:n)
!!
!! !REVISION HISTORY:
!!  Original author(s): Hans Burchard, Karsten Bolding
!!
!! !LOCAL VARIABLES:
       integer  :: i,j,k
!!EOP
!!-----------------------------------------------------------------------
!!BOC
        do i=1,n
           r(i)=r(i)/a(i,i)
           do j=n,i,-1
              a(i,j)=a(i,j)/a(i,i)
           end do
           do k=i+1,n
              r(k)=r(k)-a(k,i)*r(i)
              do j=i+1,n
                 a(k,j)=a(k,j)-a(k,i)*a(i,j)
              end do
           end do
        end do
     
        do i=n,1,-1
           c(i)=r(i)
           do j=i+1,n
              c(i)=c(i)-a(i,j)*c(j)
           end do
        end do
     
        return
        end subroutine matrix
!--------------------------------------------------------------------------------------
      END MODULE PZEFT
#endif
