! THis is a module based on an Optimality-based NPZD model that can couple with
! other physical (0D,1D,3D) hosts.
! The optimal phytoplankton growth part based on Pahlow et al. MEPS 2013.
! The optimal zooplankton grazing part based on Pahlow & Prowe MEPS 2010.
! Original structure taken from GOTM and adapted for FABM by Jorn Bruggeman. 
! !Add the component of optimal current feeding on 20150309 
! Modify the unimodal size functional relationship of growth rate in July, 2015 
      MODULE PZEFT
      use lambert
      implicit none
      private
      integer,parameter :: NENV=11, NPAR=33
! Define potential minimal and maximal values to confine the values of state variables
      real,parameter,public :: 
     &  NO3min    = 1E-6,
     &  NO3max    = 1d2,
     &  PHYmin    = 1E-8,
     &  PHYmax    = 2d1,
     &  ZOOmin    = 1E-10,
     &  ZOOmax    = 2d1,
     &  DETmin    = 1E-6,
     &  DETmax    = 1d3,
     &  PMUmin    = -2.1,   ! roughly ESD 0.6 µm
     &  PMUmax    = 7d0,    ! roughly ESD ? µm
     &  VARmin    = 1d-3,
     &  VARmax    = 2d1,
     &  PMUPHYmin = 1d-9,
     &  PMUPHYmax = 2d4,
     &  VARPHYmin = 1d-9,
     &  VARPHYmax = 4d2

! PUBLIC DERIVED TYPES:
      type,public :: type_PZEFT   ! Define a type
!     Variable identifiers
! Input variables: (These variables must be provided in the host model)
!-------------------------------------------------------------------------------
! Modify input data and parameters as vectors (to facilitate MCMC fitting)
!      real(8) :: NO3,PHY,ZOO,DET,PMU,VAR,CHL,par,tC,I_0,dtdays
      real :: env(NENV)
!-------------------------------------------------------------------------------
! The following parameters are internal and should be parameterized in the 'initialize' subroutine
      real :: pars(NPAR)
      integer :: nutrient_uptake,grazing_formulation
      logical :: fV_opt,CF_opt
      real :: fV_initial
! Output variables (to be used by the host model): 
!-------------------------------------------------------------------------------
      real :: Theta,ThetaAvg,QNavg,muNet,dmuNetdl,d2muNetdl2,fV,
     &        w_pAvg
      real :: AE,INGES,RES,NGR,EGES,gbar,dgdlbar,d2gdl2bar
      real :: dPMUPHYdt,dVARPHYdt,dPHYdt,dNO3dt,dZOOdt,dDETdt        ! For ODE
      real :: cc(4),pp(4,4),dd(4,4)        ! For PPDD
!-------------------------------------------------------------------------------
      contains
         procedure :: initialize
         procedure :: tf=>temp
         procedure :: growth
         procedure :: Zooplankton
         procedure :: ODE
         procedure :: modified_patankar_2
      end type type_PZEFT
      integer, parameter, public :: iNO3 =1,  iPHY = 2, iZOO = 3, 
     &                      iDET =4,
     &                      iPMUPHY =5,!iPMUPHY is the product of PHY and PMU (in biovolume)
                                      !that will enter diffusion and advection 
                                      !in the 3D physical model
     &                      iVARPHY =6,!iVARPHY is the product of PHY and cell variance
     &                      iCHL =7,  ipar = 8,
     &                      itC  =9,  iI_0 =10, idtdays=11

      integer, parameter, public :: imu0 =1,  ialphamu=2,
     &                           iaI0    =3,  ialphaI =4,
     &                           iA0     =5,  ialphaA =6,
     &                           iK0     =7,  ialphaK =8,
     &                           iV0     =9,  ialphaV =10,
     &                           iQ0     =11, ialphaQ =12,
     &                           iw_p    =13, iw_d    =14,
     &                           ikc     =15, igmax   =16,
     &                           ikp     =17, irdn    =18,
     &                           irpdu   =19, irpdl   =20,
     &                           iEa     =21, imz     =22, 
     &                           izetaN  =23, iRMchl  =24,
     &                           izetaChl=25, ibeta   =26,
     &                           iCa     =27, iCf     =28,
     &                           iEmax   =29, iRMzoo  =30,
     &                           ialphaG =31, ialphaW =32,ibetamu=33

      integer, parameter :: stdout= 6
      real,    parameter :: errr  = 1d-5, pi=3.141592653529D0 
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

      w31_p  = 0.1   !Stoke's law: V = 2/9*(rho_p - rho_f)*g*R**2/mu (mu: dynamic viscosity, which is function
                     !of temperature.) w_p0 = w(31)*(0.451 + 0.0178*T) (Gregg et al. 2008). 
                     !We also set this to be size dependent
                     !such that w_p = w_p0* Wstar * exp(logsize * alphaW)
      w31_d  = 1.0
      dtdays = self%env(idtdays) 
      temp   = self%env(itC    )
      self%pars(imu0    )    = 2.5D0*dtdays
      self%pars(ialphamu)    = 1.2d-1
      self%pars(ibetamu )    = -0.02
! Initial slope of the P-I curve (mol C gChl-1 (W m-2)-1 d-1). 
! 1 W m-2 = 4.6 uE m-2 s-1. 
! Ward et al. (L&O) used the size independent slope of 3.83E-7 mol C (gChl-1) (uE m-2)-1. 
! Edwards et al. L&O (2015) estimated a size dependent value of 10**(-1.36) d-1 (uE m-2 s-1)-1 (= 5.05E-7 (uE m-2)-1)
! So I choose a value of 10**(-1.36)*4.6/(50/12)= 0.05 

      self%pars(iaI0   )    = 2D-1 *dtdays
      self%pars(ialphaI)    = 0.0            ! -0.13 From Edwards et al. L&O 2015
      self%pars(iA0    )    = 2D1 *dtdays     ! From Pahlow et al. MEPS 2013    
      self%pars(ialphaA)    = -0.289
      self%pars(iV0    )    = 3d0*dtdays
      self%pars(ialphaV)    = -0.27

!Normalization constant for half-saturation constant for nitrate (mmolN m-3), from Litchman et al. 2007
      self%pars(iK0    )    = 0.17    !Not estimated
      self%pars(ialphaK)    = 0.27    !Not estimated
      self%pars(ialphaG)    = 1.10    !Not estimated
      self%pars(iQ0    )    = 0.07    ! 0.07 in Litchman et al. 2007; Ward et al. L&O 2012 
      self%pars(ialphaQ)    = -0.17   ! -0.17 in Litchman et al. 2007; Ward et al. L&O 2012 
      self%pars(iw_p   )    = w31_p*dtdays*(0.451+0.0178*temp) !Temperature dependent sinking rate of phytoplankton, POSITIVE!!
      self%pars(iw_d   )    = w31_d*dtdays*(0.451+0.0178*temp)!Sinking rate of detritus
      self%pars(ialphaW)    = 0.39       ! Laws 1975
      self%pars(ikc    )    = 0.025      !Light attenuation of chlorophyll (mg Chl m-2)
      self%pars(igmax  )    = 3.0*dtdays !Maximal biomass-specific ingestion rate of zooplankton at 15 ºC
      self%pars(ikp    )    = 8D-1        !half-saturation constant for zooplankton grazing
      self%pars(irdn   )    = 0.03*dtdays!The conversion rate from detritus to nutrient at reference temperature
      self%pars(irpdu  )    = 0.10*dtdays!Phytoplankton mortality rate at reference temperature in the euphotic zone
      self%pars(irpdl  )    = 0.02*dtdays!Phytoplankton mortality rate at reference temperature below the euphotic zone 
      self%pars(iEa     )   = 0.3             !Activation energy (eV)
      self%pars(imz     )   = 0.25*dtdays    !Coefficient for zooplankton mortality at 15 C
      self%pars(izetaN  )   = 7.29d-1             !Cost of DIN uptake (mol C/mol N)   
      self%pars(iRMchl  )   = 4.15d-5*dtdays     !Cost of chl maintenance (d-1) 
      self%pars(izetaChl)   = 1.5D0            !Cost of photosynthesis coefficient (mol C/g Chl)
      self%pars(ibeta   )   = 0.2             !Assimilation coefficient 
      self%pars(iCa     )   = 0.3             !Cost of assimilation coefficient
      self%pars(iCf     )   = 0.3             !Effective cost of foraging coefficient  
      self%pars(iEmax   )   = 0.99            !Maximal assimilation efficiency
      self%pars(iRMzoo  )   = 0.16*dtdays     !Specific maintenance respiration (d-1)
! End of pars definition.

! Choose Different models:
      self%fV_opt    = .true.     !logic variable indicating whether the model should be optimality-based
      self%nutrient_uptake=1      !The model used for nutrient uptake
      self%grazing_formulation=3
      self%CF_opt     = .true.
      self%fV_initial = 0.25 !Initial value for the fraction of energy devoted to nutrient uptake

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
       tf=exp(-(Ea/kb)*(1.0/(273.15 + tC)-1.0/(273.15 + Tr)))
     
       end function temp
!EOC
 
!-----------------------------------------------------------------------
! !IROUTINE:Computing phytoplankton growth components  
!
! !INTERFACE:
      subroutine growth(self)
!INPUT PARAMETERS:
      class (type_PZEFT),intent(inout) :: self
!LOCAL VARIABLES:
      real :: I_zero, tf, PMU,VAR,par,NO3,zetaChl,alphaI,alphamu,
     &        betamu
      real :: dmu0hatdl, d2mu0hatdl2
      real :: V0hat,Kn,Lmin,A0hat,dA0hatdl,d2A0hatdl2,fA
      real :: VNhat,dVNhatdl,d2VNhatdl2 ! Nutrient uptake variables
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
      PMU    = log(self%env(iPMUPHY)/self%env(iPHY)) ! Convert the unit to log(mu**3)
      PMU    = min(max(PMU,PMUmin),PMUmax)
      VAR    = self%env(iVARPHY)/self%env(iPHY)
      VAR    = min(max(VAR,VARmin),VARmax)
      par    = self%env(ipar)
      NO3    = self%env(iNO3)
      NO3    = min(max(NO3,NO3min),NO3max)
      zetaChl= self%pars(izetaChl)
      alphaI = self%pars(ialphaI )
      alphamu= self%pars(ialphamu)
      betamu = self%pars(ibetamu )

      mu0hat = tf*self%pars(imu0)
     &   *exp(alphamu * PMU + betamu * PMU**2)

      dmu0hatdl = mu0hat*(alphamu + 2D0 * betamu * PMU)

      d2mu0hatdl2=mu0hat * 2D0* betamu
     & + mu0hat*(alphamu + 2D0* betamu * PMU)**2

      V0hat  = ScaleTrait(PMU, tf*self%pars(iV0),  self%pars(ialphaV) )
      ! Initial slope of P-I curve
      aI     = ScaleTrait(PMU, self%pars(iaI0), alphaI             )
      Qs     = ScaleTrait(PMU, self%pars(iQ0 ), self%pars(ialphaQ))/2.
      RMchl  = tf*self%pars(iRMChl)
   
      ! Threshold irradiance and RMchl is set temperature dependent
      I_zero = zetaChl*RMchl/aI  
      ! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)
       IF( self%fV_opt .EQ. .true. ) THEN 
            if( par .gt. I_zero ) then
        
        larg1 = exp(1.0 + aI*par/(mu0hat*zetaChl))
        larg  = (1.0 + RMchl/mu0hat)*larg1   
        
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

       W_Y= WAPR(larg,0,0)

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
! First derivative of ThetaHat (has not been corrected for unimodal relationship of u0hat ~ PMU)
!       dThetaHatdl = (self%alphamu - self%alphaI)
!     &  *         (ThetaHat - (1/self%zetaChl)/(1 + WAPR(larg,0,0)))    
!     &  +        self%alphamu*( WAPR(larg,0,0)/(1 + WAPR(larg,0,0)))   
!     &  *         (mu0hat/(aI*self%par))*RMchl/(mu0hat + RMchl)   
! Second derivative of ThetaHat 
!         d2ThetaHatdl2 = (self%alphamu-self%alphaI) * dThetaHatdl
!     &  -          ( (aI*self%par/mu0hat)                               
!     &  *            ( (self%alphaI-self%alphamu) /self%zetaChl)**2 )   
!     &  *             WAPR(larg,0,0)/(1.+WAPR(larg,0,0))**3
!     &  -             2.*self%alphamu
!     &  *            ( (self%alphaI-self%alphamu)/self%zetaChl)       
!     &  *            ( RMchl/(RMchl + mu0hat) )                     
!     &  *            ( mu0hat/(aI*self%par) )
!     &  *              WAPR(larg, 0, 0)/(1 + WAPR(larg, 0, 0))**3
!     &  +              self%alphamu                            
!     &  *            ( (self%alphamu - self%alphaI))
!     &  -               self%alphamu
!     &  *            ( mu0hat/(RMchl + mu0hat) )  
!     &  *            (  RMchl/(RMchl + mu0hat) )
!     &  *            ( mu0hat / (aI*self%par) )  
!     &  *            WAPR(larg,0,0)/(1.+WAPR(larg, 0, 0))    
            else
               ThetaHat= 0.01  !  a small positive value 
               dThetaHatdl=0.
               d2ThetaHatdl2=0.
            endif
         ELSE
            ThetaHat = 0.8
            dThetaHatdl=0.
            d2ThetaHatdl2=0.
         ENDIF

        !The light limitation index (SI)
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
      muIhat   = max((mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat),
     &            errr) ! Net specific carbon uptake

       dmu0hatSIdl = mu0hatSI*(alphamu+2.*betamu*PMU) + mu0hat*dSIdl

       d2mu0hatSIdl2 = d2mu0hatdl2*SI+2.*dmu0hatdl*dSIdl+mu0hat*d2SIdl2

       dmuIhatdl = (1.- zetaChl*ThetaHat)
     & * dmu0hatSIdl  
     & - dThetaHatdl*zetaChl*(mu0hatSI+RMchl) 
   
       d2muIhatdl2=d2mu0hatSIdl2 - zetaChl 
     &  *(ThetaHat*d2mu0hatSIdl2+2.*dThetaHatdl*dmu0hatSIdl
     &  +mu0hatSI*d2ThetaHatdl2)
     &  -zetaChl*RMchl*d2ThetaHatdl2
    
      ! Define vNhat: the saturation function of ambient nutrient concentration
      select case(self%nutrient_uptake)  
      ! case 1: Classic Michaelis Menton 
        case(1)
      ! Potential maximal nutrient-limited uptake
          alphaK = self%pars(ialphaK) 
          alphaV = self%pars(ialphaV)
          ! Half-saturation constant of nitrate uptake
          Kn = ScaleTrait(PMU,self%pars(iK0), alphaK ) 
          VNhat = V0hat*NO3/(NO3 + Kn)
       dVNhatdl = -VNhat*alphaK*Kn/(NO3+Kn) + alphaV*VNhat

          d2VNhatdl2 = -alphaK*(Kn/(NO3+Kn)*dVNhatdl    
     &   - Kn*VNhat*alphaK* NO3/(NO3+alphaK)**2)
     &   + alphaV*dVNhatdl

      ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
        case(2)
         alphaA= self%pars(ialphaA)
         alphaV= self%pars(ialphaV)

          Lmin = log(0.2**3/6*pi) + log(1+Penfac)/( Penfac*alphaA) 

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
      end select  
      
      ! Define optimal fV
       IF( self%fV_opt .EQ. .TRUE. ) THEN
           ! Define ZINT
           alphaQ = self%pars(ialphaQ)
           zetaN  = self%pars(izetaN )
     
           ZINT = Qs*(muIhat/vNhat + zetaN)
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
            fV = (-1.0 + sqrt(1.0 + 1.0/ZINT))*Qs*muIhat/vNhat
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
            fV = 0.01
          dfVdl= 0.
         d2fVdl2 = 0.
         end if
   
         ELSE
           fV   = self%fV_initial 
           Qopt = 16.0/106.0 
           dfVdl= 0.
          d2fVdl2=0.
          dQoptdl=0.
          d2Qoptdl2=0.
       ENDIF
! Net growth rate (d-1) of phytoplankton at the average size as a function of temperature, light, and nutrient   
      self%muNet = max(muIhat*(1.0-fV-Qs/Qopt)
     &     - zetaN*fV*vNhat,errr)
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
          self%fV    = fV
          self%Theta = Theta
! Calculate the mean chl:C ratio of phytoplankton (averaged over all sizes) 
! using the second derivative of the cell quota w.r.t. log size. (Why negative?)
            self%ThetaAvg = Theta + VAR/2d0          
     &    *  ( d2ThetaHatdl2*(1.0 - fV - Qs/Qopt)         
     &    -      ThetaHat   *(d2fVdl2 + Qs*((dQoptdl**2)*(2./Qopt)
     &    - d2Qoptdl2)/(Qopt**2)))
            
            if (self%ThetaAvg .lt. 0) self%ThetaAvg = Theta
   
! Calculate the mean N:C ratio of phytoplankton (averaged over all sizes) 
! using the second derivative of the cell quota w.r.t. log size. (How to derive?)
       self%QNavg = Qopt/(1.+((2/Qopt)*dQoptdl**2 - d2Qoptdl2)
     &    * VAR/(2*Qopt))  
          
! Total Chl a (unit: mg Chl m-3)
          self%env(iCHL) = self%env(iPHY)
     &    /self%QNavg*self%ThetaAvg

! Calculate the community sinking rate of phytoplankton
      ! Phytoplankton sinking rate at the average size
      w_p = ScaleTrait(PMU,self%pars(iw_p),self%pars(ialphaW))
      d2wpdl2 = self%pars(ialphaW)**2*w_p
      self%w_pAvg = w_p+0.5*VAR*d2wpdl2 ! Sinking rate of phytoplankton

      end subroutine growth
!-----------------------------------------------------------------------
   
!---------------------------------------------------------
! The -1 branch of the Lambert-W function
      PURE REAL FUNCTION w1(x)
         IMPLICIT NONE
         REAL(8), PARAMETER  :: m1=0.3361, m2=-0.0042, m3=-0.0201
         REAL(8), INTENT(IN) :: x
         REAL(8) :: s
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
!-------------------------------------------------------------------
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
       integer, parameter :: N=10**6
       real :: x(N), Ptot=0.
       real, parameter :: lmin = log(pi/6D0*0.5**3),
     &    lmax =  log(pi/6D0*500**3)

          tf  = self%tf()
          PHY = self%env(iPHY)
          PMU = log(self%env(iPMUPHY)/self%env(iPHY)) ! Convert the unit to log(mu**3)
          PMU = min(max(PMU,PMUmin),PMUmax)
          VAR = self%env(iVARPHY)/PHY
!        ! Add a factor that makes ghat decrease with increasing phyto size
!         rhol = exp(-self%pars(iS)*PMU)
!        ! Calculate the integration  
!         x(1) = max(lmin,PMU-6D0*sqrt(VAR))  ! Minimal size
!         x(N) = min(lmax,PMU+6D0*sqrt(VAR))  ! Maximal size
!         dx   = (x(N) - x(1))/(N-1)
!         do i=2,N-1
!           x(i) = x(i-1) + dx 
!         enddo
!        ! Calculate total phytoplankton palatable prey:
!         do i=1,N-1
!           Ptot = Ptot + normal(x(i))*dx*exp(-self%pars(iS)*x(i))
!         enddo
         
         

          kp  = self%pars(ikp)

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
            ghat = 1.0 - exp(-LOG(2D0)*PHY/kp)  
          END SELECT

    
        ! Calculate feeding threshold of ghat (dimensionless)
          g_th = self%pars(iCf)/self%pars(iEmax)/(1.0 - self%pars(iCa))
   
        ! Calculate total foraging activity (d-1)
          aAt = tf*self%pars(igmax)/self%pars(ibeta)*      
     &      (-1.0 - w1(-(1.0-g_th)/EXP(1.0+self%pars(ibeta))))
   
        ! Compute the optimal foraging activity to maximize net growth of zooplankton
          IF (ghat .GT. g_th) THEN
            if (self%CF_opt .eq. .true.) then
           ! Optimal foraging, aAf: Specific foraging activity (d-1) = 1/h (h: biomass specific handling time)
           !With increasing food conc., optimal aAf slightly increases with food conc.
              aAf = aAt*self%pars(ibeta)/(-1.0 - w1(-(1.0-g_th/ghat) 
     &         /EXP(1.0 + self%pars(ibeta))))
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
            AE  = self%pars(iEmax)*
     &     (1.0 - exp(-self%pars(ibeta)*(aAt/aAf - 1.0)))  
             
            ELSE
   
            aAf = 0.0
            AE  = self%pars(iEmax)
          ENDIF
        
        !Ingestion rate (d-1) of zooplankton (gmax*Q(PT)), now should be a function of size 
          self%INGES = aAf * ghat 
          
          ! Zooplankton respiration (-> Nut)
          self%RES = self%pars(iRMzoo)*tf
     & +             self%pars(iCa)*AE*self%INGES
     & +             self%pars(iCf)*aAf    
   
! Grazing rate on the mean size (specific to N-based Phy biomass, unit: d-1) (Eq. 17)
          self%gbar = self%INGES*self%env(iZOO)
     & /             PHY*sqrt(self%pars(ialphaG))

          ! Zooplankton net growth (To force it positive)
          self%NGR = max(AE*self%INGES - self%RES,errr)
   
          ! ZOOPLANKTON EGESTION (-> Detritus)
          self%EGES = max(self%INGES - self%NGR - self%RES,errr)

         
          self%dgdlbar = 0.  ! for all cases (when evaluated at the mean size; Eq. B-3)
          
          self%d2gdl2bar = -(self%pars(ialphaG)-1.d0)*
     &         self%gbar/VAR  ! Eq. B-6 (always negative)   

      end subroutine Zooplankton
!--------------------------------------------------------------------------------------
      real function normal(mean,var,l)
        implicit none
        real, intent(in) :: mean,var,l
        normal = 1D0/sqrt(2D0*var*pi)*exp(-(l-mean)**2/var/2D0) 
      end function
!--------------------------------------------------------------------------------------
      subroutine ODE(self)
       implicit none
       class(type_PZEFT),intent(inout),target :: self
       real I_min,rpd,tf
       integer, parameter :: numc=4
       integer :: i,j
       real :: PP_nz,PP_nd,PP_dp,PP_dz,PP_pn,PP_zp,PP_np,
     &         PHY,ZOO,DET,VAR,PMU,dPMUdt,dVARdt

      tf = self%tf()

      PHY= self%env(iPHY)
      PHY= min(max(PHY,PHYmin),PHYmax)
      VAR= self%env(iVARPHY)/PHY  !self%env(iVARPHY) is the product of PHY and VAR
      VAR= min(max(VAR,VARmin),VARmax)
      PMU= log(self%env(iPMUPHY)/PHY)
      PMU= min(max(PMU,PMUmin),PMUmax)
      DET= self%env(iDET)
      DET= min(max(DET,DETmin),DETmax)
      ZOO= self%env(iZOO)
      ZOO= min(max(ZOO,ZOOmin),ZOOmax)
      
      ! Define I_min as 1% of I_0 (bottom of euphotic layer)
      I_min = 0.01*self%env(iI_0)
   
      ! Loss rate of phytoplankton to detritus depends on local light intensity.
      if (self%env(ipar) .ge. I_min) then
         rpd = tf*self%pars(irpdu)
      else
         rpd = tf*self%pars(irpdl)
      endif
   
      self%dPHYdt = 0.5*PHY*(2D0*(self%muNet - rpd) 
     &  + VAR*self%d2muNetdl2)
     &  - ZOO*self%INGES  ! The total amount of phytoplankton grazed by zoo is independent of average size    

      ! Nutrient changes (Eq. 23)
      self%dNO3dt = self%RES*ZOO + self%pars(irdn)*DET*tf
     &  - PHY*(self%muNet + 0.5*VAR*self%d2muNetdl2)

! dz = zooplankton growth(excluding excretion and egestion) - zooplankton mortality (molN d-1)
      self%dZOOdt = ZOO*self%NGR - ZOO**2*self%pars(imz)*tf
! dD = phytoplankton mortality + zooplankton mortality + zooplankton egestion - conversion-to-nutrient (molN d-1)
      self%dDETdt = -(self%dNO3dt+self%dPHYdt+self%dZOOdt)   

! Since we treat the product of (PMU*PHY) as a single unit in the differential eq.
! d(PMUPHY)dt = PMU*dPHYdt + PHY*dPMUdt (here the PMU is exp(PMU)), the same applies to d(VARPHY)dt
      dPMUdt = VAR*(self%dmuNetdl - self%dgdlbar)
      dVARdt = VAR**2*(self%d2muNetdl2 - self%d2gdl2bar)  !Eq. 19 
!      if (self%d2muNetdl2 > 0) write(6,*) 'd2muNetdl2 > 0!'
      self%dPMUPHYdt=exp(PMU)*(self%dPHYdt + PHY*dPMUdt)
      self%dVARPHYdt=VAR*self%dPHYdt + PHY*dVARdt

! For PPDD matrix:
      PP_nz = self%RES*ZOO
      PP_nd = self%pars(irdn)*DET*tf
      PP_dp = PHY*rpd
      PP_dz = ZOO**2*self%pars(imz)*tf + ZOO*self%EGES
      PP_zp = ZOO*self%INGES
      PP_pn = PHY*(self%muNet + 0.5*VAR*self%d2muNetdl2)
      PP_np = 0.
      if (PP_pn .lt. 0) then
!        write(6,*) 'PP_pn < 0!'
        PP_pn = PHY*self%muNet 
        PP_np = -PHY*0.5*VAR*self%d2muNetdl2
      endif
! Write out production/destruction matrix
      ! 1: NO3; 2: PHY; 3: ZOO; 4: DET
      self%pp      = 0.
      self%dd      = 0.
      self%pp(1,2) = PP_np
      self%pp(1,3) = PP_nz
      self%pp(1,4) = PP_nd
      self%pp(2,1) = PP_pn
      self%pp(3,2) = PP_zp
      self%pp(4,2) = PP_dp
      self%pp(4,3) = PP_dz
      
      do i = 1,numc
        do j = 1,numc
          self%dd(j,i)=self%pp(i,j) 
        end do
      end do 

      end subroutine ODE
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------
      subroutine modified_patankar_2(self)
!
! !DESCRIPTION:
! Here, the second-order Modified Patankar-Runge-Kutta (MPRK2) scheme is coded,
! with two evaluations of the right hand sides per time step:  (Burchard et al. 2003)
!
! \begin{equation}\label{eq:am:MPRK}
!   \left. \begin{array}{rcl}
!     c_i^{(1)} &=&
! \displaystyle
! c_i^n  + \Delta t
! \left\{
! \sum\limits_{\stackrel{j=1}{j \not= i}}^I p_{i,j}\left(\underline{c}^n\right)
! \dfrac{c_j^{(1)}}{c_j^n}
! + p_{i,i}\left(\underline{c}^n\right)
! - \sum_{j=1}^I d_{i,j}\left(\underline{c}^n\right)
! \dfrac{c_i^{(1)}}{c_i^n}
! \right\},
! \\ \\
! c_i^{n+1} &=&
! \displaystyle
! c_i^n  + \dfrac{\Delta t}{2}
!                   \left\{
!                     \sum\limits_{\stackrel{j=1}{j \not= i}}^I
!                       \left(p_{i,j}\left(\underline{c}^n\right)
!                           + p_{i,j}\left(\underline{c}^{(1)}\right)
!                       \right) \dfrac{c_j^{n+1}}{c_j^{(1)}}
!                       + p_{i,i}\left(\underline{c}^n\right)
!                       + p_{i,i}\left(\underline{c}^{(1)}\right)
! \right.\\ \\
!               & &
! \displaystyle
! \left.\phantom{c_i^n  + \dfrac{\Delta t}{2} }
!                   - \sum_{j=1}^I
!                       \left(d_{i,j}\left(\underline{c}^n\right)
!                           + d_{i,j}\left(\underline{c}^{(1)}\right)
!                       \right) \dfrac{c_i^{n+1}}{c_i^{(1)}}
!                   \right\}.
!   \end{array}
!   \right\}
! \end{equation}
!
! !USES:
        IMPLICIT NONE
!
! !INPUT PARAMETERS:
       integer, parameter  :: numc=4,stdout=6
       class(type_PZEFT),intent(inout),target :: self
!
! !INPUT/OUTPUT PARAMETER:
        REAL :: cc(1:numc)
     
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
       REAL :: pp2(1:numc,1:numc),dd2(1:numc,1:numc)
       REAL :: pp1(1:numc,1:numc),dd1(1:numc,1:numc)
       REAL :: a(1:numc,1:numc),r(1:numc)
       REAL :: cc1(1:numc)
       real :: dPMUPHYdt1,dVARPHYdt1,PMUPHY,VARPHY
       REAL, parameter :: eps = 1D-8
       integer  :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
! Filter out NaN  values and restore them to a small positive number!
      if(self%env(iPMUPHY) /= self%env(iPMUPHY)) self%env(iPMUPHY)=eps 
      if(self%env(iVARPHY) /= self%env(iVARPHY)) self%env(iVARPHY)=eps
      if(self%env(iNO3) /= self%env(iNO3)) self%env(iNO3) =eps
      if(self%env(iPHY) /= self%env(iPHY)) self%env(iPHY) =eps
      if(self%env(iZOO) /= self%env(iZOO)) self%env(iZOO) =eps
      if(self%env(iDET) /= self%env(iDET)) self%env(iDET) =eps 

        cc(1)=self%env(iNO3)
        cc(2)=self%env(iPHY)
        cc(3)=self%env(iZOO)
        cc(4)=self%env(iDET)
        PMUPHY = self%env(iPMUPHY)
        VARPHY = self%env(iVARPHY)
!  absolutely essential since not all elements are calculated
        pp2=0. ; pp1=0.
        dd2=0. ; dd1=0.
     
        call self%initialize
        call self%growth
        call self%Zooplankton
        call self%ODE
        
         
        pp1=self%pp
        dd1=self%dd
        dPMUPHYdt1 = self%dPMUPHYdt    
        dVARPHYdt1 = self%dVARPHYdt

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
      
      self%env(iNO3)=cc1(1)
      self%env(iPHY)=cc1(2)
      self%env(iZOO)=cc1(3)
      self%env(iDET)=cc1(4) 
      self%env(iPMUPHY)=self%env(iPMUPHY) + self%dPMUPHYdt
      self%env(iVARPHY)=self%env(iVARPHY) + self%dVARPHYdt

        call self%growth
        call self%Zooplankton
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
     
      self%env(iNO3)=cc(1)
      self%env(iPHY)=cc(2)
      self%env(iZOO)=cc(3)
      self%env(iDET)=cc(4) 
      self%env(iPMUPHY)=PMUPHY 
     & + 0.5*(self%dPMUPHYdt + dPMUPHYdt1)
      self%env(iVARPHY)=VARPHY
     & + 0.5*(self%dVARPHYdt + dVARPHYdt1)

      end subroutine modified_patankar_2
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Matrix solver
!
! !INTERFACE:
        subroutine matrix(n,a,r,c)
!
! !DESCRIPTION:
! This is a Gaussian solver for multi-dimensional linear equations.
!
! !USES:
        IMPLICIT NONE
!
! !INPUT PARAMETERS:
        integer, intent(in)                 :: n
!
! INPUT/OUTPUT PARAMETERS:
       REAL                             :: a(1:n,1:n),r(1:n)
!
! OUTPUT PARAMETERS:
       REAL, intent(out)                :: c(1:n)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
       integer  :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
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
