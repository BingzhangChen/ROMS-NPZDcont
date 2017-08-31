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
      MODULE PZEFTFe
      use lambert
      implicit none
      private
      integer,parameter,public :: NVAR = 4,
     &                            iNO3 = 1,
     &                            iPHY = 2,
     &                            iZOO = 3,
     &                            iDET = 4
! PUBLIC DERIVED TYPES:
      Type     :: PHY
         !State variables:
         real  :: Nconc,LNV,VAR,
         !Input parameters:
     &            mu0,aI0,alphaI,V0N,alphaV,K0N,alphaK,W0,alphaW,
     &            Ep,alphamu,betamu,Q0N,alphaQ,Penfac,zetaN,RMchl,
     &            A0N,alphaA,K0Fe,tf,LNVmax,VARmax,zetaChl,Lref,alphaFe
       integer :: nutrient_uptake
         !Output:
         real  :: ThetaAvg,QNavg,muNet,w_pAvg,fV,QN,muIhat,VNhat,
     &            dmuNetdl,d2muNetdl2,dwdl,d2wdl2,
     &            dlnVdt,dVARdt,dPHY

      end type PHY

      Type     :: ZOO
          real :: Nconc,tf,Ez,
     &            gmax,kp,alphaG,ga,gb,gc,
     &            beta,Ca,Cf,Emax,RMzoo,
     &            GGE,unass,mz
       integer :: grazing_formula
          !Output:
          real :: gbar,dgdlbar,d2gdl2bar,dZOO
     &            AE,INGES,RES,NGR,EGES
      end Type ZOO

      Type     :: NO3
          real :: Nconc,dNO3
      end Type NO3

      Type     :: DET
          real :: Nconc,rDN,dDET
      end Type DET

      Type     :: IRON
          real :: conc,lFe,KScm,Ksc,Fe_N,Femin
      !Output:
     &           ,Fe_scav,dFe
      end Type IRON

      Type     :: ENVR
          real :: temp,dtdays,light
      end Type ENVR

      Type,public :: NPZDFe 
        type(PHY) :: p      !A continous size class of phytoplankton
        type(NO3) :: n
        type(ZOO) :: z      !Two size classes of zooplankton
        type(DET) :: D
        type(IRON):: Fe
        type(ENVR):: env
        real      :: pp(NVAR,NVAR),vars(NVAR)      ! For PPDD
!-------------------------------------------------------------------------------
      contains
         procedure :: initialize
         procedure :: phyto
         procedure :: IRONCYCLE
         procedure :: Zooplankton
         procedure :: ODE
         procedure :: modified_patankar_2
      END TYPE NPZDFe

      real, parameter :: eps  = 1d-9, 
     &                   PI   = 3.141592653529D0

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
      class (NPZDFE), intent(inout), target :: self 
      real :: dtdays,tC
      real :: w31_p, w31_d,E_p,E_z
!Stokes law: V = 2/9*(rho_p - rho_f)*g*R**2/mu (mu: dynamic viscosity, which is function
!of temperature.) w_p0 = w(31)*(0.451 + 0.0178*T) (Gregg et al. 2008). 
!We also set this to be size dependent
!such that w_p = w_p0* Wstar * exp(logsize * alphaW) (Acevedo-Trejos et al. 2015)
!This is quite critical in modeling size and size variance!

      w31_p            = 0.02   
      w31_d            = 5.0 
      dtdays           = self%env%dtdays 
      tC               = self%env%temp
      self%p%mu0       = 
#if defined NPZDSIMPLE || defined NPZDSIZE
     &                   1.5  *dtdays
#else
     &                   5.0  *dtdays      !Smith et al. 2015
#endif
! Note:
! Initial slope of the P-I curve (mol C gChl-1 (W m-2)-1 d-1). 
! 1 W m-2 = 4.6 uE m-2 s-1 = 1.84 E m-2 d-1
! Pahlow et al. (2013) estimated the aI0 0.5~3.7 m2 mol C gChl-1 E-1
! corresponding to 0.27~2.01 mol C gChl-1 (W m-2)-1 d-1.
! Ward et al. (L&O) used the size independent slope of 3.83E-7 mol C (gChl-1) (uE m-2)-1. 
! Edwards et al. L&O (2015) estimated a size dependent value of 10**(-1.36) d-1 (uE m-2 s-1)-1 (= 5.05E-7 (uE m-2)-1)
! So I choose a value of 10**(-1.36)*4.6/(50/12)= 0.05 
      self%p%aI0       = 0.25 *dtdays      ! Fennel et al. (2006)
      self%p%V0N       = 5d0  *dtdays      ! Smith et al. 2015
      self%p%K0N       = 
#if defined NPZDSIMPLE || defined NPZDSIZE
     &           0.02    !Half saturation constant for growth
#else
     &           0.17    !Half saturation constant for nitrogen uptake (Litchman et al. 2007)
#endif
      !sinking rate (m) of phytoplankton, POSITIVE
      self%p%W0        = w31_p*dtdays
      !Sinking rate of POM
      self%D%wDET      = w31_d*dtdays
      !Maximal biomass-specific ingestion rate of zooplankton at 15 ºC
      self%z%gmax      = 
#ifdef OCF
     &                   2.5*dtdays   !Pahlow et al. MEPS (2010)
#else 
     &                   1.0*dtdays   !Fennel et al. GBC (2006)
#endif

      !Chai et al. (2002),half-saturation constant for zooplankton grazing
      self%z%kp        = 1.0
      !The conversion rate from DET to DIN at reference temperature (Ward et al. 2012)
      self%D%rDN       = 0.05*dtdays
      !Phytoplankton activation energy (eV)
      self%p%Ep        = 0.4          
      !Zooplankton activation energy (eV)
      self%z%Ez        = 0.6 
      !Coefficient for zooplankton mortality at 15 C
      self%z%mz        = 0.05*dtdays  !Chai et al. (2002) 
      self%p%alphamu   = 0.0
      self%p%betamu    = 0.0
      self%p%alphaI    = 0.0     ! Smith et al. 2015,but -0.13 Edwards et al. L&O 2015
      self%p%A0N       = 5.0*dtdays     ! From Smith et al. JPR 2015
      self%p%alphaA    = -0.3
      self%p%alphaV    = 0.0  

!Normalization constant for half-saturation constant for nitrate (mmolN m-3) (Litchman et al. 2007)
      self%p%alphaK    = 0.27    
      self%z%alphaG    = 1.2
      self%p%alphaW    = 0.39      !Smayda 1970; Laws 1975
      ! Assume that ZOO grazing does not affect size
      self%z%ga        = 0.0
      self%z%gb        = 0.0
      self%z%gc        = 0.0
      self%p%K0Fe      = 0.08      !Unit: nM. Gregg (2003) used ~0.1 nM for phyto
      self%p%alphaFe   = 0.14      !Gregg 2003
      self%p%Q0N       = 0.1       ! 0.07 in Litchman et al. 2007; Ward et al. L&O 2012 
      self%p%Penfac    = 1.0
      self%p%alphaQ    = -0.17     ! -0.17 in Litchman et al. 2007; Ward et al. L&O 2012 
      self%p%Lref      = 0.5       ! Reference size from which growth rate starts to decline
      self%p%zetaN     = 0.6       ! Smith et al. 2015;Cost of DIN uptake (mol C/mol N)   
      self%p%RMchl     = 0.1*dtdays!Cost of chl maintenance (d-1) 
      self%p%zetaChl   = 0.8       !Cost of photosynthesis coefficient (mol C/g Chl)
      self%z%beta      = 0.2             !Assimilation coefficient 
      self%z%Ca        = 0.3   !Cost of assimilation coefficient
      self%z%Cf        = 0.3   !Effective cost of foraging coefficient  
      self%z%Emax      = 0.99            !Maximal assimilation efficiency
      self%z%RMzoo     = 0.16*dtdays !Specific maintenance respiration (d-1)
      self%z%GGE       = 0.33
      self%z%unass     = 0.0
      self%Fe%KScm     = 3D-5*dtdays  !Minimal scavenging rate of Fe (d-1)
      self%Fe%Ksc      = 3D-2*dtdays  !Particle dependent scavening rate of Fe (umolN-1 d-1)
      self%Fe%lFe      = 0.6          ! Iron ligand concentration (nM). (TOM10 P. 19)
      self%Fe%Femin    = 0.02         !minimal Fe conc. 0.02 nM (Aumont et al. 2003) 
      self%Fe%Fe_N     = 0.0265       !Fe:Nitrogen molar ratio
! End of pars definition.

! Choose Different models:
      self%p%nutrient_uptake= 2      !The model used for nutrient uptake
      self%z%grazing_formula= 3      !Choose different grazing functional groups

!Calculate temperature factors for phytoplankton and others:
      E_p         = self%p%Ep
      E_z         = self%z%Ez
      self%p%tf   = temp(E_p,tC)
      self%z%tf   = temp(E_z,tC)
      self%p%LNVmax = 20.
      self%p%VARmax = 100.

      end subroutine initialize
   
      
!---------------------------------------------------------
      pure real function temp(Ea,tC)
!IROUTINE: Temperature dependence of plankton rates. Added by Bingzhang Chen on 20150213
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
       real, intent (in)  :: Ea, tC
       real, parameter    :: kb = 8.62d-5, ! boltzman constant constant [ eV /K ]
     &                       Tr = 15.
      
      temp = exp(-(Ea/kb)*(1D0/(273.15 + tC)-1D0/(273.15 + Tr)))
      
      return 
      end function temp
!EOC
 
!-----------------------------------------------------------------------
!==========================================
! !IROUTINE:Computing phytoplankton growth components 
!  and first and second derivatives at the mean size  
!
! !INTERFACE:
      subroutine phyto(self)
      implicit none
!INPUT PARAMETERS:
      class (NPZDFE),intent(inout) :: self
!LOCAL VARIABLES:
      real :: I_zero, tf, PMU,VAR,par,NO3,zetaChl,alphaI,alphamu,
     &    betamu,PHY,X,B,dXdl,dBdl,alphaFe,KFe,Fe,dmu0hatdl, d2mu0hatdl2
      real :: V0hat,Kn,Lmin,A0hat,dA0hatdl,d2A0hatdl2,fA
      real :: VNhat,dVNhatdl,d2VNhatdl2, fN, d2fNdl2 ! Nutrient uptake variables
      real :: mu0hat,muIhat,mu0hatSI,dmu0hatSIdl,d2mu0hatSIdl2
      real :: dmuIhatdl,d2muIhatdl2! Growth rate variables
      real :: fV,dfVdl,d2fVdl2  !fV
      real :: ZINT,dZINdl,d2ZINdl2 !ZINT
      real :: aI,SI,dSIdl,d2SIdl2         !Light dependent component
      real :: RMchl,Theta,ThetaHat,dThetaHatdl,d2ThetaHatdl2 !Chl related variables
      real :: QN,Qs,dQNdl,d2QNdl2  ! cell quota related variables
      real    :: larg1,larg,w_p,dmu0hat_aIdl,d2mu0hat_aIdl2,dlargdl,
     &        d2largdl2,W_Y,dWYYdl,daI_mu0hatdl,d2aI_mu0hatdl2,
     &        alphaK,alphaV,alphaA,alphaQ,zetaN,d2wpdl2
      real :: Penfac,Lref,dV0hatdl,d2V0hatdl2
      real,parameter :: Wmax = 50. !Maximal sinking rate
     &                 ,thetamin = 0.02
     &                 ,QNmin = 0.01
   
!-----------------------------------------------------------------------
! Phytoplankton codes below mostly come from S. Lan Smith
! calculate the trait values at the mean size, lmean 
! multiplying RATES only by tC (Temperature dependence)

! Only for pot. max. affinity, apply a penalty function (to avoid unrealistic increase at v. small cell size). 
! To get the max. affinity at ESD = 0.2 micrometers: 
! The maximal growth rate (mu0hat) under temperature tC 
      par = max(self%env%light,eps)
      NO3 = max(self%N%Nconc,eps)
      PHY = max(self%p%Nconc,eps)
      tf  = self%p%tf
! Reset the values to be calculated:
      
      alphaQ  = self%p%alphaQ
      alphaI  = self%p%alphaI 
      alphamu = self%p%alphamu
      betamu  = self%p%betamu
      Penfac  = self%p%Penfac      
      PMU     = self%p%lnV
      VAR     = self%p%VAR

#ifdef IRON
      ! Fe related:
      Fe      = max(self%Fe%conc,self%Fe%Femin)  !Dissolved Fe concentration
      alphaFe = self%p%alphaFe
      KFe     = ScaleTrait(PMU, self%p%K0Fe,alphaFe) !The half saturation constant of Fe at average size
#endif

      Qs      = ScaleTrait(PMU, self%p%Q0N, alphaQ)/2.
    
      mu0hat = tf*self%p%mu0*exp(alphamu*PMU + betamu*PMU*PMU)
#ifdef IRON
     & * Fe/(Fe + KFe)
    
      X =alphaFe*KFe/(Fe + KFe) 
#endif
      dmu0hatdl  = mu0hat*(alphamu + 2D0*betamu*PMU 
#ifdef IRON
     &  - X 
#endif
     &  )

      d2mu0hatdl2= dmu0hatdl*(alphamu + 2D0*betamu*PMU
#ifdef IRON
     &  - X
#endif
     &  ) + mu0hat * (2D0*betamu
#ifdef IRON
     &  - alphaFe*Fe*X/(Fe+KFe)
#endif
     &  )

      zetaChl= self%p%zetaChl
      zetaN  = self%p%zetaN
      alphaV = self%p%alphaV

      ! Iron limits nitrogen uptake
      V0hat  = ScaleTrait(PMU, tf*self%p%V0N, alphaV)
#ifdef IRON
     &  *Fe/(Fe + KFe)
#endif

      dV0hatdl   = V0hat*(alphaV
#ifdef IRON
     &  - X
#endif
     &  )
      d2V0hatdl2 = dV0hatdl*(alphaV 
#ifdef IRON
     &  - X
#endif
     &  ) - V0hat
#ifdef IRON
     &  * alphaFe*Fe*X/(Fe+KFe) 
#endif

      ! Initial slope of P-I curve
      aI     = ScaleTrait(PMU, self%p%aI0, alphaI)
      ! Cost of photosynthesis
      RMchl  = tf*self%p%RMchl
      ! Threshold irradiance and RMchl is set temperature dependent
      I_zero = zetaChl*RMchl/aI  

      !Define VNhat: the saturation function of ambient nutrient concentration
       selectcase(self%p%nutrient_uptake)  
      ! case 1: Classic Michaelis Menton 
        case(1)
      ! Potential maximal nutrient-limited uptake
          alphaK = self%p%alphaK 
          ! Half-saturation constant of nitrate uptake
          Kn     = ScaleTrait(PMU,self%p%K0N, alphaK) 
          VNhat  = V0hat*NO3/(NO3 + Kn)

       dVNhatdl  = -VNhat*alphaK*Kn/(NO3+Kn) + NO3/(NO3 + Kn)*dV0hatdl

          d2VNhatdl2 = -alphaK*(VNhat * alphaK * NO3/(NO3+Kn)**2*Kn
     &   + Kn/(NO3 + Kn)*dVNhatdl)
     &   + NO3/(NO3+ Kn)*d2V0hatdl2
     &   - dV0hatdl*NO3/(NO3 + Kn)**2*alphaK*Kn

      ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
        case(2)

         alphaA= self%p%alphaA
          Lref = self%p%Lref
          Lmin = log(Lref**3/6.*pi) + log(1.+Penfac)/( Penfac*alphaA) 

          A0hat = PenAff(PMU, alphaA, Penfac, Lmin )
     &        * ScaleTrait(PMU, tf*self%p%A0N, alphaA)

          A0hat = max(A0hat,eps*self%env%dtdays)

       dA0hatdl = alphaA*A0hat 
     &        - self%p%A0N*exp(PMU*alphaA)*Penfac*alphaA
     &        * exp(Penfac*alphaA *(PMU-Lmin))

       d2A0hatdl2 = alphaA*dA0hatdl 
     &        - Penfac*alphaA*exp(alphaA*((1.+Penfac)*PMU-Penfac*Lmin)) 
     &        * (dA0hatdl + A0hat*alphaA*(1.+Penfac))  
   
          !Define fA
          fA = 1.0/( 1.0 + sqrt(A0hat * NO3/V0hat) ) 

       VNhat = (1.0-fA)*V0hat*fA*A0hat*NO3/
     &        ((1.0-fA)*V0hat + fA*A0hat*NO3) 
   
       !X: temporary variable
       X    = V0hat/A0hat + 2.*sqrt(V0hat*NO3/A0hat) + NO3       

       !B: d(V0/A0)dl
       B    = dV0hatdl/A0hat - V0hat/A0hat**2*dA0hatdl

       dXdl = B*(1. + sqrt(NO3 * A0hat / V0hat))

       dBdl = d2V0hatdl2/A0hat - dV0hatdl*dA0hatdl/A0hat**2
     &   - (V0hat/A0hat**2*d2A0hatdl2
     &   +  dA0hatdl*(dV0hatdl/A0hat**2 - 2.*V0hat*dA0hatdl/A0hat**3))
       
       dVNhatdl = NO3*(dV0hatdl/X 
     &   - V0hat/X**2*B*( 1.+ sqrt(NO3*A0hat/V0hat) ) )

       d2VNhatdl2 = NO3*(d2V0hatdl2/X - dV0hatdl*dXdl/X**2
     &   - (V0hat/X**2*B*(-sqrt(NO3)*.5*(A0hat/V0hat)
     &   * sqrt(A0hat/V0hat) * B)
     &   + B*(1. + sqrt(NO3*A0hat/V0hat) )
     &   * (dV0hatdl/X**2 - 2.*V0hat*dXdl/X**3)
     &   + V0hat/X**2*(1.+sqrt(NO3*A0hat/V0hat))*dBdl  )  )

       case default
        write(*,*) 'Error: Incorrect option for nutrient uptake!'
      ENDSELECT  
! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)
! Only calculate within the euphotic zone, otherwise many numerical problems.
      if( par .gt. I_zero ) then
        
        larg1 = exp(1. + min(aI*par/(mu0hat*zetaChl),600.))

        larg  = (1. + RMchl/mu0hat)*larg1   
        
        dmu0hat_aIdl   = (dmu0hatdl - alphaI*mu0hat)/aI

        d2mu0hat_aIdl2 = d2mu0hatdl2/aI -alphaI/aI*dmu0hatdl 
     &  - aI*dmu0hat_aIdl

        daI_mu0hatdl = -(aI/mu0hat)**2*dmu0hat_aIdl

        d2aI_mu0hatdl2 = -((aI/mu0hat)**2*d2mu0hat_aIdl2
     &  - 2./(mu0hat/aI)**3*(dmu0hat_aIdl)**2)

        dlargdl = -RMchl*larg1/(mu0hat**2)*dmu0hatdl
     &  + (1.+RMchl/mu0hat)*larg1 * par/zetaChl*daI_mu0hatdl
        
        d2largdl2 = -RMchl*(larg1*mu0hat**(-2)*d2mu0hatdl2
     &  + larg1*par/zetaChl*daI_mu0hatdl*mu0hat**(-2)*dmu0hatdl
     &  + larg1*dmu0hatdl*(-2.*mu0hat**(-3)*dmu0hatdl))      ! Corrected on Dec. 15 2015
     &  + par/zetaChl*((1+RMchl/mu0hat)*larg1*d2aI_mu0hatdl2
     &  + (1.+RMchl/mu0hat)*larg1*par/zetaChl*daI_mu0hatdl*daI_mu0hatdl
     &  + RMchl*(-mu0hat**(-2)*dmu0hatdl)*larg1*daI_mu0hatdl)

       W_Y      = WAPR(larg,0,0)
       ThetaHat = 1.0/zetaChl + (1.0- W_Y)*mu0hat/(aI * par)
       ThetaHat = max(ThetaHat,thetamin) 

       dThetaHatdl = 1./par
     &  *  (-W_Y/larg/(1.+W_Y)*dlargdl*mu0hat/aI
     &  +  (1.-W_Y)*dmu0hat_aIdl)
       
       dWYYdl = dlargdl*(-W_Y**2/larg**2/(1.+W_Y)**3
     &  -  W_Y/larg**2/(1.+W_Y) + W_Y/(larg*(1.+W_Y))**2)

       d2ThetaHatdl2 = 1/par*(-(W_Y/larg/(1.+W_Y)*dlargdl*dmu0hat_aIdl
     &  +  W_Y/larg/(1.+W_Y)*d2largdl2*mu0hat/aI
     &  +  dWYYdl*dlargdl*mu0hat/aI)
     &  -  W_Y/larg/(1.+W_Y)*dlargdl * dmu0hat_aIdl
     &  +  (1.-W_Y)*d2mu0hat_aIdl2)

        SI = 1. - max(exp(-aI*par*ThetaHat/mu0hat),0.)

        dSIdl = ( (alphaI- dmu0hatdl/mu0hat)
     &   * ThetaHat + dThetaHatdl) * (1.-SI)*aI*par/mu0hat    !confirmed

       d2SIdl2 = par*(- dSIdl*aI/mu0hat*(ThetaHat*alphaI 
     & - ThetaHat/mu0hat*dmu0hatdl + dThetaHatdl) + (1.-SI)
     & * (ThetaHat*alphaI- ThetaHat/mu0hat*dmu0hatdl + dThetaHatdl)
     & * daI_mu0hatdl + (1.-SI)*aI/mu0hat*(
     & - (d2mu0hatdl2/mu0hat - dmu0hatdl**2/mu0hat**2)*ThetaHat
     & + (alphaI-dmu0hatdl/mu0hat)*dThetaHatdl + d2ThetaHatdl2)  )

       ! Light dependent growth rate 
       ! (needs to take into account the cost of dark and light-dependent chl maintenance)
       mu0hatSI = mu0hat*SI  ! Gross specific carbon uptake (photosynthesis)
       muIhat   = mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat ! Net specific carbon uptake
       muIhat   = max(muIhat,eps*self%env%dtdays)

       dmu0hatSIdl = SI*dmu0hatdl + mu0hat*dSIdl

       d2mu0hatSIdl2 = d2mu0hatdl2*SI+2.*dmu0hatdl*dSIdl+mu0hat*d2SIdl2 !Correct

       dmuIhatdl = (1D0- zetaChl*ThetaHat)*dmu0hatSIdl  
     & - dThetaHatdl*zetaChl*(mu0hatSI+RMchl) !Correct

       d2muIhatdl2=d2mu0hatSIdl2 - zetaChl 
     &  *(ThetaHat*d2mu0hatSIdl2+2.*dThetaHatdl*dmu0hatSIdl
     &  +mu0hatSI*d2ThetaHatdl2)
     &  -zetaChl*RMchl*d2ThetaHatdl2   !Correct

       ZINT   = Qs*(muIhat/VNhat + zetaN)

       dZINdl = Qs*(dmuIhatdl/VNhat - muIhat*dVNhatdl/VNhat**2)
     &           +  alphaQ*ZINT    
   
   
       d2ZINdl2 = Qs/VNhat*((alphaQ-dVNhatdl/VNhat)*dmuIhatdl 
     &  + d2muIhatdl2) - Qs/VNhat**2*(muIhat*d2VNhatdl2 
     &  + dVNhatdl*(dmuIhatdl+alphaQ*muIhat-2.*muIhat
     &  / VNhat*dVNhatdl)) + alphaQ*dZINdl
     

        fV = (-1.0 + sqrt(1.0 + 1.0/ZINT))*Qs*muIhat/VNhat
        fV = max(fV,0.001)
!
       else
! Under the conditions of no light:
          ThetaHat      = thetamin  !  a small positive value 
          dThetaHatdl   = 0.
          d2ThetaHatdl2 = 0.

          ZINT          = Qs*zetaN
          dZINdl        = alphaQ*ZINT
          d2ZINdl2      = alphaQ*dZINdl
          fV            = 0.01
          muIhat        = 0.
          dmuIhatdl     = 0.
          d2muIhatdl2   = 0.
       endif

       ! Optimal nutrient quota:
       QN = (1.+ sqrt(1.+1./ZINT))*Qs

       dQNdl  = alphaQ*QN-dZINdl*Qs/(2.*ZINT*sqrt(ZINT*(1.+ZINT))) !confirmed  

        d2QNdl2 = alphaQ*dQNdl - Qs/(2.*ZINT*sqrt(ZINT*(ZINT+1.)))
     &*(d2ZINdl2+alphaQ*dZINdl-(2.*ZINT+1.5)/(ZINT*(ZINT+1.))*dZINdl**2)      ! Confirmed

        dfVdl = alphaQ*Qs*(1/QN+2.*zetaN)-(zetaN+Qs/QN**2)*dQNdl  !Confirmed
!
        d2fVdl2 = (alphaQ**2)*Qs*(1/QN + 2*zetaN) 
     &     -  2.*alphaQ*Qs*dQNdl/QN**2    
     &     +  2.*(dQNdl**2)*Qs/QN**3           
     &     -     (zetaN + Qs/QN**2) * d2QNdl2  ! Confirmed

        self%p%fV = fV

        if (par .gt. I_zero) then
        self%p%muIhat= muIhat
        self%p%QN    = QN
        self%p%VNhat = VNhat
!Net growth rate (d-1) of phytoplankton at the average size
        self%p%muNet = muIhat*(1.-fV-Qs/QN) - zetaN*fV*VNhat
!Here the derivative of muNet includes respiratory costs of both N Assim and Chl maintenance       
        self%p%dmuNetdl = muIhat*(Qs/(QN**2)*dQNdl
     & -  alphaQ*Qs/QN-dfVdl) +  (1.-fV-Qs/QN)*dmuIhatdl 
     & -  zetaN*(fV*dVNhatdl+VNhat*dfVdl)

         self%p%d2muNetdl2 = (Qs/(QN*QN))*dQNdl*dmuIhatdl
     &  + muIhat*(Qs/QN**2)*(dQNdl*(alphaQ - 2.*dQNdl/QN) + d2QNdl2)   
     &  - (alphaQ*(Qs/QN*(dmuIhatdl + alphaQ*muIhat) 
     &  - muIhat*Qs/(QN**2)*dQNdl))  
     &  - (muIhat*d2fVdl2+dfVdl*dmuIhatdl)  
     &  + dmuIhatdl*(Qs/(QN**2)*dQNdl-alphaQ*Qs/QN-dfVdl) 
     &  + (1.-fV-Qs/QN)*d2muIhatdl2     
     &  - zetaN*(fV*d2VNhatdl2 + 2.*dfVdl*dVNhatdl + VNhat*d2fVdl2)  !dC

        else
         self%p%muNet      = 0.
         self%p%dmuNetdl   = 0.
         self%p%d2muNetdl2 = 0.
        endif
!  chl:C ratio [ g chl / mol C ] of the whole cell, at the mean cell size 
          Theta = ThetaHat*(1. - fV - Qs/QN)
! Calculate the mean chl:C ratio of phytoplankton (averaged over all sizes) 
! using the second derivative of the cell quota w.r.t. log size. (Why negative?)
          self%p%ThetaAvg = max(Theta + VAR*.5          
     & * (d2ThetaHatdl2*(1.0 - fV - Qs/QN)         
     & - ThetaHat*(d2fVdl2 + Qs*(dQNdl**2*(2./QN)-d2QNdl2)/(QN**2))),
     &    thetamin)
            
! Calculate the mean N:C ratio of phytoplankton (averaged over all sizes) 
! using the second derivative of the cell quota w.r.t. log size. (How to derive?)
       self%p%QNavg = max(QN/(1.+((2/QN)*dQNdl**2 - d2QNdl2)
     &    * VAR/(2.*QN)),QNmin)  
      ! Calculate the community sinking rate of phytoplankton
      ! Phytoplankton sinking rate at the average size
      w_p         = ScaleTrait(PMU,self%p%W0,self%p%alphaW)
      ! Constrain the sinking rate not too large for large cells (Smayda 1970)
      w_P         = min(w_P,self%env%dtdays*Wmax)
      self%p%dwdl   = self%p%alphaW*w_p
      self%p%d2wdl2 = self%p%alphaW*self%p%dwdl
      ! Sinking rate (m) of phytoplankton
      self%p%w_pAvg =min(w_p+0.5*VAR*self%p%d2wdl2,self%env%dtdays*Wmax)
#ifdef DEBUG_BIO2
          if (self%p%dmuNetdl .ne. self%p%dmuNetdl) then
             write(*,*) 'dmuNetdl is NA!'
             write(*,*) 'PMU    = ',PMU
             write(*,*) 'muIhat = ',muIhat
             write(*,*) 'NO3    = ',NO3
             write(*,*) 'A0hat  = ',A0hat
             write(*,*) 'fA     = ',fA
             write(*,*) 'VNhat  = ',VNhat
             write(*,*) 'Qs     = ',Qs
             write(*,*) 'ZINT   = ',ZINT
             write(*,*) 'dmuIhatdl =',dmuIhatdl
             write(*,*) 'QN     =',QN
             write(*,*) 'dQNdl  = ',dQNdl
             write(*,*) 'fV     = ',fV
             write(*,*) 'dfVdl  = ',dfVdl
             write(*,*) 'VNhat  =' ,VNhat
             write(*,*) 'dVNhatdl =',dVNhatdl
             write (*,*) 'NO3 = ',self%n%Nconc
             write (*,*) 'PAR = ',self%env%light
             write (*,*) 'Temp= ',self%env%temp

             stop
          endif
          if (self%p%d2muNetdl2 .ne. self%p%d2muNetdl2) then
             write(*,*) 'd2muNetdl2 is NA!'
             write(*,*) 'muIhat = ',muIhat
             write(*,*) 'dmuIhatdl =',dmuIhatdl
             write(*,*) 'd2muIhatdl2 =',d2muIhatdl2
             write(*,*) 'QN =',QN
             write(*,*) 'dQNdl = ',dQNdl
             write(*,*) 'fV = ',fV
             write(*,*) 'dfVdl = ',dfVdl
             write(*,*) 'VNhat =' ,VNhat
             write(*,*) 'dVNhatdl =',dVNhatdl
             write(*,*) 'VAR =',VAR
             write (*,*) 'NO3 = ',self%n%Nconc
             write (*,*) 'PAR = ',self%env%light
             write (*,*) 'Temp= ',self%env%temp

             stop
          endif
#endif

      end subroutine phyto
!==============================================================
!The iron subroutine must be called after the ODE subroutine
      subroutine IRONCYCLE(self)
       implicit none
       class(NPZDFE),intent(inout),target :: self
       real :: Tem,Fe_N,PON,DFe,lFe
       real*16 :: keq
!Iron:Nitrogen ratio, Aumont et al. (2003) set Fe/C = 4D-6 (mol:mol). 
!I assume redfield ratio of C/N. Times 1000 to convert umol N to nmol Fe

      PON  = self%PON%Nconc
      Tem  = max(self%Fe%conc,self%Fe%Femin)
      Fe_N = self%Fe%Fe_N

!dFedt = -phytouptake + Zoo excretion - scavenging + remineralization + dust deposition
! (From TOM10 and PISCES)
!     Fe remineralization
      DFe  = Tem+(self%pp(iNO3,iPHY)+self%pp(iNO3,iDET))*Fe_N

!     Phytoplankton uptake
      DFe  = DFe/(1. + self%pp(iPHY,iNO3)*Fe_N/Tem)   

      !lFe: total ligand conc. (Nikelsen et al. Geosci. Model. Dev. 2015)
      lFe  = self%Fe%lFe
!     When iron concentration is above 0.6 nM, it is scavenged by DET
      !Following TOM10Appendix, Eq. 46
      !The equilibrium constant between free iron and ligands and organic complexes. 
      Tem  = self%env%temp    !Discard original Fe conc (Tem)
      keq  = 10**(17.27-1565.7/(273.15 + Tem))
      !lFe: 0.6 nM
      !Following TOM10Appendix, Eq. 45
      !Iron scavenging rate = (Basal scavenging rate + particle asorbtion)*FEprime
      Tem  = 1D0+(lFe-DFe)*keq !Discard temperature
      Tem  = self%Fe%KScm + self%Fe%Ksc*PON*self%z%tf
     &  * (-Tem + sqrt(Tem*Tem + 4.*DFe*keq))/2D0/keq  
      
      self%Fe%conc = DFe/(1. + Tem)

      !Write out Fe fluxes:
      self%Fe%Fe_scav = Tem
      end subroutine
!----------------------------------------------------------------
!-----------------------------------------------------------------------
      pure real function PenAff(logsize,alpha,Pfac,lmin) 
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
! of micro- and mesozooplankton given the total amount of prey available 
! Mesozooplankton preferentially feed on large phytoplankton and microzooplankton
! Change this subroutine to include prey size dependence of food palatability 
! Based on Smith and Vallina (in prep.)
! Also ignore the optimal feeding model by Pahlow (2010)
       implicit none
       class (NPZDFE), intent(inout), target :: self
       integer, parameter :: M=100    !discretize the continous normal distribution
       ! Declare local variables
       real :: tf,PHY,VAR,kp,PMU,Cf,GGE,unass 
     &        ,gmax,beta,Emax,RMzoo,ga,gb,gc,aG,MIC,ZOO,x(M)
     &        ,Ptot,Ptot1,Ptot2,dx,INGES,RES,EGES,gbar,rl,Pl
       integer :: i,j,gform
       logical :: sizepref

       tf   = self%z%tf
       PMU  = self%p%LNV
       VAR  = self%p%VAR
       PHY  = max(self%p%Nconc ,eps)
       ZOO  = max(self%z%Nconc ,eps)
       ga   = self%z%ga
       gb   = self%z%gb
       gc   = self%z%gc
       sizepref = .true.
       if (ga .eq.0. .and. gb.eq.0. .and. gc.eq.0) sizepref = .false.

       ! Calculate the integration  
       if (sizepref .eq. .true.) then
           x(1)  = PMU-3D0*sqrt(VAR)  ! Minimal size
           x(M)  = PMU+3D0*sqrt(VAR)  ! Maximal size
           dx    = (x(M) - x(1))/float(M-1)
           do i  = 2,M-1
             x(i) = x(i-1) + dx 
           enddo
       endif

       kp    =self%z%kp
       gmax  =self%z%gmax
       gform =self%Z%grazing_formula
       aG    =self%Z%alphaG
       ZOO   =max(self%z%Nconc,eps)
#ifdef OCF
       Cf    =self%Z%Cf
       beta  =self%Z%beta
       Emax  =self%Z%Emax
       RMzoo =self%Z%RMzoo
#else
        !For control model:
       unass =self%Z%unass
       GGE   =self%Z%GGE
#endif
        ! Calculate total phytoplankton palatable prey (Ptot,Eq. 13):
  
       if (sizepref .eq. .true.) then
        ! Ptot roughly equals to Phy
         Ptot  = 0.
         Ptot1 = 0.
         do i=1,M-1
           !Probability function
           Pl=PHY*normal(PMU,VAR,x(i))
           !Feeding preference
           rl=rhol(ga,gb,gc,x(i))
           Ptot = Ptot +Pl*dx*rl
          ! Considering killing-the-winner strategy
           Ptot1= Ptot1+Pl**aG*dx*rl
         enddo
       else
         Ptot = PHY
       endif
        ! The grazing dependence on total prey (originally ghat)
         gbar = grazing(gform,kp,Ptot)

#ifdef OCF
        
            call Opt_feed(gmax,Cf,Emax,Ca,beta,
     &           gbar,RMzoo,INGES,RES,EGES)

#else
          !Zooplankton mass-specific ingestion rate (= m*Phy/Zoo in dilution exp., unit: d–1)

          !Zooplankton total ingestion rate
          INGES = tf*gmax*gbar
          !Zooplankton excretion rate (-> DOM)
            RES = INGES*(1.-GGE-unass)
          !ZOOPLANKTON EGESTION (-> POM)
           EGES = INGES*unass
#endif /* OCF  */
          self%z%INGES = max(INGES,0.) 
          self%z%EGES  = max(EGES ,0.)
          self%z%RES   = max(RES  ,0.)

! Grazing rate on the mean size of PHY (specific to N-based Phy biomass, unit: d-1) (Eq. 12)
          if (sizepref .eq. .true.) then
             rl   = rhol(ga,gb,gc,PMU)
             Pl   = normal(PMU,VAR,PMU)*PHY
             gbar = INGES*ZOO*rl*Pl**(aG-1d0)/Ptot1 
             gbar = max(gbar,0.)
             Cf   = 2.*ga*PMU+gb
! First derivative of grazing rate evaluated at the mean size (negative)
             self%z%dgdlbar  = Cf*gbar 
             self%z%d2gdl2bar= ((1.-aG)/VAR+Cf*Cf+2.*ga)*gbar

          else
             gbar = INGES*ZOO/Ptot*sqrt(aG)
             self%z%dgdlbar = 0.
             self%z%d2gdl2bar = (1.-aG)/VAR*gbar
          endif
          self%z%gbar = gbar
       enddo
        
      end subroutine Zooplankton
!--------------------------------------------------------------------------------------
      pure real function grazing(Hollingtype, Ksat, Prey)
        implicit none
        real,   intent(in) :: Ksat, Prey
        integer,intent(in) :: Hollingtype
        ! kp relates to the capture coefficient
        SELECT CASE(Hollingtype)
          ! Holling Type I
          case (1)
            grazing = min(Prey/2.0/Ksat,1.0)
          ! Holling Type II
          case (2)
            grazing = Prey/(Ksat + Prey)  
          ! Holling Type III
          case (3) 
            grazing = Prey*Prey/(Ksat*Ksat + Prey*Prey)
         ! Ivlev
          case (4)
            grazing = 1. - exp(-LOG(2D0)*Prey/Ksat) ! To be consistent with other functions  
      
        END SELECT
        return
      end function
!--------------------------------------------------------------------------------------
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
      pure real function rhol(a,b,c,l)
        implicit none
        real, intent(in) :: a,b,c,l
        rhol = exp(a*l*l+b*l+c)  !c is positive
        return
      end function rhol
!======================================================================
      subroutine ODE(self)
       implicit none
       class(NPZDFE),intent(inout),target :: self
       integer, parameter :: NVAR=4
       integer :: i,j
       real :: PP_nz,PP_nd,PP_dp,PP_dz,PP_pn,PP_zp,PP_np,
     &         NO3,PHY,ZOO,DET,VAR,PMU,Zmort,tf_z,tf_p,
     &         INGES

      tf_z  = self%z%tf
      tf_p  = self%p%tf
      NO3   = max(self%n%Nconc,eps)
      DET   = max(self%d%Nconc,eps)
      PHY   = max(self%p%Nconc,eps)
      ZOO   = max(self%z%Nconc,eps)
      Zmort = ZOO*ZOO*self%z%mz*tf_z  !Mortality term for ZOO
      INGES = self%z%INGES
      ! Save original Size and VAR
      PMU   = self%p%LNV
      VAR   = self%p%VAR
      
      call self%phyto
      call self%Zooplankton

      if (self%z%d2gdl2bar .lt. 0.) then
         self%p%VAR = self%p%VAR*(1.-VAR*self%z%d2gdl2bar)
      else                                               
         self%p%VAR = self%p%VAR/(1.+VAR*self%z%d2gdl2bar)
      endif

      ! Update PMU and VAR:    
      if (self%p%d2muNetdl2 .gt. 0.) then
!      self%dVARdt=VAR*VAR*(self%d2muNetdl2- self%d2gdl2bar-self%d2wdl2)  !Eq. 19 
         self%p%VAR =self%p%VAR*(1.+VAR*self%p%d2muNetdl2)
      else    
         self%p%VAR =self%p%VAR/(1.-VAR*self%p%d2muNetdl2)
      endif
      self%p%VAR    =self%p%VAR/(1.+VAR*self%p%d2wdl2)

      !Eq. 18, Update PMU:
!      self%dPMUdt = self%p%VAR*(self%dmuNetdl - self%dgdlbar - self%dwdl)
      VAR=0.5*(self%p%VAR+VAR)
      if (self%p%dmuNetdl .gt. 0.) then
         self%p%LNV = PMU  + VAR * self%p%dmuNetdl
      else
         self%p%LNV = PMU/(1.-VAR/PMU*self%p%dmuNetdl)
      endif

      if (self%z%dgdlbar .lt. 0.) then
         self%p%LNV = self%p%LNV - VAR*self%z%dgdlbar
      else                                                 
         self%p%LNV = self%p%LNV/(1.+VAR/PMU*self%z%dgdlbar)
      endif

      self%p%LNV    = self%p%LNV/(1.+ VAR/PMU*self%p%dwdl)
#ifdef DEBUG_BIO2
      if ( self%p%LNV .ne. self%p%LNV) then
         write (*,*) 'PMU is NA!'
         write (*,*) 'dmuNetdl = ',self%p%dmuNetdl
         write (*,*) 'dgdlbar  = ',self%z%dgdlbar
         write (*,*) 'dwdl     = ',self%p%dwdl
         write (*,*) 'NO3 = ',self%n%Nconc
         write (*,*) 'PAR = ',self%env%light
         write (*,*) 'Temp= ',self%env%temp
         stop
      endif
      if ( self%p%VAR .ne. self%p%VAR) then
         write (*,*) 'VAR is NA!'
         write (*,*) 'd2muNetdl2 = ',self%p%d2muNetdl2
         write (*,*) 'd2gdl2bar  = ',self%z%d2gdl2bar
         write (*,*) 'd2wdl2     = ',self%p%d2wdl2
         write (*,*) 'NO3 = ',self%n%Nconc
         write (*,*) 'PAR = ',self%env%light
         write (*,*) 'Temp= ',self%env%temp
         stop
      endif
#endif
      !Contrain the PMU and VAR: 
      self%p%LNV = min(max(self%p%LNV,eps),self%p%LNVmax)
      self%p%VAR = min(max(self%p%VAR,eps),self%p%VARmax)
      
! For production/destruction matrix:
      PMU = self%p%LNV
      VAR = self%p%VAR

      self%pp(:,:) = 0.
      self%pp(iNO3,iDET)=max(self%D%rDN*DET*tf_z      ,0.)
      self%pp(iNO3,iZOO)=max(ZOO*self%Z%RES           ,0.)
      self%pp(iDET,iZOO)=max(ZOO*self%Z%EGES+Zmort    ,0.)
      self%pp(iZOO,iPHY)=max(ZOO*self%Z%INGES         ,0.)

      PP_pn = PHY*(self%p%muNet+0.5*self%p%VAR*self%p%d2muNetdl2)
      self%pp(iPHY,iNO3)=max(0.5*( PP_pn + abs(PP_pn)), 0.)
      self%pp(iNO3,iPHY)=max(0.5*(-PP_pn + abs(PP_pn)), 0.)

      self%D%Nconc = (DET + self%pp(iDET,iZOO))/(1.+self%D%rDN*tf_z)

      self%n%Nconc = (NO3 + self%pp(iNO3,iDET)+self%pp(iNO3,iZOO)
     &   + self%pp(iNO3,iPHY))/(1. + self%pp(iPHY,iNO3)/NO3)

      self%p%Nconc = (PHY+self%pp(iPHY,iNO3))
     &   / (1.+(self%pp(iZOO,iPHY) + self%pp(iNO3,iPHY))/PHY)

      self%z%Nconc = (ZOO+self%pp(iZOO,iPHY))/(1.
     &   +  self%Z%EGES + ZOO*self%z%mz*tf_z + self%Z%RES)


      end subroutine ODE
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
      subroutine modified_patankar_2(self)

        IMPLICIT NONE
!
! !INPUT PARAMETERS:
       class(NPZDFE),intent(inout),target :: self
!
! !INPUT/OUTPUT PARAMETER:
       REAL :: cc(1:NVAR)
!     
!! !LOCAL VARIABLES:
       REAL :: pp2(1:NVAR,1:NVAR)
       REAL :: pp1(1:NVAR,1:NVAR)
!       REAL :: a(1:NVAR,1:NVAR),r(1:NVAR)
       REAL :: cc1(1:NVAR),cc2(1:NVAR)
       real :: dFe(2),Fe(0:2),PHY,PMU(0:2),VAR(0:2),Theta(2),QNavg(2),
     &         fV(2),dmudl(2),dgdl(2),dwdl(2),d2mudl2(2),d2gdl2(2),
     &         d2wdl2(2),Wpavg(2)
       integer  :: i,j
!!EOP
!!-----------------------------------------------------------------------
!!BOC
      cc(iNO3)   = max(self%n%Nconc,eps)
      cc(iPHY)   = max(self%p%Nconc,eps)
      cc(iZOO)   = max(self%Z%Nconc,eps)
      cc(iDET)   = max(self%D%Nconc,eps)
      !To avoid unrealistic small and large values due to advection and diffusion:
      self%p%LNV = max(min(self%p%LNV,self%p%LNVmax),eps)
      self%p%VAR = max(min(self%p%VAR,self%p%VARmax),eps)
      PMU(0)     = self%p%LNV
      VAR(0)     = self%p%VAR
      Fe(0)      = self%Fe%conc
!  absolutely essential since not all elements are calculated
      pp2(:,:)=0. ; pp1(:,:)=0.
     
      call self%ODE

#ifdef IRONCYCLE
      call self%IRONCYCLE 
#endif

      pp1        = self%pp
      PMU(1)     = self%p%LNV
      VAR(1)     = self%p%VAR
      Theta(1)   = self%p%ThetaAvg
      QNavg(1)   = self%p%QNavg
      fV(1)      = self%p%fV
      dmudl(1)   = self%p%dmuNetdl
      dwdl(1)    = self%p%dwdl
      d2mudl2(1) = self%p%d2muNetdl2
      d2wdl2(1)  = self%p%d2wdl2
      dgdl(1)    = self%z%dgdlbar
      d2gdl2(1)  = self%z%d2gdl2bar
      Wpavg(1)   = self%p%w_pAvg
      
      call self%ODE

#ifdef IRONCYCLE
      call self%IRONCYCLE
#endif

      Fe(2)      = self%Fe%conc
      pp2        = 0.5*(self%pp+pp1)
      PMU(2)     = self%p%LNV
      VAR(2)     = self%p%VAR
      Theta(2)   = self%p%ThetaAvg
      QNavg(2)   = self%p%QNavg
      fV(2)      = self%p%fV
      dmudl(2)   = self%p%dmuNetdl
      dwdl(2)    = self%p%dwdl
      d2wdl2(2)  = self%p%d2wdl2
      d2mudl2(2) = self%p%d2muNetdl2
      dgdl(2)    = self%z%dgdlbar
      d2gdl2(2)  = self%z%d2gdl2bar
      Wpavg(2)   = self%p%w_pAvg

      self%Fe%dFe         = 0.5*(Fe(2) - Fe(0))
      self%Fe%conc        = Fe(0) + self%Fe%dFe
      self%p%dlnVdt       = 0.5*(PMU(2) - PMU(0))
      self%p%dVARdt       = 0.5*(VAR(2) - VAR(0))
      self%p%LNV          = PMU(0)+self%p%dlnVdt
      self%p%VAR          = VAR(0)+self%p%dVARdt
      self%pp             = pp2      !Averaged fluxes
      self%p%ThetaAvg     = 0.5*(Theta(1)+Theta(2))
      self%p%QNavg        = 0.5*(QNavg(1)+Theta(2))
      self%p%fV           = 0.5*(fV(1)   +   fV(2))
      self%p%dmuNetdl     = 0.5*(dmudl(1)+dmudl(2))
      self%p%d2muNetdl2   = 0.5*(d2mudl2(1)+d2mudl2(2))
      self%z%dgdlbar      = 0.5*(dgdl(1)   +   dgdl(2))
      self%z%d2gdl2bar    = 0.5*(d2gdl2(1) + d2gdl2(2))
      self%p%dwdl         = 0.5*(dwdl(1)+dwdl(2))
      self%p%d2wdl2       = 0.5*(d2wdl2(1)+d2wdl2(2))
      self%p%w_pAvg       = 0.5*(wpavg(1)+wpavg(2))
      self%p%dPHY         = 0.5*(self%p%Nconc - cc(iPHY))
      self%n%dNO3         = 0.5*(self%n%Nconc - cc(iNO3))
      self%z%dZOO         = 0.5*(self%z%Nconc - cc(iZOO))
      self%D%dDET         = 0.5*(self%D%Nconc - cc(iDET))
      self%p%Nconc        = cc(iPHY)+self%p%dPHY
      self%n%Nconc        = cc(iNO3)+self%n%dNO3
      self%z%Nconc        = cc(iZOO)+self%z%dZOO
      self%D%Nconc        = cc(iDET)+self%D%dDET
      end subroutine modified_patankar_2
!--------------------------------------------------------------------------------------
      END MODULE PZEFTFe
#endif
