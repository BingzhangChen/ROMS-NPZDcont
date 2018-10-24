      Module BIO_MOD
      !The biological module that provides common parameters for all biological models
      implicit none
      integer, parameter :: nutrient_uptake=2, grazing_formulation=3
      real,    parameter :: RMCHL0   =0.1
      real,    parameter :: GMAX     =1.0
      real,    parameter :: KP       =0.5
      real,    parameter :: MZ       =0.1
      real,    parameter :: RDN      =0.1
      real,    parameter :: WDET     =2D1
      logical, parameter :: DO_IRON        =.TRUE.
      logical, parameter :: kill_the_winner=.TRUE.
      real,    parameter :: thetamin = 0.02
      real,    parameter :: thetamax = 0.62
      real,    parameter :: mu0a     = 0.8   !For the size model
      real,    parameter :: mu0b     = 0.3   !For the three trait model
      real,    parameter :: Qmin     = 0.06 
      real,    parameter :: Qmax     = 0.18
      real,    parameter :: alphaI   = -0.13,aI0_C   =0.03
      real,    parameter :: alphaK   = 0.27, K0N     =0.2
      real,    parameter :: Femin    = 0.01   ! Minimal dissolved Fe conc.
      real,    parameter :: alphamu  = 0.2,  alphaG  =1.0d0
      real,    parameter :: betamu   = -0.01

      !Temperature senstivity tuned by the algorithm
      real,    parameter :: Ep     =0.41, Ez=0.65
      real,    parameter :: K0Fe   =0.02  !Unit: nM. Gregg (2003) used ~0.1 nM for phyto
      real,    parameter :: alphaFe=0.27
      real,    parameter :: zetaChl=0.6, zetaN =0.8
      real,    parameter :: VTR    =0.1    ! 0.1 for highest TD
      real               :: dtdays
      real               :: dx
      real               :: tf_p
      CONTAINS
      !========================================================
      real function TEMPBOL(Ea,tC)
      implicit none
      !DESCRIPTION:
      !The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
      ! tC: in situ temperature
      ! Tr: reference temperature
      !
      !INPUT PARAMETERS:
      real, intent (in) :: Ea, tC
      ! boltzman constant constant [ eV /K ]
      real, parameter   :: kb = 8.62d-5, Tr = 15D0
      
      TEMPBOL = exp(-(Ea/kb)*(1D0/(273.15 + tC)-1D0/(273.15 + Tr)))
      return 
      end function TEMPBOL
      !====================================================
      real function ScaleTrait( logsize, star, alpha ) 
      implicit none
      real, intent(IN) :: logsize, star, alpha
      
      ! Calculate the size-scaled value of a trait
      ! for the given log (natural, base e) of cell volume as pi/6*ESD**3 (micrometers). 
      
      ScaleTrait = star * exp( alpha * logsize )
      
      return
      end function ScaleTrait
      !====================================================
      real function PenAff( logsize, alpha, Pfac, lmin ) 
      implicit none
      real, intent(IN) :: logsize, alpha, Pfac, lmin 
      
      !A 'penalty' function to reduce the value of affinity for nutrient at very small cell sizes
      !in order to avoid modeling unrealistically small cell sizes.  This is needed because affnity
      !increases with decreasing cell size, which means that under low-nutrient conditions, without
      !such a penalty, unrealistically small cell sizes could be predicted.
      !This penalty function becomes zero at logsize = lmin.   
         
        PenAff = 1.0 - exp(Pfac*alpha*(logsize - lmin))
      end function PenAff
      !====================================================
      real function grazing(Hollingtype, Ksat, Prey)
      implicit none
      real,    intent(in) :: Ksat, Prey
      integer, intent(in) :: Hollingtype
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
          grazing = min(Prey*Prey/(Ksat*Ksat + Prey*Prey), 1D0)
       ! Ivlev
        case (4)
       !To be consistent with other functions  
          grazing = 1d0-exp(-log(2d0)*Prey/Ksat) 
      
      END SELECT
      return
      end function grazing


      subroutine Calculate_PAR(I_0, nlev_, Hz, Chl, PAR)
        !top level is nlev_, bottom layer is 1, following ROMS convention
        implicit none
        real   , intent(in) :: I_0
        integer, intent(in) :: nlev_    ! Total number of vertical layers
        real   , intent(in) :: Hz(nlev_), Chl(nlev_)
        real   , intent(out):: PAR(nlev_)
        integer             :: i
        real                :: par0, attn    ! Scratch variable
        real   , parameter  :: kw = 0.04, kc = 0.025
      
        par0 = I_0   !Light at the grid surface
        do i = nlev_,1,-1
           attn   = exp(-0.5*Hz(i)*(kw + kc*Chl(i)))
           PAR(i) = par0*attn
           par0   = PAR(i)*attn
        enddo
      end subroutine Calculate_PAR

      subroutine NPZDPhy_size(PMU, NO3, tC, par_, muNet, SI, Lno3)
      implicit none
      
      !INPUT PARAMETERS:
      real, intent(in)    :: PMU, NO3, tC, par_
      
      !Output parameters:
      real, intent(out)   :: muNet, SI, Lno3
      
      !Local variable:
      real                :: Kn, mu0hat
      real, parameter     :: alphamu=0.08
      real, parameter     :: alphaKN=0.26
      real, parameter     ::     mu0=2.5d0
      real, parameter     ::     K0N=0.11d0
      real, parameter     ::     aI0=0.01d0
      
      tf_p    = TEMPBOL(Ep,tC)
      mu0hat  = dtdays*tf_p * mu0 * exp(alphamu*PMU)
      
      Kn      = K0N * exp( alphaKN * PMU)
      
      ! Effect of light limitation
      SI      = 1d0 - max(exp(-aI0 * par_ ),0d0)
      Lno3    = NO3/(NO3+Kn)
      muNet   = mu0hat * NO3 / (NO3 + Kn) * SI
      return
      end subroutine NPZDPhy_size
      !===================================================
      END Module BIO_MOD
