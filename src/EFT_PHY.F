      subroutine EFT_PHY(NO3, tC, par_,DFe,mu0,KN,A0NMAX,KFe,
     &                   muNet,QN,Theta,LSI,Lno3)
      ! This subroutine use PAHLOW model to calculate phytoplankton growth rate
      ! N:C ratio, and Chl:C ratio based on environmental conditions
      use bio_MOD, only: TEMPBOL,dtdays,tf_p,zetaChl,zetaN,Ep,DO_IRON
     &                  ,nutrient_uptake,RMCHL0
      use lambert, only: WAPR
      implicit none
      !Environmental conditions:
      real, intent(in)    :: NO3, tC, par_,DFe

      !Model parameters:
      real, intent(in)    :: mu0, KN, A0NMAX,KFe
      
      !Output parameters:
      real, intent(out)   :: muNet, QN, Theta, LSI, Lno3
      
      !LOCAL VARIABLES of phytoplankton:
      real :: I_ZERO
      real :: larg   !Environmental variables
      real :: V0hat, A0hat,fA,muIhat
      real :: VNhat
      real :: mu0hat,mu0hatSI,RMchl
      real :: fV, SI
      real :: ZINT
      real :: aI
      real :: ThetaHat
      real :: larg1, W_Y
      real,  parameter :: Qs=0.038   !0.076/2, half of minimal N:C ratio
      real,  parameter :: aI0=0.2D0
      !-----------------------------------------------------------------------
      !Warning: All rates should be multiplied by dtdays to get the real rate
      tf_p  = TEMPBOL(Ep,tC)
      ! Fe related:
      
      mu0hat= dtdays*tf_p*mu0
      
      if (do_IRON) then
        mu0hat= mu0hat *DFe/(DFe + KFe)
      endif
        
        ! Iron limits nitrogen uptake
      V0hat  = mu0hat  
        
      ! Initial slope of P-I curve
      aI     = dtdays*aI0
      ! Cost of photosynthesis
      RMchl  = tf_p*RMCHL0*dtdays
      ! Threshold irradiance and RMchl is set temperature dependent
      I_ZERO = zetaChl*RMchl/aI  
        
      !Define VNhat: the saturation function of ambient nutrient concentration
      SELECTCASE(nutrient_uptake)  
      ! case 1: Classic Michaelis Menton 
        case(1)
      ! Potential maximal nutrient-limited uptake
          ! Half-saturation constant of nitrate uptake
        VNhat  = V0hat*NO3/(NO3 + KN)
      
      ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
        case(2)
      
        A0hat = dtdays*tf_p*A0NMAX
      
        !Define fA
        fA    = 1D0/( 1D0 + sqrt(A0hat*NO3/V0hat) ) 

        VNhat = (1D0-fA)*V0hat*fA*A0hat*NO3
     &        / ((1D0-fA)*V0hat+fA*A0hat*NO3) 
          
        case default
         write(6,*) 'Error: Incorrect option for nutrient uptake!'
         stop
      ENDSELECT  
      
      ! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)
      ! Only calculate within the euphotic zone, otherwise many numerical problems.
      
      if(par_ .gt. I_ZERO) then    !To avoid large C:Chl ratios
        
        larg1   = exp(1d0 + min(aI*par_/(mu0hat*zetaChl),6d2))
        larg    = (1d0 + RMchl/mu0hat)*larg1   
        W_Y     = WAPR(larg,0,0)
       ThetaHat = 1d0/zetaChl + (1d0- W_Y)*mu0hat/(aI * par_)
      
        ! Effect of light limitation
        SI = 1d0 - max(exp(-aI*par_*ThetaHat/mu0hat),0d0)
      
      ! Light dependent growth rate 
      ! (needs to take into account the cost of dark and light-dependent chl maintenance)
        mu0hatSI = mu0hat*SI  ! Gross specific carbon uptake (photosynthesis)
        muIhat   = mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat ! Net specific carbon uptake
        muIhat   = max(muIhat,1D-10*mu0hat)
        LSI      = 1D0-SI
        ZINT     = Qs*(muIhat/VNhat + zetaN)
        fV       = (-1d0 + sqrt(1d0 + 1d0/ZINT))*Qs*muIhat/VNhat
        fV       = max(fV,1D-3)
          
      else
        ! Under the conditions of no light:
      
        ThetaHat      = 0.1  !  Minimal Theta
        ZINT          = Qs*zetaN
        fV            = 1d-3
        muIhat        = 0d0
        LSI           = 1D0
      endif
      
        ! Nutrient limitation index:
        Lno3 =1d0/(1d0 + sqrt(1D0 +1D0/ZINT)) 
        ! Optimal nutrient quota:
        QN   = Qs/Lno3
        
        if (par_ .gt. I_ZERO) then
      !Net growth rate (d-1) of phytoplankton at the average size
          muNet = muIhat*(1d0-2d0*Lno3)
        else
          muNet = 0d0
        endif
      ! Chl:C ratio [ g chl / mol C ] of the whole cell
      Theta  = ThetaHat*(1D0-fV-Qs/QN)
      RETURN   
      END subroutine EFT_PHY
