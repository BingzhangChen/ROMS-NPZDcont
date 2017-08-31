      !-----------------------------------------------------------------------
      !BOP
      !
      ! !MODULE: Add the variable of Chl by Bingzhang Chen on 20150219. 
      ! !Modified NPZD biogeochemical model based upon
      ! Fennel et al. GBC 2006.
      ! Original structure taken from GOTM and adapted for FABM by Jorn Bruggeman. 
      ! !Add temperature components by Bingzhang Chen on 20150213
      !
      ! !INTERFACE:
      module npzdchl
      
      ! !USES:
      use lambert
      
      implicit none
      
      !  default: all is private.
      private
      !
      ! !PUBLIC DERIVED TYPES:
         type,extends(type_base_model),public :: type_gotm_npzdchl
      !     Variable identifiers
            type (type_state_variable_id)        :: id_n,id_p,id_z,id_d,id_chl
            type (type_state_variable_id)        :: id_dic
            type (type_dependency_id)            :: id_par,id_tC  
            type (type_horizontal_dependency_id) :: id_I_0
            type (type_diagnostic_variable_id)   :: id_phygrowth,id_theta,id_grazing,id_dPAR,id_fpar,id_rhochl,id_PP
            type (type_horizontal_diagnostic_variable_id) :: id_dI_0
      
      !     Model parameters
            real(rk) :: w_p,w_d,kc,rmax,gmax,kn,kp,alpha,beta,rdnu,rdnl,rpdu,rpdl,Ea,GGE,unass,mz,thetm,Tr
            real(rk) :: dic_per_n,phycn
         contains
            procedure :: initialize
            procedure :: do
            procedure :: do_ppdd
            procedure :: get_light_extinction
         end type
      !EOP
      !-----------------------------------------------------------------------
      
         contains
      
      !-----------------------------------------------------------------------
      !BOP
      !
      ! !IROUTINE: Initialise the NPZDCHL model
      !
      ! !INTERFACE:
         subroutine initialize(self,configunit)
      !
      ! !DESCRIPTION:
      !  Here, the npzdchl namelist is read and the variables exported
      !  by the model are registered with FABM.
      !
      ! !INPUT PARAMETERS:
         class (type_gotm_npzdchl), intent(inout), target :: self
         integer,                   intent(in)            :: configunit
      !
      ! !LOCAL VARIABLES:
         real(rk)          :: n_initial
         real(rk)          :: p_initial
         real(rk)          :: z_initial
         real(rk)          :: d_initial
         real(rk)          :: chl_initial
         real(rk)          :: w_p
         real(rk)          :: w_d
         real(rk)          :: kc
         real(rk)          :: rmax
         real(rk)          :: gmax
         real(rk)          :: kn
         real(rk)          :: kp
         real(rk)          :: alpha
         real(rk)          :: beta
         real(rk)          :: rdnu
         real(rk)          :: rdnl
         real(rk)          :: rpdu
         real(rk)          :: rpdl
         real(rk)          :: dic_per_n
         real(rk)          :: Ea
         real(rk)          :: GGE
         real(rk)          :: unass
         real(rk)          :: mz
         real(rk)          :: thetm
         real(rk)          :: phycn
         real(rk)          :: Tr
      
         character(len=64) :: dic_variable
      
         real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
         namelist /gotm_npzdchl/ n_initial,p_initial,z_initial,d_initial,phycn,chl_initial,&
                                  w_p,w_d,kc,rmax,gmax,kn,kp,alpha,beta,  &
                                  rdnu,rdnl,rpdu,rpdl,dic_variable,dic_per_n,Ea,&
                                  GGE,unass,mz,thetm,Tr
      !EOP
      !-----------------------------------------------------------------------
      !BOC
         n_initial = 4.5_rk
         p_initial = 0.1_rk
         z_initial = 0.1_rk
         d_initial = 4.5_rk      
         phycn     = 6.625_rk         !C:N molar ration in phytoplankton 
         chl_initial = 0.15_rk        !Initial Chl concentration (mg m-3)
         w_p       = -0.1_rk          !Sinking rate of phytoplankton
         w_d       = -1.0_rk          !Sinking rate of detritus
         kc        = 0.025_rk         !Light attenuation of chlorophyll (mg Chl m-2)
         rmax      = 1.0_rk           !Maximal growth rate of phytoplankton at 15 ºC
         gmax      = 2.0_rk           !Maximal biomass-specific ingestion rate of zooplankton at 15 ºC
         kn        = 0.1_rk           !half-saturation constant of nutrient uptake of phytoplankton 
         kp        = 0.2_rk           !half-saturation constant for zooplankton grazing
         alpha     = 0.125_rk         !Initial slope of the P-I curve (mol C gChl-1 (W m-2)-1 d-1).
         beta      = 0.000125_rk      !photoinhibition parameter in Platt (1980) model, same unit with alpha
         rdnu      = 0.03_rk          !The conversion rate from detritus to nutrient at 15 C below the euphotic zone
         rdnl      = 0.01_rk          !The conversion rate from detritus to nutrient at 15 C in the euphotic zone
         rpdu      = 0.10_rk          !Phytoplankton mortality rate at 15 C in the euphotic zone
         rpdl      = 0.02_rk           !Phytoplankton mortality rate at 15 C below the euphotic zone 
         dic_per_n = 106.0_rk/16.0_rk  ! Redfield C:N
         dic_variable = ''
         Ea        = 0.65_rk           !Activation energy
         GGE       = 0.30_rk      !GGE of zooplankton, independent of temperature, taken from Buitenhuis et al. GBC (2010)
         unass     = 0.20_rk      !The fraction of unassimilated food by zooplankton egested -> detritus (Buitenhuis et al. 2010) 
         mz        = 0.025_rk     !Coefficient for zooplankton mortality at 15 C
         thetm     = 0.6_rk          !maximal chlorophyll (g/cell) -to-phytoplankton (molC/cell) ratio, taken from Fennel et al. GBC (2006)
         Tr        = 4.0_rk
      
         ! Read the namelist
         if (configunit>0) read(configunit,nml=gotm_npzdchl,err=99,end=100)
      
         ! Store parameter values in our own derived type
         ! NB: all rates must be provided in values per day,
         ! and are converted here to values per second.
         call self%get_parameter(self%kc,'kc','m2 mgchl-1','specific light extinction of Chl',default=kc)
         call self%get_parameter(self%rmax,'rmax','d-1','maximum specific growth rate of phytoplankton',default=rmax,scale_factor=d_per_s)
         call self%get_parameter(self%gmax,'gmax','d-1','maximum specific ingestion rate of zooplankton',default=gmax,scale_factor=d_per_s)
      
         call self%get_parameter(self%kp,'kp','(mmol N m-3)^2','grazing half-saturation constant',default=kp)
         call self%get_parameter(self%kn,'kn','mmol N m-3','half-saturation nutrient concentration for phytoplankton',default=kn)
         call self%get_parameter(self%alpha,'alpha','mol C (gChl)-1 d-1 (W m-2)-1','initial slope of the P-I curve',default=alpha,scale_factor=d_per_s)
         call self%get_parameter(self%beta,'beta','mol C (gChl)-1 d-1 (W m-2)-1','photoinhibition parameter',default=beta,scale_factor=d_per_s)
         call self%get_parameter(self%rdnu,'rdnu','d-1','detritus remineralization rate below euphotic zone',default=rdnu,scale_factor=d_per_s)
         call self%get_parameter(self%rdnl,'rdnl','d-1','detritus remineralization rate in euphotic zone',default=rdnl,scale_factor=d_per_s)  
         call self%get_parameter(self%rpdu,'rpdu','d-1','phytoplankton mortality in euphotic zone',default=rpdu,scale_factor=d_per_s)
         call self%get_parameter(self%rpdl,'rpdl','d-1','phytoplankton mortality below euphotic zone',default=rpdl,scale_factor=d_per_s)
         call self%get_parameter(self%mz,'mz','d-1 (mmol N m-3)-1','zooplankton mortality',default=mz,scale_factor=d_per_s)
         call self%get_parameter(self%dic_per_n,'dic_per_n','-','C:N ratio of biomass',default=dic_per_n)
         call self%get_parameter(w_p,'w_p','m d-1','vertical velocity of phytoplankton (<0 for sinking)',default=w_p, scale_factor=d_per_s)
         call self%get_parameter(w_d,'w_d','m d-1','vertical velocity of detritus  (<0 for sinking)',default=w_d,scale_factor=d_per_s)
         call self%get_parameter(self%Ea,'Ea','eV','activation energy',default=Ea)
         call self%get_parameter(self%GGE,'GGE','-','Gross growth efficiency of zooplankton',default=GGE)
         call self%get_parameter(self%unass,'unass','-','fraction of egested',default=unass)
         call self%get_parameter(self%thetm,'thetm','gChl/molC','maximal chlorophyll-to-phytoplankton ratio',default=thetm)
         call self%get_parameter(self%phycn,'phycn','molC/molN','phytoplankton C:N ratio',default=phycn)
         call self%get_parameter(self%Tr,'Tr','ºC','reference temperature',default=Tr)    
      
         ! Register state variables
         call self%register_state_variable(self%id_n,'nut','mmolN/m**3','nutrients',    n_initial,minimum=0.0_rk,no_river_dilution=.true.)
         call self%register_state_variable(self%id_p,'phy','mmolN/m**3','phytoplankton',p_initial,minimum=0.0_rk,vertical_movement=w_p)
         call self%register_state_variable(self%id_z,'zoo','mmolN/m**3','zooplankton',  z_initial,minimum=0.0_rk)
         call self%register_state_variable(self%id_d,'det','mmolN/m**3','detritus',     d_initial,minimum=0.0_rk,vertical_movement=w_d)
         call self%register_state_variable(self%id_chl,'chl','mg/m**3','chlorophyll', chl_initial,minimum=0.0_rk,vertical_movement=w_p,no_river_dilution=.true.)
      
         ! Register the contribution of all state variables to total nitrogen
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_n)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_p)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_z)
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_d)
      
         ! Register link to external DIC pool, if DIC variable name is provided in namelist.
         call self%register_state_dependency(self%id_dic,'dic','mmol/m**3','total dissolved inorganic carbon',required=.false.)
         if (dic_variable/='') call self%request_coupling(self%id_dic,dic_variable)
      
         ! Register diagnostic variables
         call self%register_diagnostic_variable(self%id_dPAR,'PAR','W/m**2','photosynthetically active radiation',&
                           output=output_time_step_averaged)
         call self%register_diagnostic_variable(self%id_phygrowth,'Phygrowth','per d','Phytoplankton growth rate',&
                           output=output_time_step_averaged)
         call self%register_diagnostic_variable(self%id_grazing,'grazing','per d','Zooplankton grazing rate',     &
                           output=output_time_step_averaged)
         call self%register_diagnostic_variable(self%id_theta,'theta','gChl/molC','phytoplankton chl-to-carbon ratio',               output=output_time_step_averaged)
         call self%register_diagnostic_variable(self%id_fpar,'fpar','','light limitation index',     & 
                           output=output_time_step_averaged)
         call self%register_diagnostic_variable(self%id_rhochl,'rhochl','','chl synthesis factor', output=output_time_step_averaged)
         call self%register_HORIZONTAL_DIAGNOSTIC_variable(self%id_dI_0,'I_0','W/m**2','photosynthetically active radiation',&
                           output=output_time_step_averaged)
         call self%register_diagnostic_variable(self%id_PP,'PP','mgC/m**3','primary production', output=output_time_step_averaged)
        
       
      
         ! Register environmental dependencies
         call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
         call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
         call self%register_dependency(self%id_tC,  standard_variables%temperature)
      
         return
      
      99 call self%fatal_error('gotm_npzdchl_init','Error reading namelist gotm_npzdchl.')
      
      100 call self%fatal_error('gotm_npzdchl_init','Namelist gotm_npzdchl was not found.')
      
         end subroutine initialize
      !EOC
      
      !-----------------------------------------------------------------------
      !BOP
      !
      ! !IROUTINE: Right hand sides of NPZD model
      !
      ! !INTERFACE:
         subroutine do(self,_ARGUMENTS_DO_)
      !
      !
      ! !INPUT PARAMETERS:
         class (type_gotm_npzdchl),intent(in) :: self
         _DECLARE_ARGUMENTS_DO_
      !
      ! !REVISION HISTORY:
      !  Original author(s): Hans Burchard, Karsten Bolding
      !
      ! !LOCAL VARIABLES:
         real(rk)                   :: n,p,z,d,chl,par,I_0,tC,theta,rmax_T
         real(rk)                   :: I_min,rpd,rdn,fpar,growth,rhochl,graz_zoo,dn,dp,dz,dd,dchl
         real(rk), parameter        :: secs_pr_day = 86400.0
      !EOP
      !-----------------------------------------------------------------------
      !BOC
         ! Enter spatial loops (if any)
         _LOOP_BEGIN_
      
         ! Retrieve current (local) state variable values.
         _GET_(self%id_n,n) ! nutrient
         _GET_(self%id_p,p) ! phytoplankton
         _GET_(self%id_z,z) ! zooplankton
         _GET_(self%id_d,d) ! detritus
         _GET_(self%id_chl,chl) ! chlorophyll
      
         ! Retrieve current environmental conditions.
         _GET_(self%id_par,par)             ! local photosynthetically active radiation
         _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation
         _GET_(self%id_tC,tC)               ! local temperature in degree celcius !added by Bingzhang 
      
         ! Deal with zero I_0 on the first day
         if (par .le. 0) par = par + 0.001  
      
         ! Define I_min as 1% of I_0 (bottom of euphotic layer)
         I_min = 0.01*I_0
         
         ! Loss rate of phytoplankton to detritus depends on local light intensity.
      
         if (par .ge. I_min) then   ! euphotic layer
            rpd = self%rpdu
            rdn = self%rdnl    ! conversion rate of detritus to nutrient depends on light (light inhibition of nitrification)
         else
            rpd = self%rpdl
            rdn = self%rdnu
         end if
      
         ! Define phytoplankton Chl-to-carbon ratio (gChl/molC)
         theta = chl/(p*self%phycn)
      
         ! The maximal growth rate (rmax_T) under temperature tC 
         rmax_T = self%rmax*temp_f(self,tC)
      
         !The light limitation index (fpar)
         fpar = (1-exp(-self%alpha*par*theta/rmax_T))*exp(-self%beta*par*theta/rmax_T)
      
         ! Growth rate (d-1) is the product function of temperature, light, and nutrient   
         growth = rmax_T*(n/(n+self%kn))*fpar
      
         ! define rhochl (gChl/molC;the fraction of phytoplankton carbon production that is devoted to Chl synthesis)
         rhochl = self%thetm*growth/(self%alpha*par*theta)
      
         ! The total amount of phytoplankton grazed by zooplankton (molN;gmax is the maximal specific ingestion rate!)
         graz_zoo = z*self%gmax*temp_f(self,tC)*p**2/(self%kp+p**2)
         
         ! Phytoplankton are assumed to release no nitrogen. Some of them are only converted to detritus (rpd).
         ! rpd depends on light because the nitrification rate is inversely related with light level 
         ! dn = -phytoplankton uptake + zooplankton excretion (80% of (1-GGE)*total_grazed) 
         !      + detritus conversion (detritus to NH4 and nitrification)
         dn = -p*growth + graz_zoo*(1-self%GGE)*(1-self%unass) + rdn*d*temp_f(self,tC)
      
         ! dp = phytoplankton production - zooplankton grazing - phytoplankton->detritus (molN d-1)
         dp = p*growth  - graz_zoo - rpd*p*temp_f(self,tC)
      
         ! dz = zooplankton growth(excluding excretion and egestion) - zooplankton mortality (molN d-1)
         dz = graz_zoo*self%GGE - z**2*self%mz*temp_f(self,tC)
      
         ! dD = phytoplankton mortality + zooplankton mortality + zooplankton egestion - conversion-to-nutrient (molN d-1)
         dd = rpd*p*temp_f(self,tC) + z**2*self%mz*temp_f(self,tC) + graz_zoo*(1-self%GGE)*self%unass - rdn*d*temp_f(self,tC) 
         ! dChl = chl synthesis - zooplankton grazing on Chl - Chl->detritus (Fennel et al. GBC 2006)
         dchl = rhochl*growth*p*self%phycn - graz_zoo*chl/p - rpd*chl*temp_f(self,tC)
       
         ! Set temporal derivatives
         _SET_ODE_(self%id_n,dn)
         _SET_ODE_(self%id_p,dp)
         _SET_ODE_(self%id_z,dz)
         _SET_ODE_(self%id_d,dd)
         _SET_ODE_(self%id_chl,dchl)
      
         ! If an externally maintained DIC pool is present, change the DIC pool according to the
         ! the change in nutrients (assuming constant C:N ratio)
         if (_AVAILABLE_(self%id_dic)) _SET_ODE_(self%id_dic,self%dic_per_n*dn)
      
         ! Export diagnostic variables
         _SET_DIAGNOSTIC_(self%id_dPAR,          par)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dI_0,I_0)
         _SET_DIAGNOSTIC_(self%id_theta,         theta)
         _SET_DIAGNOSTIC_(self%id_rhochl,        rhochl)
         _SET_DIAGNOSTIC_(self%id_phygrowth,     growth*secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_grazing,       (graz_zoo/p)*secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_PP,            growth*secs_pr_day*p*12*self%phycn)
      
         ! Leave spatial loops (if any)
         _LOOP_END_
      
         end subroutine do
      !-----------------------------------------------------------------------
      
      
      
      !-----------------------------------------------------------------------
      !BOP
      !
      ! !IROUTINE: Get the light extinction coefficient due to biogeochemical
      ! variables
      ! !Revised by Bingzhang Chen on 20150220 to only include the extinction of chl
      ! !INTERFACE:
         subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
      !
      ! !INPUT PARAMETERS:
         class (type_gotm_npzdchl), intent(in) :: self
         _DECLARE_ARGUMENTS_GET_EXTINCTION_
      !
      ! !REVISION HISTORY:
      !  Original author(s): Jorn Bruggeman
      !
      ! !LOCAL VARIABLES:
         real(rk)                     :: chl
         ! Enter spatial loops (if any)
         _LOOP_BEGIN_
      
         ! Retrieve current (local) state variable values.
         _GET_(self%id_chl,chl) ! chlorophyll
      
      
      
         ! Self-shading with explicit contribution from background phytoplankton concentration.
         _SET_EXTINCTION_(self%kc*chl)
      
         ! Leave spatial loops (if any)
         _LOOP_END_
      
         end subroutine get_light_extinction
      !EOC
      !-----------------------------------------------------------------------
      
      
      
      !-----------------------------------------------------------------------
      !BOP
      !
      ! !IROUTINE: Right hand sides of NPZDCHL model exporting production/destruction matrices
      !
      ! !INTERFACE:
         subroutine do_ppdd(self,_ARGUMENTS_DO_PPDD_)
      
      !! !INPUT PARAMETERS:
         class (type_gotm_npzdchl),intent(in) :: self
         _DECLARE_ARGUMENTS_DO_PPDD_
      !
      ! !LOCAL VARIABLES:
         real(rk)                   :: n,p,z,d,chl,par,I_0,tC,theta,rmax_T
         real(rk)                   :: I_min,rpd,rdn,fpar,growth,rhochl,graz_zoo,dn
         real(rk), parameter        :: secs_pr_day = 86400.
      !EOP
      !-----------------------------------------------------------------------
      !BOC
         ! Enter spatial loops (if any)
         _LOOP_BEGIN_
      
         ! Retrieve current (local) state variable values.
         _GET_(self%id_n,n) ! nutrient
         _GET_(self%id_p,p) ! phytoplankton
         _GET_(self%id_z,z) ! zooplankton
         _GET_(self%id_d,d) ! detritus
         _GET_(self%id_chl,chl) ! chlorophyll
         ! Retrieve current environmental conditions.
         _GET_(self%id_par,par)  ! local photosynthetically active radiation
         _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation
         _GET_(self%id_tC,tC)     !local temperature
      
         if (par .le. 0) par = par + 0.001 
      
         ! Define I_min as 1% of I_0 (bottom of euphotic layer)
         I_min = 0.01*I_0
      
         ! Loss rate of phytoplankton to detritus depends on local light intensity.
         if (par .ge. I_min) then
            rpd = self%rpdu
            rdn = self%rdnl    ! conversion rate of detritus to nutrient depends on light (light inhibition of nitrification)
         else
            rpd = self%rpdl
            rdn = self%rdnu
         end if
      
         ! Define phytoplankton Chl-to-carbon ratio (gChl/molC)
         theta = chl/(p*self%phycn) 
      
         ! The maximal growth rate (rmax_T) under temperature tC 
         rmax_T = self%rmax*temp_f(self,tC)
      
         !The light limitation index (fpar)
         fpar = (1.-exp(-self%alpha*par*theta/rmax_T))*exp(-self%beta*par*theta/rmax_T)
      
         ! Growth rate (d-1) is the product function of temperature, light, and nutrient   
         growth = rmax_T*(n/(n+self%kn))*fpar
      
         ! define rhochl (gChl/molC;the fraction of phytoplankton carbon production that is devoted to Chl synthesis)
         rhochl = (self%thetm)*growth*p*(self%phycn)/(self%alpha*par*chl)
         ! The total amount of phytoplankton grazed by zooplankton (molN;gmax is the maximal specific ingestion rate!)
         graz_zoo = z*self%gmax*temp_f(self,tC)*p**2./(self%kp+p**2.)
      
         ! Nutrient changes. 
         dn = -p*growth + graz_zoo*(1.-self%GGE)*(1.-self%unass) + rdn*d*temp_f(self,tC)
         ! Assign destruction rates to different elements of the destruction matrix.
         ! By assigning with _SET_DD_SYM_(i,j,val) as opposed to _SET_DD_(i,j,val),
         ! assignments to dd(i,j) are automatically assigned to pp(j,i) as well.
         _SET_DD_SYM_(self%id_n,self%id_p,p*growth)                                ! snp
         _SET_DD_SYM_(self%id_p,self%id_z,graz_zoo)                                ! spz
         _SET_DD_SYM_(self%id_z,self%id_n,graz_zoo*(1.-self%GGE)*(1.-self%unass))  ! szn
         _SET_DD_SYM_(self%id_d,self%id_n,rdn*d*temp_f(self,tC))                   ! sdn
         _SET_DD_SYM_(self%id_p,self%id_d,rpd*p*temp_f(self,tC))                   ! spd
         _SET_DD_SYM_(self%id_z,self%id_d,z**2.*self%mz*temp_f(self,tC)+graz_zoo*(1.-self%GGE)*self%unass)       ! szd
         _SET_PP_(self%id_chl,self%id_chl,rhochl*growth*p*self%phycn)
         _SET_DD_(self%id_chl,self%id_chl,graz_zoo*chl/p)
         _SET_DD_(self%id_chl,self%id_chl,rpd*chl*temp_f(self,tC))
      
      	! If an externally maintained DIC pool is present, change the DIC pool according to the
      	! the change in nutrients (assuming constant C:N ratio)
         if (_AVAILABLE_(self%id_dic)) _SET_PP_(self%id_dic,self%id_dic,self%dic_per_n*dn)
      
      	! Export diagnostic variables
      	_SET_DIAGNOSTIC_(self%id_dPAR,                              par)
              _SET_DIAGNOSTIC_(self%id_fpar,               fpar)
      	_SET_DIAGNOSTIC_(self%id_phygrowth,          growth*secs_pr_day)
              _SET_DIAGNOSTIC_(self%id_grazing,      (graz_zoo/p)*secs_pr_day)
              _SET_DIAGNOSTIC_(self%id_theta,                           theta)
              _SET_DIAGNOSTIC_(self%id_PP,            growth*secs_pr_day*p*12.*self%phycn)
      ! Leave spatial loops (if any)
      	_LOOP_END_
      
         end subroutine do_ppdd
      !EOC
      !--------------------------------------------------------------------------------------------------
       
      !---------------------------------------------------------
       pure real(rk) function temp_f(self,tC)
      ! !IROUTINE: Temperature dependence of plankton rates. Added by Bingzhang Chen on 20150213
      ! !DESCRIPTION:
      ! Here, the temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
      ! tC: in situ temperature
      ! Tr: reference temperature
      !
      ! !INPUT PARAMETERS:
         type (type_gotm_npzdchl), intent(in) :: self
         real(rk), intent(in)         :: tC
         real(rk), parameter          :: kb = 8.62E-5 ! boltzman constant constant [ eV /K ]
         
         temp_f =  exp(-(self%Ea/kb)*(1.0_rk/(273.15_rk+tC) - 1.0_rk/(273.15_rk + self%Tr)))
      
      
         end function temp_f
      !EOC
      
      !-----------------------------------------------------------------------
      
         end module gotm_npzdchl
      
      !-----------------------------------------------------------------------
      ! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
      !-----------------------------------------------------------------------
