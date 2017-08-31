#include "cppdefs.h"
#if defined BIOLOGY && defined BIO_FLEX

      subroutine biology_tile (Istr,Iend,Jstr,Jend)
      USE bio_MOD
      USE forces, only: srflx
!
! Compute biological forcing functions
!
! In this particular implementation there is 4 compartments:
! NO3, PHYTOplankton, ZOOplanknton, DETritus.
!
      implicit none
      integer Istr,Iend,Jstr,Jend
#include "param.h"
#include "grid.h"
#include "ocean3d.h"
#include "ocean2d.h"
#include "diagnostics.h"
#include "scalars.h"
#include "mixing.h"

      integer, parameter :: nsink   = NumVSinkTerms
      real,    parameter :: kw      = 0.04, 
     &                      kc      = 0.025, !Light attenuation of chlorophyll (mg Chl m-2)
     &                      mon2sec = 2592D3,!How many seconds that one month contains
     &                      eps     = 1D-20
      real,    parameter :: d_per_s = 864d2   ! how many seconds in one day

#ifdef DIAGNOSTICS_BIO
      integer l, iflux
      real ThisVSinkFlux(N, NumVSinkTerms),  ! [mmol m-2 s-1], upward flux is positive
     &     ThisFlux(N, NumFluxTerms)
#endif

      integer i,j,k, iB
      real    NO3(N), PHY(N), ZOO(N), DET(N), TEMP(N), FC(0:N),
     &        PARsup,par, I_0, theta, rmax_T, thetaHat, dtdays,attn,
     &        Q, larg, I_min, Ith, rpd,fpar,rIhat,vNhat,fA,ZINT,fV,
     &        Rchl,mu,dn,
     &        dp,dd,dz,ghat,g_th,aAt,aAf,INGES,AE,RES,NGR,EGES,
     &        SB(N,nsink),dSB(0:N,nsink),wSB(nsink),cff,cff1,cff2,cff6,
     &        TwPHYT, TwDET, Trmax,Tgmax,Talpha,Trdn,Trpdu,Trpdl,Tmz,
     &        TRMchl,TRMzoo,TV0,TA0,aR,aL,cu
# if defined DIAGNOSTICS_BIO
     &      , LastVSinkFlux,ColumnMassOld(NumVSinkTerms)
     &      , ColumnMassNew(NumVSinkTerms)
# endif /* DIAGNOSTICS_BIO */

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
        NO3(k) =max(t(i,j,k,nnew,ino3_)  ,0.)   ! Nitrate
        PHYT(k)=max(t(i,j,k,nnew,iphy1)  ,0.)   ! Phytoplankton
        ZOO(k) =max(t(i,j,k,nnew,izoo1)  ,0.)   ! Zooplankton
        DET(k) =max(t(i,j,k,nnew,idet1)  ,0.)   ! Detritus
        TEMP(k)=max(t(i,j,k,nnew,itemp),  0.)   ! temperature; [deg. C]
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
        I_0 = max(I_0,1E-6)  ! Force I_0 to be positive
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

          attn=exp(-0.5*(kw + kc*PHYT(k))*     ! Attenuation in this grid
     &            (z_w(i,j,k)-z_w(i,j,k-1)))   ! z_w: depth of the grid (negative values!).

          PARsup = par * attn       ! Light level at the middle of the grid

         ! Loss rate of phytoplankton to detritus depends on local light intensity.
          if (PARsup .ge. I_min) then
              rpd = Trpdu
           else
              rpd = Trpdl
          endif
       ! The maximal growth rate (rmax_T) under temperature temp(k) 
            rmax_T = Trmax*tf(Ea,TEMP(k))
              
       ! Threshold irradiance and RMchl is set temperature dependent
            Ith = zetaChl * TRMchl * tf(Ea,TEMP(k))/Talpha  

       ! Calculate optimal thetahat (ThetaHat) under nutrient replete conditions!
            if(fV_opt .eq. .true. ) then 
                 if( PARsup .gt. Ith ) then
                    larg = (1.0 + TRMchl/Trmax) * 
     &                     exp(1.0 + Talpha*PARsup/(rmax_T*zetaChl))      

                    ThetaHat = 1.0/zetaChl +
     &                     (1.0-WAPR(larg,0,0))*tf(Ea,TEMP(k))
     &                     * rmax/(alpha*PARsup)
                 else
                    ThetaHat = 1e-6  !  a small positive value 
                 endif
              else
                 ThetaHat = 0.6
              endif
             !The light limitation index (fpar)
             fpar = 1.-exp(-alpha*PARsup*ThetaHat/rmax_T)
             
             ! Light dependent growth rate (minus the cost of chl maintenance)
             ! Also assume that chl maintenance is temperature dependent
             rIhat = rmax_T*fpar - (rmax_T*fpar+RMchl*tf(Ea,temp(k)))
     &               *zetaChl*ThetaHat
             ! Define vNhat: the saturation function of ambient nutrient concentration
             if (nutrient_uptake == 1) then 
             ! case 1: Classic Michaelis Menton 
                 fA = 0.5
             ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
               else
                 !Define fA
                 fA = 1.0/(1.0+sqrt(TA0*NO3(k)/(tf(Ea,temp(k))*TV0))) 
             endif
         
             ! Calculate vNhat based on optimal allocation between nutrient uptake and assimilation
             vNhat=(1.0-fA)*TV0*fA*TA0*NO3(k)/((1.0-fA)
     &             *TV0+fA*TA0*NO3(k)) 

                ! Define ZINT (see Smith et al. JPR 2015)
             ZINT = 0.5*Q0*(rIhat/vNhat + zetaN)
                
                ! Estimate optimal Q
              if(fV_opt .eq. .true.) then
                   Q = (1.0 + sqrt(1.0 + 1.0/ZINT))*(Q0/2.0)
                else
                   Q = 16.0/106.0 
              endif
             
             ! Define optimal fV
              if( fV_opt .eq. .true. ) then
                if (PARsup .gt. Ith) then   
                  fV=(-1.0+sqrt(1.0+1.0/ZINT))*(Q0/2.0)*rIhat/vNhat
                  elseif (PARsup .le. Ith) then
                  fV=0.01
                endif
                else
                  fV = fV_ini 
              endif
             
             ! Calculate real theta from thetahat
             theta=ThetaHat * (1.0 - fV - Q0/2.0/Q)
             ! Growth rate (d-1) of phytoplankton is a function depending on
             ! temperature, light, and nutrient   
             mu=rIhat*(1.0-fV-Q0/(2.0*Q))-zetaN*fV*vNhat
             ! END OF COMPUTING PHYTOPLANKTON GROWTH RATE.  
!Debug:
      if((k .eq. N) .and. (i .eq. 50)) then
        MPI_master_only write(stdout,*) 'NO3 = ',NO3(k), 'mu = ', mu 
      endif

             
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
            PHYT(k) = PHYT(k) + dp
             NO3(k) =  NO3(k) + dn
             ZOO(k) =  ZOO(k) + dz 
             DET(k) =  DET(k) + dd
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

          Phyt(k) = SB(k,2)

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

