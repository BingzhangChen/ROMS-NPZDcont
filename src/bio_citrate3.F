#include "cppdefs.h"
#if defined BIOLOGY && defined PZEFTFE
      SUBROUTINE biology_tile (Istr,Iend,Jstr,Jend)
!
! Compute biological forcing functions
!
! In this particular implementation there are 7 state variables:
! NO3, Fe, PHYTOplankton, ZOOplankton, DET, PMU (Avg. log size), VAR (Variance of log size)
!
      USE forces,     only: dust,srflx
      USE NPZDCONTFE, only: thetamin
      USE NPZDCONTFE, only: dtdays,Fe_N,Femin,IRONCYCLE,citrate_PHY3
      USE NPZDCONTFE, only: gmax,grazing,eps,WDET,d2Y_Xdl2,TEMPBOL,mz
      USE NPZDCONTFE, only: RDN,Ez,KP,grazing_formulation,Calculate_PAR
      implicit none
      integer Istr,Iend,Jstr,Jend,IMAX,JMAX,it
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
     &                      mon2sec = 2592D3  !How many seconds that one month contains
      real,    parameter :: d_per_s = 864d2   ! how many seconds in one day
      real,    parameter :: PMU0    = log(1d1)! Positivize PMU
      real,    parameter :: PMUmax  = 15.
      real,    parameter :: VARmax  = 25.
      real,    parameter :: MTomax  = 40.
      real,    parameter :: VTomax  = 4D2
      real,    parameter :: MIomax  = log(1D3)
      real,    parameter :: VIomax  = 15.
      real,    parameter :: nu_L    = 0.1
      real,    parameter :: nu_T    = 0.1
      real,    parameter :: nu_I    = 0.1
      real :: ThisFlux(N, NumFluxTerms) = 0. ! [mmol m-2 s-1], upward flux is positive
      real :: muNet = 0.
      real :: Fe_scav
      real :: theta, QN
      real :: dthetadl,d2thetadl2,dQNdl,d2QNdL2,d2theta_QNdl2
      real :: Temp(N), PAR(N), PARavg 
      real :: Hz_(N)  ! Vertical distance of local grid
      real :: PMU1, VAR1, PHY1, PMU, VAR
      real :: VIo1, VTo1, MIo1, MTo1

      integer i,j,k, iB, itrc, iflux, N_MLD
      real NO3(N),PHY(N),ZOO(N),DFE(N),DETFe(N),CHL_(N),
     &     DET(N),FC(0:N),PPMU(N),PVAR(N),PMTo(N),PVTo(N),PMIo(N),PVIo(N),
     &     I_0,I_min,
     &     SB(N,nsink),dSB(0:N,nsink),wSB(nsink),cff,cff1,cff2,cff6,
     &     aR,aL,cu,tf_z,gbar,INGES,EGES,RES,ZMORT,PP_ND,PP_NZ,PP_DZ,
     &     PP_ZP,
     &     dmudL,d2mudL2, dmudx, d2mudx2, dMudZ, d2mudZ2
      real :: MTo, VTo, MIo, VIo
      real :: dthetadx, d2thetadx2, dthetadZ, d2thetadZ2
      real :: d2theta_QNdx2
      real :: d2theta_QNdZ2
      real :: d2muNet_QNdl2
      real :: d2muNet_QNdx2
      real :: d2muNet_QNdZ2
      real            :: par_,PP_PN
      real, parameter :: GGE=0.3, unass  = 0.24
      logical         :: hel_found       = .FALSE., ALOHA = .FALSE.
# include "compute_auxiliary_bounds.h"

     !the fraction of a time step in one day
      dtdays = dt/d_per_s
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
        bioFlux(i,j,:,:) =0.0
#endif

#if defined MASKING
! Do not calculate on landmask points:
        if(rmask(i,j) .LT. 0.5) then 
           t(i,j,:,nnew,:)=0.          
        else 
!
! Extract biological variables from tracer arrays; place them into
! scratch variables; restrict their values to be positive definite.
! t: tracer array!

      EXTRACT: do k=1,N      ! N: the number of vertical layers
                  do it = 5,NT
                    !Check if there are NaN values:
                    if (t(i,j,k,nnew,it) .ne. t(i,j,k,nnew,it)) then
                     write(stdout,99) it, -z_r(i,j,k)
                     write(stdout,111) ' Lon  ='
     &                  ,lonr(i,j), 'Lat = ', latr(i,j)
                     write(stdout,102) t(i,j,k,nnew,iNO3_),
     &                                 t(i,j,k,nnew,iPHY1),
     &                                 t(i,j,k,nnew,iZoo1)

                     stop
                    endif
                  enddo
        NO3(k) =max(t(i,j,k,nnew,iNO3_), 1D-6)  ! Nitrate
      
        !Calculate total CHL:
        PHY(k) =min(max(t(i,j,k,nnew,iPhy1),eps), 200.)  ! Phytoplankton
        ZOO(k) =max(t(i,j,k,nnew,iZOO1),eps)  ! Zooplankton
        DET(k) =min(max(t(i,j,k,nnew,iDET1),eps), 200.)  ! Detritus
        CHL_(k)=max(CHL(i,j,k),eps)         ! CHL
        PPMU(k)=max(t(i,j,k,nnew,iPMU1),eps)  ! Avg. log size
        PVAR(k)=max(t(i,j,k,nnew,iVAR1),eps)  ! Variance of log size
        PMTo(k)=max(t(i,j,k,nnew,iMTo1),eps)  ! Avg. mean Topt
        PVTo(k)=max(t(i,j,k,nnew,iVTo1),eps)  ! VAR. of Topt
        PMIo(k)=max(t(i,j,k,nnew,iMIo1),eps)  ! Avg. mean Iopt
        PVIo(k)=max(t(i,j,k,nnew,iVIo1),eps)  ! VAR. of Iopt
        DFE(k) =max(t(i,j,k,nnew,iFer1),Femin)   ! Dissolved iron
        DETFe(k)=max(t(i,j,k,nnew,iDETFe1),eps)  ! Detrital iron
        Temp(k) =t(i,j,k,nnew,iTemp)  ! temperature; [deg. C]

      enddo EXTRACT

#  ifdef DEBUG_BIO2
      ALOHA=.FALSE.
       if (lonr(i,j) .gt. -159 .and. lonr(i,j) .le. -157 
     &    .and. latr(i,j) .gt. 22.5 .and. latr(i,j) .le. 24) ALOHA=.TRUE.
       if (ALOHA) then
           write(6,112) 'NO3 = ',NO3(N)
           write(6,112) 'DFE = ',DFE(N)
           write(6,112) 'DETFe = ',DETFe(N)
           write(6,112) 'PHY = ',PHY(N)
           write(6,112) 'PPMU = ',PPMU(N)
           write(6,112) 'PVAR = ',PVAR(N)
           write(6,112) 'PMTo = ',PMTo(N)
           write(6,112) 'PVTo = ',PVTo(N)
           write(6,112) 'PMIo = ',PMIo(N)
           write(6,112) 'PVIo = ',PVIo(N)
           write(6,112) 'Temp = ',Temp(N)
       endif
#  endif

! Calulate PAR: Set Photosynthetically Available Radiation (PAR) at
! surface from solar radiation x 0.43.
      I_0 = srflx(i,j)*rho0*Cp*0.43    ! Surface short wavelength radiation
                                       ! and convert to the unit W m-2 
      I_min = 0.01*I_0    ! Light of bottom of euphotic layer

! convert back from Celsius meter second-1 to Watts meter-2
! srflx obtained from the Routine set_forces.F
! rho0: average density of seawater
! Cp  : specific heat of seawater

      ! Update vertical light level:
      ! Obtain local Hz:
      Hz_(:) = Hz(i,j,:)

      call Calculate_PAR(I_0, N, Hz_, CHL_(:), PAR)

#if defined LMD_SKPP || defined LMD_BKPP
      ! Calculate the vertical grid index (N_MLD) at the bottom of MLD:
      ! kbl : Index of first grid level below surface boundary layer.
      N_MLD = min(kbl(i,j) + 1, N)  ! For calculating average PAR within MLD

   !  Calculate average PAR within the surface mixed layer (from nlev to N_MLD):
      PARavg=0d0
      do k=N,N_MLD,-1
         PARavg = PARavg+PAR(k)*Hz_(k)
      enddo 
      PARavg = PARavg/abs(Z_w(i,j,N_MLD-1))
#   ifdef DEBUG_BIO2
      if (ALOHA) then
         write(stdout,*)  'MLD = ', abs(Z_w(i,j,N_MLD-1))
         write(stdout,*)  'PARavg = ', PARavg
      endif 
#   endif
#endif
!
! Iron atmospheric deposition:
! Soluble iron deposition unit: kg/m2/s.
! so need to convert into nM at each time step
! Deposition = Dust*10^12/56*dtsec*surface_area/surface_grid_volume
       cff      = dust(i,j)*1D9/55.85/Hz_(N)*dt
#   ifdef DEBUG_BIO2
       if (ALOHA) write(stdout,112) 'Iron depo. d-1 = ', cff/dt*d_per_s
#   endif
       DFE(N)   = DFE(N) + cff

! Determine the depth of Euphotic layer:
      hel_found = .FALSE.

      DO k=N,1,-1     !<-- irreversible

#if defined LMD_SKPP || defined LMD_BKPP
        ! Check whether in the MLD or not
        if (k .lt. N_MLD) then
           par_ = PAR(k)
        else
           par_ = PARavg
        endif
#else
        par_    = PAR(k)
#endif
        IF (par_ < I_0*1D-3 .or. z_r(i,j,k) < -250.) THEN
         muNet  =0.
         dmudL  =0.
         d2mudL2=0.
         dmudx  =0.
         d2mudx2=0.
         dmudZ  =0.
         d2mudZ2=0.
         ThisFlux(k,oPPt) = 0.
        ELSE
         ! Convert to real PMU and VAR:
         PMU=PPMU(k)/PHY(k)
         VAR=PVAR(k)/PHY(k) - PMU**2
         PMU=PMU - PMU0         
         PMU=min(PMU,          PMUmax)
         VAR=min(max(VAR,1D-3), VARmax)
         MTo=PMTo(k)/PHY(k)
         VTo=PVTo(k)/PHY(k) - MTo**2
         MTo=min(MTo, MTomax)
         VTo=min(max(VTo,1D-3), VTomax)
         MIo=PMIo(k)/PHY(k)
         VIo=PVIo(k)/PHY(k) - MIo**2
         MIo=min(MIo, MIomax)
         VIo=min(max(VIo,1D-3), VIomax)

         call citrate_PHY3(NO3(k),par_,Temp(k),DFE(k), 
     &           PMU,MIo,MTo, muNet,
     &           dmudL,d2mudL2, dmudx, d2mudx2, dMudZ, d2mudZ2,
     &           theta,QN,dthetadl,d2thetadL2,dQNdL, d2QNdL2, 
     &           dthetadZ, d2thetadZ2,
     &           dthetadx,d2thetadx2)
         
#  ifdef DEBUG_BIO2
        if (ALOHA .and. k.eq.N) then
           write(stdout,112) 'par_ = ', par_
           write(stdout,112) 'PMU = ', PMU
           write(stdout,112) 'MTo = ', MTo
           write(stdout,112) 'VTo = ', VTo
           write(stdout,112) 'MIo = ', MIo
           write(stdout,112) 'VIo = ', VIo
           write(stdout,112) 'muNet = ', muNet
           write(stdout,112) 'dmudL = ', dmudL
           write(stdout,112) 'd2mudL2 = ', d2mudL2
           write(stdout,112) 'dmudx = ', dmudx
           write(stdout,112) 'd2mudx2 = ', d2mudx2
           write(stdout,112) 'dmudZ = ', dmudZ
           write(stdout,112) 'd2mudZ2 = ', d2mudZ2
           write(stdout,112) 'theta = ', theta 
           write(stdout,112) 'dthetadl = ', dthetadl 
           write(stdout,112) 'd2thetadl2 = ', d2thetadl2 
           write(stdout,112) 'dQNdL = ', dQNdL 
           write(stdout,112) 'QN = ', QN 
        endif
#  endif
        d2theta_QNdl2 = d2Y_Xdl2(theta,QN,dthetadl, dQNdL,
     &                           d2thetadL2, d2QNdL2) 
        d2theta_QNdx2 = d2thetadx2/QN**2
        d2theta_QNdZ2 = d2thetadZ2/QN**2

        Qbulk(i,j,k)  = QN + 0.5*VAR*d2QNdL2

        !Chl:C ratio of the whole community:
        cff1 = theta/QN + .5*VAR*d2theta_QNdl2
     &                  + .5*VIo*d2theta_QNdx2
     &                  + .5*VTo*d2theta_QNdZ2
        cff1 = max(cff1,thetamin)
        CHL(i,j,k) = PHY(k)*cff1

        d2muNet_QNdl2 = d2Y_Xdl2(muNet,QN,dmudL, dQNdL,d2mudL2,d2QNdL2) 
        d2muNet_QNdx2 = d2mudx2/QN**2
        d2muNet_QNdZ2 = d2mudZ2/QN**2
        ! Primary production (ug C L-1 d-1)
        cff1 = 12.*PHY(k)*(muNet/QN + 0.5*VAR*d2muNet_QNdl2 
     &                              + 0.5*VIo*d2muNet_QNdx2
     &                              + 0.5*VTo*d2muNet_QNdZ2) 
        ThisFlux(k,oPPt)=max(cff1,0.)
#   ifdef DEBUG_BIO2
        if (ALOHA .and. k .eq. N) then
           write(stdout,112) 'N:C = ', Qbulk(i,j,k)
           write(stdout,112) 'Chl:C = ', CHL(i,j,k)/PHY(k)
           write(stdout,112) 'NPP = ', ThisFlux(k,oPPt)
        endif
#   endif

       ENDIF ! <== PAR_ < I_0/1000
!      ZOOplankton section:
        tf_z = TEMPBOL(Ez,Temp(k))
         
        ! The grazing dependence on total prey
        gbar = grazing(grazing_formulation,KP,PHY(k))

        !Zooplankton per capita total ingestion rate
        INGES = tf_z*gmax*gbar

        !Zooplankton excretion rate (-> DOM)
        RES = INGES*(1.-GGE-unass)

        !ZOOPLANKTON EGESTION (-> POM)
        EGES = INGES*unass
!      !End of zooplankton section
!=============================================================
!! Solve ODE functions:
       Zmort = ZOO(k)*ZOO(k)*mz*tf_z  !Mortality term for ZOO

       ! Update trait variance:
       VAR1  = VAR+dtdays*(VAR*VAR*d2mudL2 + 2.*nu_L*muNet)
       VTo1  = VTo+dtdays*(VTo*VTo*d2mudZ2 + 2.*nu_T*muNet)
       VIo1  = VIo+dtdays*(VIo*VIo*d2mudx2 + 2.*nu_I*muNet)
  
       Thisflux(k, omuNet)    = muNet
       Thisflux(k, odmuNetdl) = dmudL
       Thisflux(k, od2mudl2)  = d2mudL2
       Thisflux(k, odmudx)    = dmudx
       Thisflux(k, od2mudx2)  = d2mudx2
       Thisflux(k, odmudZ)    = dmudZ
       Thisflux(k, od2mudZ2)  = d2mudZ2
       
       !Update mean trait:
       PMU = PMU + PMU0  ! Restore to positive values
       PMU1= PMU + dtdays*VAR*dmudL 
       MTo1= MTo + dtdays*VTo*dmudZ
       MIo1= MIo + dtdays*VIo*dmudx
#  ifdef DEBUG_BIO2
       if (ALOHA .and. k .eq. N) then
          write(stdout,112) 'VAR = ', VAR
          write(stdout,112) 'VAR1 = ', VAR1
          write(stdout,112) 'MTo1 = ',MTo1
          write(stdout,112) 'VTo = ',VTo
          write(stdout,112) 'VTo1 = ',VTo1
          write(stdout,112) 'VIo = ',VIo
          write(stdout,112) 'VIo1 = ',VIo1
       endif
#  endif
       ! For production/destruction matrix:
       PP_ND= RDN*DET(k)*tf_z   
       PP_NZ= ZOO(k)*RES        
       PP_DZ= ZOO(k)*EGES+Zmort 
       PP_ZP= ZOO(k)*INGES      
       PP_PN= PHY(k)*(muNet+.5*(VAR*d2mudL2+VTo*d2mudZ2+VIo*d2mudx2))
       PP_PN= max(PP_PN,0.)
       DET(k)  = DET(k) + dtdays*(PP_DZ         - PP_ND)
       NO3(k)  = NO3(k) + dtdays*(PP_ND + PP_NZ - PP_PN)
       PHY1    = PHY(k) + dtdays*(PP_PN         - PP_ZP)
       ZOO(k)  = ZOO(k) + dtdays*(PP_ZP - PP_NZ - PP_DZ)
       PPMU(k) = PPMU(k)+ PHY(k)*(PMU1-PMU)+PMU*(PHY1-PHY(k))
       PMTo(k) = PMTo(k)+ PHY(k)*(MTo1-MTo)+MTo*(PHY1-PHY(k))
       PMIo(k) = PMIo(k)+ PHY(k)*(MIo1-MIo)+MIo*(PHY1-PHY(k))

       PVAR(k) = PVAR(k)+ PHY(k)*(VAR1-VAR+2.*PMU*(PMU1-PMU))
     &         + (VAR+PMU**2)*(PHY1-PHY(k))

       PVIo(k) = PVIo(k)+ PHY(k)*(VIo1-VIo+2.*MIo*(MIo1-MIo))
     &         + (VIo+MIo**2)*(PHY1-PHY(k))

       PVTo(k) = PVTo(k)+ PHY(k)*(VTo1-VTo+2.*MTo*(MTo1-MTo))
     &         + (VTo+MTo**2)*(PHY1-PHY(k))
       PHY(k)  = PHY1
#ifdef IRON
      ! Calculate dissolved Iron:    
#  ifdef DEBUG_BIO2
       if (ALOHA .and. (k.eq.N)) then
          write(stdout,112) 'DFE = ', DFE(k)
          write(stdout,112) 'DET = ', DET(k)
          write(stdout,112) 'PP_ND = ', PP_ND
          write(stdout,112) 'PP_NZ = ', PP_NZ
          write(stdout,112) 'PP_PN = ', PP_PN
          write(stdout,112) 'PP_ZP = ', PP_ZP
          write(stdout,112) 'PP_DZ = ', PP_DZ
       endif
#   endif

       call IRONCYCLE(Temp(k), DET(k),PP_NZ,PP_ND,PP_PN,PP_DZ,Fe_scav,
     &                DETFe(k),DFE(k))

#   ifdef DEBUG_BIO2
       if (ALOHA .and. (k.eq.N)) then
          write(stdout,112) 'Iron scavenge = ', Fe_scav
          write(stdout,112) 'After Iron cycle, DFe = ', DFe(k)
       endif
#   endif
#endif

! Calculate euphotic depth
        IF (.NOT. hel_found) THEN
           if (PAR(k) < I_min) then
               hel(i,j)  = -z_w(i,j,k-1)
               hel_found = .TRUE.
           endif
        ENDIF
      ENDDO
! ---------------------------------------------------------------------------
! Vertical sinking: Vertical advection algorithm based on monotonic,
! continuous conservative parabolic splines.
! 1: Put biological tracers into SB array.
      do k=1,N
        SB(k,1)=DET(k)
        SB(k,2)=DETFe(k)
        wSB(1) =WDET   ! Only detritus sink
        wSB(2) =WDET   ! Only detritus sink
      enddo

      do iB=1,nsink

! Part (i): Construct parabolic splines: compute vertical derivatives
! of the fields SB. The derivatives are located at W-points;
! Neumann boundary conditions are assumed on top and bottom.

! Neumann boundary condition: the 1st derivatives at the boundaries are constants.
! SB: biological tracers
! Hz: vertical level thickness
! FC: 
! dSB: Vertical derivatives of fields SB at the grid box interfaces

          dSB(0,iB)=0d0
          FC(0)=0d0
          cff6=6d0
          do k=1,N-1
            cff=1d0/(2d0*Hz(i,j,k+1)+Hz(i,j,k)*(2d0-FC(k-1)))
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
          cff      =1d0/3d0
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
            FC(k)=dtdays/Hz(i,j,k)
            aR   =dSB(k,iB)
            aL   =dSB(k-1,iB)
            cff1 =(aR-aL)*6d0*(SB(k,iB)-0.5*(aR+aL))
            cff2 =(aR-aL)**2
            if ((aR-SB(k,iB))*(SB(k,iB)-aL).lt.0.) then
              aL=SB(k,iB)
              aR=SB(k,iB)
            elseif (cff1.gt.cff2) then
              aL=3d0*SB(k,iB)-2d0*aR
            elseif (cff1.lt.-cff2) then
              aR=3d0*SB(k,iB)-2d0*aL
            endif
            cu         =wSB(iB)*FC(k)
            dSB(k-1,iB)=SB(k,iB)-(1d0-cu)*(0.5d0*(aR-aL)-(0.5d0*(aR+aL)
     &                                   -SB(k,iB) )*(1d0-2d0*cu))
          enddo
          dSB(N,iB)=0d0   ! Set no-flux boundary conditions at top.
!
! Apply fluxes:

          do k=N,1,-1     !<-- irreversible
            SB(k,iB)=SB(k,iB)+wSB(iB)*FC(k)*(dSB(k,iB)-dSB(k-1,iB))
            SB(k,iB)=max(SB(k,iB),eps)  !Positivity
          enddo

          if(iB==1) then
            !The flux of DET sunk out of the water column
            cu    =WDET*dtdays/Hz(i,j,1)*dSB(0,1)
            !DET sunk out of the water column are replenished by NO3 immediately
            NO3(1)=NO3(1)+cu
            !DFe is replenished by DET times Fe/N ratio
            DFE(1)=DFE(1)+cu*Fe_N*TEMPBOL(Ez, TEMP(1))
          endif
        enddo  ! <-- iB
#   ifdef DEBUG_BIO2
        if (ALOHA) then
           write(stdout,*) 'After sinking, bottome values'
           write(stdout,112) 'NO3 = ', NO3(1)
           write(stdout,112) 'DFE = ', DFE(1)
           stop
        endif
#   endif
        do k=1,N
          DET(k)  =SB(k,1)    ! Register new DET
          DETFe(k)=SB(k,2)
        enddo

#ifdef DIAGNOSTICS_BIO
! Transfer fluxes to global arrays at the end of each biological time step
        do k = 1, N
          do iflux = 1, NumFluxTerms
            bioFlux(i,j,k,iflux) = ThisFlux(k, iflux)
          end do
        end do
#endif
!
! Write back the data into t array:
        do k=1,N
          t(i,j,k,nnew,iNO3_)  =max(NO3(k),eps)
          t(i,j,k,nnew,iFer1)  =max(DFE(k),Femin)
          t(i,j,k,nnew,iDETFe1)=max(DETFe(k),eps)  ! Detrital iron
          t(i,j,k,nnew,iPhy1)=max(PHY(k),eps)
          t(i,j,k,nnew,iZoo1)=max(ZOO(k),eps)
          t(i,j,k,nnew,iDET1)=max(DET(k),eps)
          t(i,j,k,nnew,iPMU1)=max(PPMU(k),eps)
          t(i,j,k,nnew,iVAR1)=max(PVAR(k),eps)
          t(i,j,k,nnew,iVTo1)=max(PVTo(k),eps)
          t(i,j,k,nnew,iVIo1)=max(PVIo(k),eps)
          t(i,j,k,nnew,iMIo1)=max(PMIo(k),eps)
          t(i,j,k,nnew,iMTo1)=max(PMTo(k),eps)
        enddo
      endif   ! <== The end of rmasking
#endif /* MASKING */
        enddo
      enddo
102      format('NO3 = ',E15.3,' PHY = ',E20.3,' MIC =',E20.3,
     &          ' MES = ',E20.3,'  Fe = ',E20.3)

   98    format('Error: Abnormal values occur in',1x,a,
     &                 1x,'at Depth',1x,F6.1)

   99    format('Error: NaN values occur in',1x,I2,
     &                 1x,'at Depth',1x,F6.1)

 111     format('2(A6,F12.3)')       
 112     format(A10,2x, 1pe12.2)       
 1001    format('NO3=',1x,F6.2,1x,'PHY=',
     &                 1x,F6.2,1x,'ZOO=',
     &                 1x,F6.2,1x,'DFE=',
     &                 1x,F6.2,1x,'at Depth',1x,F6.1)

#else
      subroutine biology_empty ()
#endif
      return
      END subroutine
