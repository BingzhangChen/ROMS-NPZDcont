#include "cppdefs.h"
#if defined BIOLOGY && defined BIO_NPZDDISC
      SUBROUTINE biology_tile (Istr,Iend,Jstr,Jend)
!
! Compute biological forcing functions
!
! In this particular implementation there are 7 state variables:
! NO3, Fe, PHYTOplankton (20 size classes), ZOOplankton, DET, DETFe
!
      USE forces,     only: dust,srflx
      USE NPZDCONTFE, only: thetamin,VTR,alphaG
      USE NPZDCONTFE, only: dtdays,Fe_N,Femin,IRONCYCLE,MONOD_size
      USE NPZDCONTFE, only: gmax,grazing,eps,WDET,TEMPBOL,mz
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
     &                      kc      = 0.025   !Light attenuation of chlorophyll (mg Chl m-2)
      real,    parameter :: d_per_s = 864d2   ! how many seconds in one day
      real :: ThisFlux(N, NumFluxTerms) = 0. ! [mmol m-2 s-1], upward flux is positive
      real :: muNet(NPHY) = 0.
      real ::    SI(NPHY) = 0., Fe_scav = 0.
      real :: theta(NPHY), QN(NPHY),Graz(NPHY)
      real ::  Lno3(NPHY) = 0.
      real :: Temp(N), PAR(N), PARavg 
      real :: Hz_(N)  ! Vertical distance of local grid

      integer i,j,k, iB, itrc, l, iflux, N_MLD
      real NO3(N),PHY(N,NPHY),ZOO(N),DFE(N),DETFe(N),CHL_(N),
     &     DET(N),FC(0:N),
     &     I_0,I_min,PHYtot,PHYC,CHLt,PHYtot2,
     &     SB(N,nsink),dSB(0:N,nsink),wSB(nsink),cff,cff1,cff2,cff6,
     &     aR,aL,cu,tf_z,gbar,INGES,EGES,RES,Zmort,pp_ND,pp_NZ,pp_DZ,
     &     pp_ZP

      real            :: par_,PP_PN,ppC_PN
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
        NO3(k) =t(i,j,k,nnew,iNO3_)  ! Nitrate
      
        !Calculate total CHL:
        do l=1,NPHY
           PHY(k,l)=max(t(i,j,k,nnew,iPhy1(l)),eps)  ! Phytoplankton
        enddo
        ZOO(k)   =max(t(i,j,k,nnew,iZOO1),eps)  ! Zooplankton
        DET(k)   =max(t(i,j,k,nnew,iDET1),eps)  ! Detritus
        CHL_(k)  =max(CHL(i,j,k),eps)         ! CHL
        DFE(k)   =max(t(i,j,k,nnew,iFer1),Femin)  ! Dissolved iron
        DETFe(k) =max(t(i,j,k,nnew,iDETFe1),eps)  ! Detrital iron
        Temp(k)  =t(i,j,k,nnew,iTemp)  ! temperature; [deg. C]

      enddo EXTRACT

#  ifdef DEBUG_BIO2
      ALOHA=.FALSE.
      if (lonr(i,j) .gt. -159 .and. lonr(i,j) .le. -157 
     &    .and. latr(i,j) .gt. 22.5 .and. latr(i,j) .le. 24) ALOHA=.TRUE.
      if (ALOHA) then
         write(6,112) 'NO3 = ',NO3(N)
         write(6,112) 'DFE = ',DFE(N)
         write(6,112) 'DETFe= ',DETFe(N)
         write(6,112) 'CHL  = ',CHL_(N)
         do l = 1, NPHY
            write(6,112) 'PHY  = ',PHY(N,l)
         enddo
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
#endif
!
#if defined IRON
! Iron atmospheric deposition:
! Dust deposition unit: kg/m2/s. 
! so need to convert into nM at each time step
! Deposition = Dust*10^12/56*dt/surface_area/surface_grid_volume/1000
      cff    = dust(i,j)*1D9/55.85/Hz_(N)*dt
#ifdef DEBUG_BIO2
      if (ALOHA) then
         write(6,112) 'DFe before dust depo.', DFE(N)
      endif
#endif
      DFE(N) = DFE(N) + cff

#ifdef DEBUG_BIO2
      if (ALOHA) then
         write(6,112) 'DFe after dust depo.', DFE(N)
      endif
#endif

#endif
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
        PHYtot  =0d0  ! Calculate total PHY biomass (N)
        PHYC    =0d0  ! Calculate total PHY biomass (C)
        CHLt    =0d0  ! Calculate total Chl a
        PHYtot2 =0d0  ! total P**alphaG
        PP_PN   =0d0  ! N-based NPP
        ppC_PN  =0d0  ! C-based NPP
      !! Phytoplankton section:
       do l = 1, NPHY

         call MONOD_size(Temp(k),par_,NO3(k),DFE(k),PMU_(l),muNet(l),
     &                    QN(l),theta(l),SI(l),Lno3(l))

#ifdef INFLEX
         QN(l)   = 16./106.  !Fix N:C to redfield ratio
#endif
         PHYtot  = PHYtot  + PHY(k,l)
         PHYtot2 = PHYtot2 + PHY(k,l)**alphaG
         PHYC    = PHYC    + PHY(k,l)/QN(l)
         CHLt    = CHLt    + PHY(k,l)/QN(l)*theta(l)
         PP_PN   = PHY(k, l)*muNet(l)       + PP_PN
         ppC_PN  = PHY(k, l)*muNet(l)/QN(l) + ppC_PN
       enddo

       ! save total NPP (carbon-based), units changed to d-1
       ThisFlux(k,oPPt)   = ppC_PN*12.
       ThisFlux(k,oPAR)   = par_
       ThisFlux(k,otheta) = CHLt/PHYC
       Qbulk(i,j,k)       = PHYtot/PHYC
       CHL(i,j,k)         = CHLt

#ifdef DEBUG_BIO2
       if (ALOHA .and. k==N) then
          write(6, 112) 'Primary production = ', ThisFlux(k,oPPt)
          write(6, 112) 'Chl:C = ', ThisFlux(k,otheta)
          write(6, 112) '  N:C = ', Qbulk(i,j,k)
          write(6, 112) '  CHL = ', CHL(i,j,k)
       endif
#endif
!      ZOOplankton section:
       tf_z = TEMPBOL(Ez,Temp(k))
         
       ! The grazing dependence on total prey
       gbar  = grazing(grazing_formulation,Kp,PHYtot)

       !Zooplankton per capita total ingestion rate
       INGES = tf_z*gmax*gbar

       !Zooplankton excretion rate (-> DOM)
       RES = INGES*(1D0-GGE-unass)

       !ZOOPLANKTON EGESTION (-> POM)
       EGES = INGES*unass

!      !End of zooplankton section
!=============================================================
!! Solve ODE functions:
       Zmort = ZOO(k)**2*mz*tf_z  !Mortality term for ZOO
       PP_ZP = ZOO(k)*INGES      


       ! For production/destruction matrix:
       PP_ND = RDN*DET(k)
       PP_NZ = ZOO(k)*RES        
       PP_DZ = ZOO(k)*EGES+Zmort 
       DET(k)= DET(k) + dtdays*(PP_DZ         - PP_ND)
       NO3(k)= NO3(k) + dtdays*(PP_ND + PP_NZ - PP_PN)

       do l = 1, NPHY
          ! Eq. 10 in Smith et al. (2016)
          Graz(l) = (INGES*ZOO(k)/PHYtot2)*PHY(k,l)**(alphaG - 1d0)

          !Add trait diffusion
          if (l == 1) then
             PHY(k,l)= PHY(k,l)*(1.+dtdays*(muNet(l)*(1.-VTR)-Graz(l)))
     &               + PHY(k,l+1)*muNet(l+1)*dtdays*VTR        
          elseif (l == NPHY) then
             PHY(k,l)= PHY(k,l)*(1.+dtdays*(muNet(l)*(1.-VTR)-Graz(l)))
     &               + PHY(k,l-1)*muNet(l-1)*dtdays*VTR        
          else
             PHY(k,l)= PHY(k,l)*(1. + dtdays
     &               * (muNet(l)*(1.-2.*VTR)-Graz(l)))
     &               + PHY(k,l-1)*muNet(l-1)*dtdays*VTR        
     &               + PHY(k,l+1)*muNet(l+1)*dtdays*VTR        
          endif
       enddo
       ZOO(k)= ZOO(k) + dtdays*(PP_ZP - PP_NZ - PP_DZ)

#ifdef DEBUG_BIO2
       if (ALOHA .and. k==N) then
          write(6,112) 'Grazing rate = ', PP_ZP/PHYtot
          do l = 1, NPHY
            write(6,112) 'After calculation, PHY  = ',PHY(N,l)
          enddo
       endif
#endif

#ifdef IRON
      ! Calculate dissolved Iron:    
       call IRONCYCLE(Temp(k), DET(k),PP_NZ,PP_ND,PP_PN,PP_DZ,Fe_scav,
     &                DETFe(k),DFE(k))
#endif

! Calculate euphotic depth
        IF (.NOT. hel_found) THEN
           if (PAR(k) < I_min) then
               hel(i,j)  = -z_w(i,j,k-1)
               hel_found = .TRUE.
           endif
        ENDIF
      ENDDO

#ifdef DEBUG_BIO2
      if (ALOHA) then
         do k=1,N
            write(6,*) 'Before sinking'
            write(6,112) 'DET   = ',DET(k)
            write(6,112) 'DETFe = ',DETFe(k)
         enddo
      endif
#endif
! Vertical sinking: Vertical advection algorithm based on monotonic,
! continuous conservative parabolic splines.
! 1: Put biological tracers into SB array.
      do k=1,N
        SB(k,1)=DET(k)
        SB(k,2)=DETFe(k)
        wSB(1) =WDET*(0.3+0.03*TEMP(N))  ! Only detritus sink
        wSB(2) =wSB(1)   ! Sinking rate of DETFe
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
          FC(0)    =0d0
          cff6     =6d0
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
            cu    =wSB(iB)*dtdays/Hz(i,j,1)*dSB(0,1)
            !DET sunk out of the water column are replenished by NO3 immediately
            NO3(1)=NO3(1)+cu
            !DFe is replenished by DET times Fe/N ratio
            DFE(1)=DFE(1)+cu*Fe_N*TEMPBOL(Ez, TEMP(1))
          endif
        enddo  ! <-- iB

        
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

#ifdef DEBUG_BIO2
      if (ALOHA) then
         do k=1,N
            write(6,*) 'After sinking'
            write(6,112) 'DET   = ',DET(k)
            write(6,112) 'DETFe = ',DETFe(k)
         enddo
         stop
      endif
#endif

!
! Write back the data into t array:
        do k=1,N
          t(i,j,k,nnew,iNO3_)  = max(NO3(k),eps)
          t(i,j,k,nnew,iFer1)  = max(DFE(k),Femin)
          t(i,j,k,nnew,iDETFe1)= max(DETFe(k),eps)  ! Detrital iron
          t(i,j,k,nnew,iZoo1)  = max(ZOO(k),eps)
          t(i,j,k,nnew,iDET1)  = max(DET(k),eps)
          do l=1,NPHY
             t(i,j,k,nnew,iPhy1(l))=max(PHY(k,l),eps)
          enddo
        enddo
      endif   ! <== The end of rmasking
#endif /* MASKING */
        enddo
      enddo
102      format('NO3 = ',E15.3,' PHY = ',E20.3,' MIC =',E20.3,
     &          ' MES = ',E20.3,'  Fe = ',E20.3)

   98    format('Error: Abnormal values occur in',1x,a,
     &                 1x,'at Depth',1x,F6.1)

   99    format('Error: NaN values occur in',1x,a,
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
