#include "cppdefs.h"
#if defined BIOLOGY && defined PZEFTFE
      SUBROUTINE biology_tile (Istr,Iend,Jstr,Jend)
!
! Compute biological forcing functions
!
! In this particular implementation there are 7 state variables:
! NO3, Fe, PHYTOplankton, ZOOplankton, DET, PMU (Avg. log size), VAR (Variance of log size)
!
      USE forces
      USE NPZDCONTFE
      implicit none
      integer Istr,Iend,Jstr,Jend,IMAX,JMAX,it
#include "grid.h"
#include "ocean3d.h"
#include "ocean2d.h"
#include "diagnostics.h"
#include "scalars.h"
#include "mixing.h"
      integer, parameter :: nsink   = NumVSinkTerms
      real,    parameter :: kw      = 0.04, 
     &                      gb      = -0.05,
     &                      kc      = 0.025, !Light attenuation of chlorophyll (mg Chl m-2)
     &                      mon2sec = 2592D3,!How many seconds that one month contains
     &                      bI0=0d0
      real,    parameter :: d_per_s = 864d2   ! how many seconds in one day
      real,    parameter :: PMU0    = log(1d1)! Positivize PMU
      real,    parameter :: PMUmax  = 20.
      real,    parameter :: VARmax  = 25.
      real,    parameter :: VTR     = 0.08
      real :: ThisFlux(N, NumFluxTerms) = 0. ! [mmol m-2 s-1], upward flux is positive
      real :: muNet = 0., muNet1
      real ::    SI = 0., Fe_scav
      real :: theta, dCHL, rl, QN
      real :: dthetadl,d2thetadl2,dQNdl,d2QNdl2,d2theta_QNdl2
      real ::  Lno3 = 0.,d2muNet_QNdl2
      !real :: x(M)
      real :: Temp(N), PAR(N), PARavg 
      real :: Hz_(N)  ! Vertical distance of local grid
      real :: PMU1, VAR1, PHY1 

      integer i,j,k, iB, itrc, ii_, iflux, N_MLD
      real NO3(N),PHY(N),ZOO(N),DFE(N),CHL_(N),
     &     DET(N),FC(0:N),PHYPMU(N),PHYVAR(N),
     &     I_0,I_min,Ptot,PMU,VAR,
     &     SB(N,nsink),dSB(0:N,nsink),wSB(nsink),cff,cff1,cff2,cff6,
     &     aR,aL,cu,tf_z,gbar,INGES,EGES,RES,ZMORT,pp_ND,pp_NZ,pp_DZ,
     &     pp_ZP, mu1, aI0, KN,d2gdl2bar,dgdlbar,Pl,d2muNetdl2,
     &     dmuNetdl,d3mudl3,d4mudl4

      real            :: par_,Ptot1,NPPt,CHLt,PP_PN,PPpn
      real, parameter :: Dust_solubility = 0.02  !Aumont et al. (2003)
     &                  ,Dust_ironfrac   = 0.035 !TOM10Appendix p. 24, unit: g/g
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
        PHY(k) =min(max(t(i,j,k,nnew,iPhy1),eps), 200.)  ! Phytoplankton
        ZOO(k) =max(t(i,j,k,nnew,iZOO1),eps)  ! Zooplankton
        DET(k) =min(max(t(i,j,k,nnew,iDET1),eps), 200.)  ! Detritus
        CHL_(k)  =max(CHL(i,j,k),eps)         ! CHL
        PHYPMU(k)=min(max(t(i,j,k,nnew,iPMU1),eps),1D3)  ! Avg. log size
        PHYVAR(k)=min(max(t(i,j,k,nnew,iVAR1),eps),1D3)  ! Variance of log size
        DFE(k) =max(t(i,j,k,nnew,iFer1),Femin)  ! Dissolved iron
        Temp(k)=t(i,j,k,nnew,iTemp)  ! temperature; [deg. C]

      enddo EXTRACT

#  ifdef DEBUG_BIO2
      ALOHA=.FALSE.
       if (lonr(i,j) .gt. -159 .and. lonr(i,j) .le. -157 
     &    .and. latr(i,j) .gt. 22.5 .and. latr(i,j) .le. 24) ALOHA=.TRUE.
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
#if defined IRON
! Iron atmospheric deposition:
! Dust deposition unit: kg/m2/month. Assume 0.035g Fe per g of dust
! so need to convert into nM at each time step
! Deposition = Dust*10^12*solubility*Fe%/56*dt/mon2sec*surface_area/surface_grid_volume/1000
         
       cff= Dust_solubility*dust(i,j)*1D9*Dust_ironfrac/55.85
     &    /abs(Hz_(N))/mon2sec*dt                     

       DFE(N) = DFE(N) + cff
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
        ! Convert to real PMU and VAR:
        PMU=PHYPMU(k)/PHY(k) - PMU0         
        VAR=PHYVAR(k)/PHY(k)
        PMU=min(PMU, PMUmax)
        VAR=min(max(VAR,eps), VARmax)

        call PHY_NPZDCONT(NO3(k),par_,Temp(k),DFE(k),PMU,
     &            muNet,dmuNetdl,d2muNetdl2,d3mudl3,d4mudl4,
     &            theta,QN,dthetadL,dQNdL,d2thetadL2,d2QNdL2)
#  ifdef DEBUG_BIO2
        if (ALOHA .and. k.eq.N) then
           write(stdout,*) 'NO3 = ', NO3(k) 
           write(stdout,*) 'par = ', par_
           write(stdout,*) 'Temp = ', Temp(k) 
           write(stdout,*) 'DFe = ', DFe(k) 
           write(stdout,*) 'PMU = ', PMU
           write(stdout,*) 'muNet = ', muNet
           write(stdout,*) 'dmuNetdl = ', dmuNetdl
           write(stdout,*) 'd2muNetdl2 = ', d2muNetdl2
           write(stdout,*) 'd3mudl3 = ', d3mudl3 
           write(stdout,*) 'd4mudl4 = ', d4mudl4 
           write(stdout,*) 'theta = ', theta 
           write(stdout,*) 'dthetadl = ', dthetadl 
           write(stdout,*) 'd2thetadl2 = ', d2thetadl2 
           write(stdout,*) 'dQNdL = ', dQNdL 
           write(stdout,*) 'QN = ', QN 
        endif
#  endif
 
       d2theta_QNdl2 = d2Y_Xdl2(theta,QN,dthetadl, dQNdL, d2thetadL2, d2QNdl2) 
       Qbulk(i,j,k)  = QN+0.5*VAR*d2QNdl2
       CHL(i,j,k)    = PHY(k)*(theta/QN+0.5*VAR*d2theta_QNdl2)

       d2muNet_QNdl2 = d2Y_Xdl2(muNet,QN,dmuNetdl, dQNdL, d2muNetdL2, d2QNdl2) 
       NPPt          = PHY(k)*(muNet/QN+0.5*VAR*d2muNet_QNdl2) 
!      !ZOOplankton section:
        tf_z = TEMPBOL(Ez,Temp(k))
         
        ! Calculate the integration  
        !x(1)  = PMU-2D0*sqrt(VAR)  ! Minimal size
        !x(M)  = PMU+2D0*sqrt(VAR)  ! Maximal size
        !dx    = (x(M) - x(1))/float(M-1)

        !do ii_  = 2,M-1
        !  x(ii_) = x(ii_-1) + dx  !Log sii_ze of each class
        !enddo

        !! Calculate total phytoplankton palatable prey (Ptot,Eq. 13):
        !Ptot  = 0D0
        !NPPt  = 0D0
        !Ptot1 = 0D0
        !CHLt  = 0D0

        !do ii_=1,M
        !  x(ii_) = max(min(x(ii_), 15.),  -3.)

        !  !Calculate Chl in this size class
        !  mu1=mu0  *exp(alphamu*x(ii_)+betamu*x(ii_)**2)
        !  aI0=aI0_C*exp(alphaI*x(ii_))
        !  KN =K0N  *exp(alphaK*x(ii_))

        !  call MONOD(Temp(k),par_,NO3(k),mu1,aI0,bI0,KN,
!     &                                   muNet1,QN,theta,SI,Lno3)

        !  !Probability function
        !  Pl=PHY(k)*normal(PMU,VAR,x(ii_))

        !  ! Chl in this size class:
        !  dCHL = Pl*dx*theta/QN

        !  ! Calculate total Chl:
        !  CHLt = CHLt + dCHL

        !  NPPt = NPPt + Pl*dx*muNet1/QN

        !  !Feeding preference
        !  rl   = exp(gb*x(i))
        !  Ptot = Ptot +Pl*dx*rl   ! Total palatable prey

        ! ! Considering killing-the-winner strategy
        !  Ptot1= Ptot1+Pl**alphaG*dx*rl
        !enddo

        ! save total NPP (carbon-based), units changed to d-1
        ThisFlux(k,oPPt) = NPPt * 12d0

        ! The grazing dependence on total prey
        gbar = grazing(grazing_formulation,Kp,PHY(k))

         !Zooplankton per capita total ingestion rate
         INGES = tf_z*gmax*gbar

         !Zooplankton excretion rate (-> DOM)
         RES = INGES*(1D0-GGE-unass)

         !ZOOPLANKTON EGESTION (-> POM)
         EGES = INGES*unass

!       Grazing rate on the mean size of PHY (specific to N-based Phy biomass, unit: d-1) (Eq. 12)
         !rl   = exp(gb * PMU)

         ! The probability density at avg. size
         !Pl   = normal(PMU,VAR,PMU)*PHY(k)
         ! Grazing rate on the mean size
      
         !gbar = INGES*ZOO(k)*rl/max(Ptot1,eps) 
         !gbar = gbar*Pl**(alphaG-1D0)

!       First derivative of grazing rate evaluated at the mean size (negative)
         !dgdlbar  = gb*gbar 
         dgdlbar  = 0.
         !d2gdl2bar= ((1D0-alphaG)/VAR+gb**2)*gbar
         d2gdl2bar= 0.
!      !End of zooplankton section
!=============================================================
!! Solve ODE functions:
       Zmort = ZOO(k)**2*mz*tf_z  !Mortality term for ZOO

       VAR1  = VAR+dtdays*(VAR*(VAR*(d2muNetdl2-d2gdl2bar+VTR*d4mudl4)
     &                              -5d0*VTR*d2muNetdl2)+2d0*VTR*muNet)
  
       Thisflux(k, omuNet)      = muNet
       Thisflux(k, odmuNetdl)   = dmuNetdl
       Thisflux(k, od2mudl2)    = d2muNetdl2
       Thisflux(k, od3mudl3)    = d3mudl3
       Thisflux(k, od4mudl4)    = d4mudl4
       
       !Eq. 18, Update PMU:
       PMU = PMU + PMU0  ! Restore to positive values
       PMU1= PMU + dtdays*(VAR*(dmuNetdl-dgdlbar+VTR*d3mudl3) 
     &                                              - 3d0*VTR*dmuNetdl)
#  ifdef DEBUG_BIO2
       if (ALOHA .and. k .eq. N) then
          write(stdout,*) 'PMU = ',PMU
          write(stdout,*) 'PMU1 = ',PMU1
          write(stdout,*) 'VAR = ',VAR
          write(stdout,*) 'VAR1 = ',VAR1
       endif
#  endif
       !Contrain the PMU and VAR: 
       PMU1 = min(max(PMU1,eps),PMUmax)
       VAR1 = min(max(VAR1,eps),VARmax)
          
       ! For production/destruction matrix:
       PP_ND= RDN*DET(k)*tf_z   
       PP_NZ= ZOO(k)*RES        
       PP_DZ= ZOO(k)*EGES+Zmort 
       PP_ZP= ZOO(k)*INGES      
       PPpn = PHY(k)*(muNet+0.5*VAR*(d2muNetdl2+VTR*d4mudl4)-1.5*VTR
     &                                                      *d2muNetdl2)
       PP_PN= max(PPpn,0D0)
       DET(k)    = DET(k) + dtdays*(PP_DZ         - PP_ND)
       NO3(k)    = NO3(k) + dtdays*(PP_ND + PP_NZ - PP_PN)
       PHY1      = PHY(k) + dtdays*(PP_PN         - PP_ZP)
       ZOO(k)    = ZOO(k) + dtdays*(PP_ZP - PP_NZ - PP_DZ)
       PHYPMU(k) = PHYPMU(k) +PHY(k)*(PMU1-PMU)+PMU*(PHY1-PHY(k))
       PHYVAR(k) = PHYVAR(k) +PHY(k)*(VAR1-VAR)+VAR*(PHY1-PHY(k))
       PHY(k)    = PHY1
#ifdef IRON
      ! Calculate dissolved Iron:    
#  ifdef DEBUG_BIO2
       if (ALOHA .and. (k.eq.N)) then
          write(stdout,*) 'DFE = ', DFE(k)
          write(stdout,*) 'DET = ', DET(k)
          write(stdout,*) 'PP_ND = ', PP_ND
          write(stdout,*) 'PP_NZ = ', PP_NZ
          write(stdout,*) 'PP_PN = ', PP_PN
          write(stdout,*) 'PP_ZP = ', PP_ZP
          write(stdout,*) 'PP_DZ = ', PP_DZ
       endif
#   endif

       call IRONCYCLE(Temp(k), DET(k),DFE(k),PP_ND,PP_NZ,PP_PN,Fe_scav)

#   ifdef DEBUG_BIO2
       if (ALOHA .and. (k.eq.N)) then
          write(stdout,*) 'Iron scavenge = ', Fe_scav
          write(stdout,*) 'After Iron cycle, DFe = ', DFe(k)
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
        wSB(1) =WDET   ! Only detritus sink
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
          dSB(0,iB)=0d0   ! Set no-flux boundary conditions at bottom.
!
! Apply fluxes:

          do k=N,1,-1     !<-- irreversible
            SB(k,iB)=SB(k,iB)+wSB(iB)*FC(k)*(dSB(k,iB)-dSB(k-1,iB))
            SB(k,iB)=max(SB(k,iB),eps)  !Positivity
          enddo
        enddo  ! <-- iB

        
        do k=1,N
          DET(k)=SB(k,1)    ! Register new DET
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
          t(i,j,k,nnew,iNO3_)=max(NO3(k),eps)
#ifdef IRON
          t(i,j,k,nnew,iFer1)=max(DFE(k),Femin)
#endif
          t(i,j,k,nnew,iPhy1)=max(PHY(k),eps)

          t(i,j,k,nnew,iZoo1)=max(ZOO(k),eps)
          t(i,j,k,nnew,iDET1)=max(DET(k),eps)
          t(i,j,k,nnew,iPmu1)=max(PHYPMU(k),eps)
          t(i,j,k,nnew,iVar1)=max(PHYVAR(k),eps)
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
 1001    format('NO3=',1x,F6.2,1x,'PHY=',
     &                 1x,F6.2,1x,'ZOO=',
     &                 1x,F6.2,1x,'DFE=',
     &                 1x,F6.2,1x,'at Depth',1x,F6.1)

#else
      subroutine biology_empty ()
#endif
      return
      END subroutine
