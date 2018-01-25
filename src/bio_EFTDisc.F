#include "cppdefs.h"
#if defined BIOLOGY && defined BIO_EFTDISC

      SUBROUTINE biology_tile (Istr,Iend,Jstr,Jend)
!
! Compute biological forcing functions
!
! In this particular implementation there are 23 state variables:
! NO3, 20 PHYTOplankton size classes, ZOOplankton, DET
!
      USE forces, only: srflx
      use bio_MOD
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
     &                      mon2sec = 2592D3,!How many seconds that one month contains
     &                      eps     = 1D-20
      real,    parameter :: d_per_s = 864d2   ! how many seconds in one day
      real :: ThisFlux(N, NumFluxTerms) = 0. ! [mmol m-2 s-1], upward flux is positive
      real :: muNet(NPHY) = 0.
      real ::  Graz(NPHY) = 0.
      real ::    SI(NPHY) = 0.
      real ::  Lno3(NPHY) = 0.
      real :: Temp(N), PAR(N), PARavg 
      real :: Hz_(N)  ! Vertical distance of local grid

      integer i,j,k, iB, itrc, ii_, iflux, N_MLD
      real NO3(N),PHY(N,NPHY),ZOO(N),CHL(N),DFE(N),
     &     DET(N),FC(0:N),
     &     I_0,I_min,
     &     SB(N,nsink),dSB(0:N,nsink),wSB(nsink),cff,cff1,cff2,cff6,
     &     aR,aL,cu,tf_z,gbar,INGES,EGES,RES,ZMORT,pp_ND,pp_NZ,pp_DZ,
     &     pp_ZP

      real            :: par_,PHYtot, PHYtot2, pp_PN, ppC_PN 
      real            :: alphaG
#if defined IRON && defined IRONCYCLE
      real, parameter :: Dust_solubility = 0.02  !Aumont et al. (2003)
     &                  ,Dust_ironfrac   = 0.035 !TOM10Appendix p. 24, unit: g/g
#endif

      real, parameter    :: GGE=0.3, unass =0.24
      logical            :: hel_found      = .FALSE.
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
        CHL(k) =0d0
        do it = 1, NPHY
           PHY(k,it) =t(i,j,k, nnew, iPHY1(it))  ! Phytoplankton
           CHL(k)    =CHL(k)+ PHY(k,it)/QN(i,j,k,it)*THETA(i,j,k,it)
        enddo

        !Calculate total CHL:

        ZOO(k) =t(i,j,k,nnew,iZOO1)  ! Zooplankton
        DET(k) =t(i,j,k,nnew,iDET1)  ! Detritus
#if defined IRON && defined IRONCYCLE
        DFE(k) =t(i,j,k,nnew,iFer1)  ! Dissolved iron
#endif
        Temp(k)=t(i,j,k,nnew,iTemp)  ! temperature; [deg. C]

      enddo EXTRACT

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
      call Calculate_PAR(I_0, N, Hz_, CHL, PAR)

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

#if defined IRON && defined IRONCYCLE
! Iron atmospheric deposition:
! Dust deposition unit: kg/m2/month. Assume 0.035g Fe per g of dust
! so need to convert into nM at each time step
! Deposition = Dust*10^12*solubility*Fe%/56*dt/mon2sec*surface_area/surface_grid_volume/1000
        DFE(N) = DFE(N)
     &    + Dust_solubility*dust(i,j)*1D9*Dust_ironfrac/55.85
     &    / z_r(i,j,N)/mon2sec*dt                     

#endif

      if (kill_the_winner) then
         alphaG=1.1d0
      else
         alphaG=1d0
      endif

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

        PHYtot  =0d0  ! Calculate total PHY biomass
        PHYtot2 =0d0  ! total P**alphaG
        pp_PN   =0d0  ! N-based NPP
        ppC_PN  =0d0  ! C-based NPP
      !! Phytoplankton section:
        do ii_ = 1, NPHY
!      EFT_size(PMU, NO3, tC, par_,muNet,QN,Theta,LSI,Lno3) 
          call EFT_size(PMU_(ii_),NO3(k),Temp(k), par_, 
     &            muNet(ii_), QN(i,j,k,ii_),THETA(i,j,k,ii_),
     &            SI(ii_),    Lno3(ii_))

        ! Adjust THETA:
          if (THETA(i,j,k,ii_) .le. 0.01 .and. k .lt. N)
     &        THETA(i,j,k,ii_) = THETA(i,j,k+1,ii_)

          PHYtot  = PHYtot +PHY(k,ii_)
          PHYtot2 = PHYtot2+PHY(k,ii_)**alphaG
          pp_PN   = PHY(k, ii_)*muNet(ii_) + pp_PN
          ppC_PN  = PHY(k, ii_)*muNet(ii_)/QN(i,j,k,ii_) + ppC_PN
        enddo

        ! save total NPP (carbon-based), units changed to d-1
        ThisFlux(k,oPPt) = ppC_PN/dtdays
!---------------------------------------------------------------
        ! ZOOplankton section:
        tf_z = TEMPBOL(Ez,Temp(k))

        ! The grazing dependence on total prey (dimensionless)
        gbar = grazing(grazing_formulation,KP,PHYtot)

        !Zooplankton total ingestion rate
        INGES = tf_z*dtdays*GMAX*gbar

        !Zooplankton excretion rate (-> DOM)
        RES  = INGES*(1d0-GGE-unass)

        !ZOOPLANKTON EGESTION (-> POM)
        EGES = INGES*unass
    
        ! Grazing rate on PHY each size class (specific to N-based Phy biomass, unit: d-1) (Eq. 12)
          ! Calculate the specific grazing rate for each size class
        do ii_ = 1, NPHY
         
          ! Eq. 10 in Smith & Sergio
          Graz(ii_) = (INGES*ZOO(k)/PHYtot2)*PHY(k,ii_)**(alphaG-1d0)
          
          PHY(k,ii_) = max( PHY(k,ii_)*(1d0+muNet(ii_))
     &      /(1d0 + Graz(ii_)), eps)
        enddo
        !!End of zooplankton section
!=============================================================
        !! Solve ODE functions:
        Zmort = ZOO(k)**2*dtdays* MZ *tf_z  !Mortality term for ZOO
        
        ! For production/destruction matrix:
        pp_ND = dtdays*RDN *DET(k)*tf_z   
        pp_NZ = ZOO(k)*RES        
        pp_DZ = ZOO(k)*EGES+Zmort 
        pp_ZP = ZOO(k)*INGES      

        DET(k)= (DET(k)+pp_DZ)/(1d0 + dtdays*RDN*tf_z)
        NO3(k)= (NO3(k)+pp_ND+pp_NZ)/(1d0+pp_PN/NO3(k))
        
        ZOO(k)= (ZOO(k)+pp_ZP)/(1d0
     &        + EGES+ZOO(k)*dtdays*MZ*tf_z+RES)
            
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
#if defined IRON && defined IRONCYCLE
          t(i,j,k,nnew,iFer1)=max(DFE(k),Femin)
#endif
          do ii_=1,NPHY
             t(i,j,k,nnew,iPhy1(ii_))=max(PHY(k,ii_),eps)
          enddo

          t(i,j,k,nnew,iZoo1)=max(ZOO(k),eps)
          t(i,j,k,nnew,iDET1)=max(DET(k),eps)
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
