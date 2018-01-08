! $Id: bio_PZEFT.F 1458 2015-04-30 15:01:25Z bzchen $
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://www.romsagrif.org
!======================================================================
!
#include "cppdefs.h"
#if defined BIOLOGY && defined BIO_PZEFT

      SUBROUTINE biology_tile (Istr,Iend,Jstr,Jend)
      USE PZEFTFe
!
! Compute biological forcing functions
!
! In this particular implementation there are 7 state variables:
! NO3, Fe, PHYTOplankton, ZOOplankton, DET, Mean size, and variance.
!
      implicit none
      integer Istr,Iend,Jstr,Jend,IMAX,JMAX,it
#include "param.h"
#include "grid.h"
#include "ocean3d.h"
#include "ocean2d.h"
#include "diagnostics.h"
#include "scalars.h"
#include "forces.h"
#include "mixing.h"
      type(NPZDFe)       :: NPZ
      integer, parameter :: ITERMAX = 1,   ! number of small implicit time step
     &                      nsink   = NumVSinkTerms
      real,    parameter :: kw      = 0.04, 
     &                      kc      = 0.025, !Light attenuation of chlorophyll (mg Chl m-2)
     &                      mon2sec = 2592000.,  !How many seconds that one month contains
     &                      eps     = 1D-9
#ifdef DIAGNOSTICS_BIO
      real trend_no3,trend_phy,trend_zoo,trend_det,somme
      real bilan_no3,bilan_phy,bilan_zoo,bilan_det, sinking_loss, 
     &     trend_total
      integer l, iflux
      real ThisVSinkFlux(N, NumVSinkTerms),  ! [mmol m-2 s-1], upward flux is positive
     &     ThisFlux(N, NumFluxTerms)
#endif

      integer i,j,k, ITER, iB, itrc
      real NO3(N,2),PHY(N,2),ZOO(N,2),CHL(N,2),DFE(N,2),
     &     DET(N,2),TEM(N),FC(0:N),APMU(N,3),AVAR(N,3),
     &     par0,I_0,I_min,attn,
     &     SB(N,nsink),dSB(0:N,nsink),wSB(nsink),cff,cff1,cff2,cff6,
     &     aR,aL,cu,PMU(N,2),VAR(N,2)
# if defined DIAGNOSTICS_BIO
     &      , LastVSinkFlux,ColumnMassOld(NumVSinkTerms)
     &      , ColumnMassNew(NumVSinkTerms),dtsec
# endif /* DIAGNOSTICS_BIO */

      real, parameter :: Dust_solubility = 0.02  !Aumont et al. (2003)
     &                  ,Dust_ironfrac   = 0.035 !TOM10Appendix p. 24, unit: g/g

# include "compute_auxiliary_bounds.h"

      NPZ%env%dtdays= dt/(24.*3600.*float(ITERMAX))  ! time step as fraction of day.
      dtsec         = dt/float(ITERMAX)              ! time step in seconds
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
#ifdef DEBUG_BIO2
#  ifdef MPI
           IMAX=i+iminmpi-1
           JMAX=j+jminmpi-1
#  else
           IMAX=i
           JMAX=j
#  endif
#endif

#if defined MASKING
! Do not calculate on landmask points:
        if(rmask(i,j) .LT. 0.5) then 
           t(i,j,:,nnew,:)=0.          
#ifdef DIAGNOSTICS_BIO
           bioFlux(i,j,:,:) =0.0
           bioVSink(i,j,:,:)=0.0
#endif /* DIAGNOSTICS_BIO */
        else 

#ifdef DIAGNOSTICS_BIO
! Reset the biogeochemical fluxes. This is necessary because the
! biological routine uses multiple. time steps for each physical time
! step.
            bioFlux(i,j,:,:) =0.0             
          do k=0,N
            bioVSink(i,j,k,:)=0.0
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

      EXTRACT: do k=1,N                                  ! N: the number of vertical layers
                  do it = 5,NT
                    !Check if there are NaN values:
                    if (t(i,j,k,nnew,it) .ne. t(i,j,k,nnew,it)) then
                     write(stdout,99) it, -z_r(i,j,k)
                     write(stdout,'(A,2I6)') '  IMAX JMAX  ='
     &                  ,IMAX,JMAX
                     write(stdout,102) t(i,j,k,nnew,iNO3_),
     &                                 t(i,j,k,nnew,iPHY1),
     &                                 t(i,j,k,nnew,iZoo1)
#if defined IRON && defined IRONCYCLE
     &                                ,t(i,j,k,nnew,iFer1)
#endif

                     stop
                    endif
                  enddo
       
        NO3(k,1) =t(i,j,k,nnew,iNO3_)  ! Nitrate
        PHY(k,1) =t(i,j,k,nnew,iPhy1)  ! Phytoplankton
        ZOO(k,1) =t(i,j,k,nnew,iZoo1)  ! Zooplankton
        DET(k,1) =t(i,j,k,nnew,iDET1)  ! Detritus
#if defined IRON && defined IRONCYCLE
        DFE(k,1) =t(i,j,k,nnew,iFer1)  ! Dissolved iron
#endif
                                       ! cell volume [µm**3]*PHY BIOMASS to facilitate advection and diffusion
        PMU(k,1) =t(i,j,k,nnew,iPmu1)  ! Mean size (log um3) * PHY BIOMASS
        VAR(k,1) =t(i,j,k,nnew,iVar1)  ! Variance (log um3)^2*PHY BIOMASS
        CHL(k,1) =t(i,j,k,nnew,iChla)
        TEM(k)   =t(i,j,k,nnew,iTemp)  ! temperature; [deg. C]

      enddo EXTRACT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      BIOITER: DO ITER=1,ITERMAX      !--> Start internal iterations to achieve
                             !    nonlinear backward-implicit solution.

#if defined IRON && defined IRONCYCLE
! Iron atmospheric deposition:
! Dust deposition unit: kg/m2/month. Assume 0.035g Fe per g of dust
! so need to convert into nM at each time step
! Deposition = Dust*10^12*solubility*Fe%/56*dt/mon2sec*surface_area/surface_grid_volume/1000
        DFE(N,1) = DFE(N,1)
     &    + Dust_solubility*dust(i,j)*1D9*Dust_ironfrac/55.85
     &    / z_r(i,j,N)/mon2sec*dtsec                     
#endif

        I_0 = srflx(i,j)*rho0*Cp*0.43    ! Surface short wavelength radiation
                                                   ! and convert to the unit W m-2 
! convert back from Celsius meter second-1 to Watts meter-2
! srflx obtained from the Routine set_forces.F
! rho0: average density of seawater
! Cp  : specific heat of seawater
        I_min= 0.01*I_0

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
        par0 = I_0      ! An intermediate PAR variable.     
        do k=N,1,-1     !<-- irreversible

! Read in data into the NPZ class
          NPZ%n%Nconc    = max(NO3(k,1),eps)
#ifdef IRON
          NPZ%Fe%conc    = 
#  ifdef IRONCYCLE
     &                     max(DFE(k,1),eps)
#  endif
#endif
          NPZ%p%Nconc    = max(PHY(k,1),eps)
          NPZ%z%Nconc    = max(ZOO(k,1),eps)
          NPZ%D%Nconc    = max(DET(k,1),eps)
          CHL(k,2)       = max(CHL(k,1),0.)
          attn=exp(-0.5*(kw + kc*CHL(k,2))*     ! Attenuation in this grid
     &            (z_w(i,j,k)-z_w(i,j,k-1)))  ! z_w: depth of the grid (negative values!).

        ! Read in light and temperature:
          NPZ%env%light = par0 * attn       ! Light level at the middle of the grid
          NPZ%env%temp  = TEM(k)
! Remove meaningless PMU and VAR:
          if (NPZ%p%Nconc .LE. eps) then
             ! Simulate very oligotrophic surface waters or deep sea where there is no phytoplankton
             NPZ%p%LNV = 0.99
             NPZ%p%VAR = 0.1   !To avoid very large PMU and VAR
             ! Update new variable of PMU and VAR in cases when PHY are changed
             PMU(k,2)  = eps*NPZ%p%LNV
             VAR(k,2)  = eps*NPZ%p%VAR
          else
             ! The unit for LNV: log(10*um**3)
             NPZ%p%LNV = max(PMU(k,1),eps)/NPZ%p%Nconc
             NPZ%p%VAR = max(VAR(k,1),eps)/NPZ%p%Nconc
             PMU(k,2)  = max(PMU(k,1),eps)
             VAR(k,2)  = max(VAR(k,1),eps)
          endif
          !Save initial absolute PMU and VAR
          APMU(k,1) = NPZ%p%LNV
          AVAR(k,1) = NPZ%p%VAR

! Define I_min as 1% of I_0 (bottom of euphotic layer)
          I_min = 0.01*I_0

          call NPZ%initialize
! Main hard work:

          call NPZ%modified_patankar_2
#ifdef DIAGNOSTICS_BIO 
          ThisFlux(k, dPHY_rhs)=NPZ%p%dPHY
          ThisFlux(k, dNO3_rhs)=NPZ%n%dNO3
          ThisFlux(k, dZOO_rhs)=NPZ%z%dZOO
          ThisFlux(k, dDET_rhs)=NPZ%D%dDET
#if defined IRON && defined IRONCYCLE
          ThisFlux(k, dFER_rhs)=NPZ%Fe%dFE
#endif
          ThisFlux(k, dLNV_rhs)=NPZ%p%dLNVdt
          ThisFlux(k, dVAR_rhs)=NPZ%p%dVARdt
          ThisFlux(k, N_P)     =NPZ%pp(iPHY,iNO3)-NPZ%pp(iNO3,iPHY)
          ThisFlux(k, P_Z)     =NPZ%pp(iZOO,iPHY)
          ThisFlux(k, Z_D)     =NPZ%pp(iDET,iZOO)
          ThisFlux(k, D_N)     =NPZ%pp(iNO3,iDET)
          ThisFlux(k, Z_N)     =NPZ%pp(iNO3,iZOO)
          ThisFlux(k, dmudl)   =NPZ%p%dmuNetdl
          ThisFlux(k, dgdl)    =NPZ%z%dgdlbar
          ThisFlux(k, dwdl )   =NPZ%p%dwdl
          ThisFlux(k, d2mudl2) =NPZ%p%d2muNetdl2
          ThisFlux(k, d2gdl2)  =NPZ%z%d2gdl2bar
          ThisFlux(k, d2wdl2)  =NPZ%p%d2wdl2
#endif
          APMU(k,2) = NPZ%p%LNV
          AVAR(k,2) = NPZ%p%VAR


 102      format('NO3 = ',E15.3,' PHY = ',E20.3,' MIC =',E20.3,
     &           ' MES = ',E20.3,'  Fe = ',E20.3)


! Calculate new PAR
      par0 = NPZ%env%light*attn
! Calculate euphotic depth
      if (NPZ%env%light .ge. I_min) then
        if (par0  .ge. I_min) then
          hel(i,j)=-z_w(i,j,k-1)
        else
          hel(i,j)=-z_r(i,j,k)
        endif
      endif

      Qbulk(i,j,k)     = NPZ%p%QNavg
      Thetabulk(i,j,k) = NPZ%p%ThetaAvg

! Read in new state variables for each depth:
      PHY(k,2) = max(NPZ%p%Nconc,   eps)
#if defined IRON && defined IRONCYCLE
      DFE(k,2) = max(NPZ%Fe%conc,   eps)
#endif
      NO3(k,2) = max(NPZ%n%Nconc,   eps)
      ZOO(k,2) = max(NPZ%z%Nconc,   eps) 
      DET(k,2) = max(NPZ%D%Nconc,   eps)
! Now the NPZ%p%Nconc has been updated
      cff      = max(PHY(k,1),eps)
      PMU(k,2) =(PMU(k,2)+cff*NPZ%p%LNV + APMU(k,1)*PHY(k,2))
     &   /    (1.+ 2.*cff*APMU(k,1)/PMU(k,2))

      VAR(k,2) =(VAR(k,2)+cff*NPZ%p%VAR + AVAR(k,1)*PHY(k,2))
     &   /    (1.+ 2.*cff*AVAR(k,1)/VAR(k,2))

      enddo

! ---------------------------------------------------------------------------
! ---------------------------------------------------------------------------
! Vertical sinking: Vertical advection algorithm based on monotonic,
! continuous conservative parabolic splines.
! 1: Put biological tracers into SB array.
      do k=1,N
        SB(k,1)=PHY(k,2)
        SB(k,2)=DET(k,2)
        SB(k,3)=PMU(k,2)
        SB(k,4)=VAR(k,2)
        wSB(1) =NPZ%p%w_pAvg
        wSB(2) =NPZ%d%wDET
        wSB(3) =wSB(1)
        wSB(4) =wSB(1)
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
            FC(k)=1./Hz(i,j,k)
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
          dSB(0,iB)=0.   ! Set no-flux boundary conditions at bottom.
!
! Apply fluxes:

          do k=N,1,-1     !<-- irreversible
               FC(k)=1./Hz(i,j,k)
            SB(k,iB)=SB(k,iB)+wSB(iB)*FC(k)*(dSB(k,iB)-dSB(k-1,iB))
            SB(k,iB)=max(SB(k,iB),eps)  !Positivity
          enddo
        enddo  ! <-- iB

#ifdef DIAGNOSTICS_BIO
        do iflux = 1, NumVSinkTerms
          ColumnMassOld(iflux) = 0.0
          ColumnMassNew(iflux) = 0.0
        end do
#endif /* DIAGNOSTICS_BIO */

        do k=1,N

#ifdef DIAGNOSTICS_BIO
! ColumnMassOld and ColumnMassNew are needed to compute the sinking flux
! into the sediment
          ColumnMassOld(1)=ColumnMassOld(1)+PHY(k,2)
          ColumnMassOld(2)=ColumnMassOld(2)+PON(k,2)
          ColumnMassOld(3)=ColumnMassOld(3)+PMU(k,2)
          ColumnMassOld(4)=ColumnMassOld(4)+VAR(k,2)
 !Only calculates losses, not the whole sinking-induced flux
          ThisVSinkFlux(k, NFlux_VSinkPHY1)=PHY(k,2)-SB(k,1) 
#endif /* DIAGNOSTICS_BIO */

          PHY(k,2) = SB(k,1)
! Total Chl a (unit: mg Chl m-3)
          CHL(k,2) = PHY(k,2)/Qbulk(i,j,k)*Thetabulk(i,j,k)
#ifdef DIAGNOSTICS_BIO
          ColumnMassNew(1)=ColumnMassNew(1)+PHY(k,2)
#endif /* DIAGNOSTICS_BIO */

          ! PON
#ifdef DIAGNOSTICS_BIO
          ThisVSinkFlux(k, NFlux_VSinkPON1)=PON(k,2)-SB(k,2)
#endif /* DIAGNOSTICS_BIO */
          PON(k,2) = SB(k,2)
          PMU(k,2) = SB(k,3)
          VAR(k,2) = SB(k,4)
#ifdef DIAGNOSTICS_BIO
          ColumnMassNew(2)=ColumnMassNew(2)+PON(k,2)
          ColumnMassNew(3)=ColumnMassNew(3)+PMU(k,2)
          ColumnMassNew(4)=ColumnMassNew(4)+VAR(k,2)
#endif /* DIAGNOSTICS_BIO */
          !Calculate new phytoplankton PMU and VAR:
          APMU(k,3) = max(min(PMU(k,2)/PHY(k,2),NPZ%p%LNVmax),eps)
          AVAR(k,3) = max(min(VAR(k,2)/PHY(k,2),NPZ%p%VARmax),eps)
          PMU(k,2)  = PHY(k,2)*APMU(k,3)
          VAR(k,2)  = PHY(k,2)*AVAR(k,3)

#ifdef DIAGNOSTICS_BIO
          ThisVSinkFlux(k, NFlux_VSinkPMU1)=APMU(k,3)-APMU(k,2)
          ThisVSinkFlux(k, NFlux_VSinkVAR1)=AVAR(k,3)-APMU(k,2)
#endif /* DIAGNOSTICS_BIO */
        enddo

#ifdef DIAGNOSTICS_BIO
! Transfer fluxes to global arrays at the end of each biological time step
! for computational efficiency, divide now by dtsec to get the correct units (s-1)
        do k = 1, N
          do iflux = 1, NumFluxTerms
            bioFlux(i,j,k,iflux) = bioFlux(i,j,k,iflux) +
     &                             ThisFlux(k, iflux)/dtsec
          end do
        end do
        do iflux = 1,NumVSinkTerms 
! Compute the vertical sinking flux into the sediment by comparing
! previous and current mass in this (i,j) column
! The flux is positive if upward, so usually it will be
! negative, i.e. into the sediment.
          LastVSinkFlux = ( ColumnMassNew(iflux) -
     &                      ColumnMassOld(iflux) )
          bioVSink(i,j,0,iflux) = (  bioVSink(i,j,0,iflux) +
     &                               LastVSinkFlux /dtsec)
          do k = 1, N
            LastVSinkFlux = LastVSinkFlux + ThisVSinkFlux(k,iflux)
            bioVSink(i,j,k,iflux) =(  bioVSink(i,j,k,iflux) +
     &                                LastVSinkFlux /dtsec  )
          end do
        end do  ! <-- iflux
#endif /* DIAGNOSTICS_BIO */
!
      END DO BIOITER ! <-- ITER
!
! Write back the data into t array:
        do k=1,N
          t(i,j,k,nnew,iNO3_)=max(NO3(k,2),eps)
#if defined IRON && defined IRONCYCLE
          t(i,j,k,nnew,iFer1)=max(DFE(k,2),NPZ%Fe%Femin)
#endif
          t(i,j,k,nnew,iPhy1)=max(PHY(k,2),eps)
          t(i,j,k,nnew,iZoo1)=max(ZOO(k,2),eps)
          t(i,j,k,nnew,iDET1)=max(DET(k,2),eps)
          t(i,j,k,nnew,iChla)=max(CHL(k,2),eps)
          t(i,j,k,nnew,iPmu1)=max(PMU(k,2),eps)
          t(i,j,k,nnew,iVar1)=max(VAR(k,2),eps)
        enddo
      endif   ! <== The end of rmasking
#endif /* MASKING */
        enddo
      enddo
   98    format('Error: Abnormal values occur in',1x,a,
     &                 1x,'at Depth',1x,F6.1)

   99    format('Error: NaN values occur in',1x,a,
     &                 1x,'at Depth',1x,F6.1)

 1001    format('NO3=',1x,F6.2,1x,'PHY=',
     &                 1x,F6.2,1x,'ZOO=',
     &                 1x,F6.2,1x,'DFE=',
     &                 1x,F6.2,1x,'at Depth',1x,F6.1)

#else
      subroutine biology_empty ()
#endif
      return
      END


