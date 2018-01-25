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
      USE PZEFT
!
! Compute biological forcing functions
!
! In this particular implementation there are 6 state variables:
! NO3, PHYTOplankton, ZOOplankton, DETritus, Mean size, and variance.
!
      implicit none
      integer Istr,Iend,Jstr,Jend,IMAX,JMAX
#include "param.h"
#include "grid.h"
#include "ocean3d.h"
#include "ocean2d.h"
#include "diagnostics.h"
#include "scalars.h"
#include "forces.h"
#include "mixing.h"
      type(type_PZEFT)   :: NPZ
      integer, parameter :: ITERMAX = 1,   ! number of small implicit time step
     &                      nsink   = NumVSinkTerms
      real,    parameter :: kw      = 0.04, 
     &                      kc      = 0.025 !Light attenuation of chlorophyll (mg Chl m-2)

#ifdef DIAGNOSTICS_BIO
      real trend_no3,trend_phy,trend_zoo,trend_det,somme
      real bilan_no3,bilan_phy,bilan_zoo,bilan_det, sinking_loss, 
     &     trend_total
      integer l, iflux
      real ThisVSinkFlux(N, NumVSinkTerms),  ! [mmol m-2 s-1], upward flux is positive
     &     ThisFlux(N, NumFluxTerms)
#endif

      integer i,j,k, ITER, iB, itrc
      real NO3(N),PHY(N),ZOO(N),DET(N),CHL(N),
     &     TEM(N),FC(0:N),
     &     par0,I_min,attn,
     &     SB(N,nsink),dSB(0:N,nsink),wSB(nsink),cff,cff1,cff2,cff6,
     &     aR,aL,cu
# if defined NPZDSIZE || defined EFTSIZE
     &      ,PMU(N),VAR(N)
# endif
# if defined DIAGNOSTICS_BIO
     &      , LastVSinkFlux,ColumnMassOld(NumVSinkTerms)
     &      , ColumnMassNew(NumVSinkTerms),dtsec
# endif /* DIAGNOSTICS_BIO */

# include "compute_auxiliary_bounds.h"

      NPZ%env(idtdays)= dt/(24.*3600.*float(ITERMAX))  ! time step as fraction of day.
# if defined DIAGNOSTICS_BIO || defined OXYGEN
      dtsec = dt / float(ITERMAX)           ! time step in seconds
# endif /* DIAGNOSTICS_BIO || OXYGEN */
      

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
        
        NO3(k) =t(i,j,k,nnew,iNO3_)   ! Nitrate
        PHY(k) =t(i,j,k,nnew,iPhy1)   ! Phytoplankton
        ZOO(k) =t(i,j,k,nnew,iZoo1)   ! Zooplankton
        DET(k) =t(i,j,k,nnew,iDet1)   ! Detritus
 

# if defined NPZDSIZE || defined EFTSIZE
! cell volume [µm**3]*PHY BIOMASS to facilitate advection and diffusion
        PMU(k) =t(i,j,k,nnew,iPmu1)   ! Mean size (log um3) * PHY BIOMASS
        VAR(k) =t(i,j,k,nnew,iVar1)   ! Variance (log um3)^2*PHY BIOMASS
# endif
        CHL(k) =t(i,j,k,nnew,iChla)
        TEM(k) =t(i,j,k,nnew,iTemp)   ! temperature; [deg. C]

      enddo EXTRACT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      BIOITER: DO ITER=1,ITERMAX      !--> Start internal iterations to achieve
                             !    nonlinear backward-implicit solution.

        NPZ%env(iI_0) = srflx(i,j)*rho0*Cp*0.43    ! Surface short wavelength radiation
                                                   ! and convert to the unit W m-2 
! convert back from Celsius meter second-1 to Watts meter-2
! srflx obtained from the Routine set_forces.F
! rho0: average density of seawater
! Cp  : specific heat of seawater
        I_min   = 0.01*NPZ%env(iI_0)

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
        par0 = NPZ%env(iI_0) ! An intermediate PAR variable.     
        do k=N,1,-1     !<-- irreversible

! Read in data into the NPZ class
          NPZ%vars(iNO3) = NO3(k)
          NPZ%vars(iPHY) = PHY(k)
          NPZ%vars(iZOO) = ZOO(k)
          NPZ%vars(iDET) = DET(k)
# if defined NPZDSIZE || defined EFTSIZE
          NPZ%vars(iPMU) = PMU(k)
          NPZ%vars(iVAR) = VAR(k)
# endif


          attn=exp(-0.5*(kw + kc*CHL(k))*     ! Attenuation in this grid
     &            (z_w(i,j,k)-z_w(i,j,k-1)))  ! z_w: depth of the grid (negative values!).

        ! Read in light and temperature:
          NPZ%env(ipar) = par0 * attn       ! Light level at the middle of the grid
          NPZ%env(itC ) = TEM(k)

! Define I_min as 1% of I_0 (bottom of euphotic layer)
      I_min = 0.01*NPZ%env(iI_0)

           call NPZ%initialize
! Main hard work:

           call NPZ%modified_patankar_2

         ! Guard against abnormal values
!        if (NPZ%vars(iPHY).gt.1000.) then
!         write(stdout,98) VarLabel(iPhy1),-z_r(i,j,k)
!         write(stdout,'(A,2I6)') '  IMAX JMAX  ='
!     &      ,IMAX,JMAX
!         write(stdout,102) NPZ%vars(iNO3),NPZ%vars(iPHY),NPZ%vars(iZOO)
!         write(stdout,103) NPZ%vars(iDET),NPZ%vars(iPMU),NPZ%vars(iVAR)
!         stop
!        endif
! 102    format('NO3 = ',E15.3,' PHY = ',E20.3,' ZOO = ',E20.3)
! 103    format('DET = ',E15.3,' PMU = ',E20.3,' VAR = ',E20.3)

#ifdef DEBUG_BIO2
!        if (IMAX.eq.30.and. JMAX.eq.19.and.mod(k,10).eq.0) then
!           write(stdout,102) NPZ%vars(iNO3),NPZ%vars(iPHY),NPZ%vars(iZOO)
!           write(stdout,103) NPZ%vars(iDET),NPZ%vars(iPMU),NPZ%vars(iVAR)
!           write(stdout,104) NPZ%pp(2,1),NPZ%pp(1,2),NPZ%pp(3,2)
!        endif


!        if (NPZ%vars(iNO3).gt.1000.) then
!         write(stdout,98) VarLabel(iNO3_),-z_r(i,j,k)
!         write(stdout,'(A,2I6)') '  IMAX JMAX  ='
!     &      ,IMAX,JMAX
!         write(stdout,102) NPZ%vars(iNO3),NPZ%vars(iPHY),NPZ%vars(iZOO)
!         write(stdout,103) NPZ%vars(iDET),NPZ%vars(iPMU),NPZ%vars(iVAR)
!         stop
!        endif

 104    format('PP_pn =',E15.3,' PP_np =',E15.3,' PP_zp= ',E15.3)
#endif
         

#ifdef DIAGNOSTICS_BIO
      !  NUT ==> PHY
      ThisFlux(k, NFlux_Cpp)    =NPZ%pp(2,1)-NPZ%pp(1,2)          !Unit: mMol N m-3 d-1
      !  PHY ==> ZOO
      ThisFlux(k, NFlux_Grazing)=NPZ%pp(3,2)
      !  PHY ==> DET
      ThisFlux(k, NFlux_Pmort)  =NPZ%pp(4,2)
      !  ZOO ==> DET1
      ThisFlux(k, NFlux_Zmort)  =NPZ%vars(iZOO)**2
     &       *                   NPZ%pars(imz)*NPZ%tf()
      ! ZOO ==> DET2
      ThisFlux(k, NFlux_Zeges)  =NPZ%vars(iZOO)*NPZ%EGES
      ! ZOO ==> NUT
      ThisFlux(k, NFlux_Zremin) =NPZ%pp(1,3)
      ! DET ==> NUT
      ThisFlux(k, NFlux_ReminD) =NPZ%pp(1,4)
#  ifdef EFTSIZE
      ThisFlux(k, Nflux_dPMUdt) =NPZ%dPMUdt
      ThisFlux(k, Nflux_dVARdt) =NPZ%dVARdt
      ThisFlux(k, Nflux_dmudl)  =NPZ%dmuNetdl
      ThisFlux(k, Nflux_dgdl )  =NPZ%dgdlbar
      ThisFlux(k, Nflux_dwdl )  =NPZ%dwdl
      ThisFlux(k, Nflux_d2mudl2)=NPZ%d2muNetdl2
      ThisFlux(k, Nflux_d2gdl2) =NPZ%d2gdl2bar
      ThisFlux(k, Nflux_d2wdl2) =NPZ%d2wdl2
#  endif
#endif /* DIAGNOSTICS_BIO */
! Calculate new PAR
      par0 = NPZ%env(ipar)*attn
! Calculate euphotic depth
      if (NPZ%env(ipar) .ge. I_min) then
        if (par0  .ge. I_min) then
          hel(i,j)=-z_w(i,j,k-1)
        else
          hel(i,j)=-z_r(i,j,k)
        endif
      endif

      Qbulk(i,j,k)     = NPZ%QNavg
      Thetabulk(i,j,k) = NPZ%ThetaAvg

!#  ifdef VERBOSE
!      if (k .eq. 35 .and. i .eq. 100 .and. j .eq. 70) then
!        write(stdout,*) 'mu =',
!     &    NPZ%muNet/NPZ%env(idtdays),'PAR= ',NPZ%env(ipar)
!      endif
!#  endif /* VERBOSE */ 
!#ifdef AVERAGES
! if the sun is down, set the logical variable "sun_down"
! to true for not taking into account this time step in the averaging
!          if ((ZEROTH_TILE).and.(srflx(Istr,Jstr).eq.0.)) then
!            sun_down=.true.
!          endif
!#endif
!          hel(i,j)=0.0
!        endif
!

! Read in new state variables for each depth:
         PHY(k) = max(NPZ%vars(iPHY),0.)
         NO3(k) = max(NPZ%vars(iNO3),0.)
         ZOO(k) = max(NPZ%vars(iZOO),0.) 
         DET(k) = max(NPZ%vars(iDET),0.)
# if defined NPZDSIZE || defined EFTSIZE
         PMU(k) = max(NPZ%vars(iPMU),0.)
         VAR(k) = max(NPZ%vars(iVAR),0.)
# endif

      enddo

! ---------------------------------------------------------------------------
! ---------------------------------------------------------------------------
! Vertical sinking: Vertical advection algorithm based on monotonic,
! continuous conservative parabolic splines.
! 1: Put biological tracers into SB array.
      do k=1,N
        SB(k,1)=PHY(k)
        SB(k,2)=DET(k)
        SB(k,3)=PMU(k)
        SB(k,4)=VAR(k)
        wSB(1) =NPZ%w_pAvg
        wSB(2) =NPZ%pars(iw_d)
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
!
!          do k=1,N-1
!            dSB(k,iB)=SB(k+1,iB)
!          enddo

          do k=N,1,-1     !<-- irreversible
            FC(k)=1./Hz(i,j,k)
            SB(k,iB)=SB(k,iB)+wSB(iB)*FC(k)*(dSB(k,iB)-dSB(k-1,iB))
!            cu=wSB(iB)*FC(k)*dSB(k,iB)  !flux from top
!            SB(k,iB)=SB(k,iB)+cu
!
!# if defined NPZDSIZE || defined EFTSIZE
!            if (iB.eq.1) then
!            !update PMU and VAR
!               PMU(k)=cu*   PMU(k+1)+SB(k,iB)         *PMU(k)
!               VAR(k)=cu*cu*VAR(k+1)+SB(k,iB)*SB(k,iB)*VAR(k)
!            endif
!# endif
!            SB(k,iB)=SB(k,iB)-wSB(iB)*FC(k)*dSB(k-1,iB)
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
          ColumnMassOld(1)=ColumnMassOld(1)+PHY(k)
 !Only calculates losses, not the whole sinking-induced flux
          ThisVSinkFlux(k, NFlux_VSinkP1)=PHY(k)-SB(k,1) 
#endif /* DIAGNOSTICS_BIO */

          PHY(k) = SB(k,1)
! Total Chl a (unit: mg Chl m-3)
          CHL(k) = PHY(k)/Qbulk(i,j,k)*Thetabulk(i,j,k)

#ifdef DIAGNOSTICS_BIO
          ColumnMassNew(1)=ColumnMassNew(1)+PHY(k)
#endif /* DIAGNOSTICS_BIO */

          ! detritus
#ifdef DIAGNOSTICS_BIO
          ColumnMassOld(2)=ColumnMassOld(2) +DET(k)
          ThisVSinkFlux(k, NFlux_VSinkD1) =DET(k)-SB(k,2)
#endif /* DIAGNOSTICS_BIO */

          DET(k) = SB(k,2)
          PMU(k) = SB(k,3)
          VAR(k) = SB(k,4)
#ifdef DIAGNOSTICS_BIO
          ColumnMassNew(2)=ColumnMassNew(2) +DET(k)
#endif /* DIAGNOSTICS_BIO */
        enddo

#ifdef DIAGNOSTICS_BIO
! Transfer fluxes to global arrays at the end of each biological time step
! for computational efficiency, divide now by dtdays to get the correct units (d-1)
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
        end do
#endif /* DIAGNOSTICS_BIO */
!
      END DO BIOITER ! <-- ITER
!
! Write back
# if defined DIAGNOSTICS_BIO && defined DEBUG_BIO1
!
      k=N
      if ((i.eq.13).and.(j.eq.15)) then
        bilan_no3 = bioFlux(i,j,k,NFlux_ReminD)
     &            + bioFlux(i,j,k,NFlux_Zremin)
     &            - bioFlux(i,j,k, NFlux_Cpp)
!
        bilan_phy = bioFlux(i,j,k,NFlux_Cpp)
     &           - bioFlux(i,j,k,NFlux_Pmort)
     &           - bioFlux(i,j,k,NFlux_Grazing)
     &           - ( bioVSink(i,j,k,NFlux_VSinkP1) 
     &               - bioVSink(i,j,k-1,NFlux_VSinkP1) )
!
        bilan_zoo = bioFlux(i,j,k,NFlux_Grazing) 
     &            - bioFlux(i,j,k,NFlux_Zeges)
     &            - bioFlux(i,j,k,NFlux_Zmort)
     &            - bioFlux(i,j,k,NFlux_Zremin)
!
        bilan_det = bioFlux(i,j,k,NFlux_Zmort)
     &            + bioFlux(i,j,k,NFlux_Zeges)
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
        ! Changes in one day:
        trend_no3 = ( (min(t(i,j,k,nnew,iNO3_),0.) +NO3(k))
     &               - t(i,j,k,nnew,iNO3_) ) / NPZ%env(idtdays)
        trend_phy = ( (min(t(i,j,k,nnew,iPhy1),0.) +Phyt(k))
     &               - t(i,j,k,nnew,iPhy1) )/ NPZ%env(idtdays)
        trend_zoo = ( (min(t(i,j,k,nnew,izoo1),0.) +Zoo(k))
     &               - t(i,j,k,nnew,izoo1) )/ NPZ%env(idtdays)
        trend_det = ( (min(t(i,j,k,nnew,idet1),0.) +Det(k))
     &       - t(i,j,k,nnew,idet1) )/ NPZ%env(idtdays)
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
        print*, '=================='
      endif
# endif


          do k=1,N
            t(i,j,k,nnew,iNO3_)=max(NO3(k),0.0)
            t(i,j,k,nnew,iPhy1)=max(PHY(k),0.0)
            t(i,j,k,nnew,iZoo1)=max(ZOO(k),0.0)
            t(i,j,k,nnew,iDet1)=max(DET(k),0.0)
            t(i,j,k,nnew,iChla)=max(CHL(k),0.0)
# if defined NPZDSIZE || defined EFTSIZE
            t(i,j,k,nnew,iPmu1)=max(PMU(k),0.0)
            t(i,j,k,nnew,iVar1)=max(VAR(k),0.0)
# endif
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
     &                 1x,F6.2,1x,'DET=',
     &                 1x,F6.2,1x,'at Depth',1x,F6.1)

#else
      subroutine biology_empty ()
#endif
      return
      END


