#include "cppdefs.h"
      Module forces
      implicit none
#include "param.h"
!--------------------------------------------------------------------
!  SURFACE MOMENTUM FLUX (WIND STRESS):
!--------------------------------------------------------------------
!  sustr |  XI- and ETA-components of kinematic surface momentum flux
!  svstr |  (wind stresses) defined at horizontal U- and V-points.
!            dimensioned as [m^2/s^2].
!
      real sustr(GLOBAL_2D_ARRAY)
      real svstr(GLOBAL_2D_ARRAY)
!
!  tsms      Time of surface momentum stresses.
!
!  sustrg |  Two-time level gridded data for XI- and ETA-componets
!  svstrg |  of kinematic surface momentum flux (wind stess).
!
!  sustrp |  Two-time level point data for XI- and ETA-componets 
!  svstrp |  of kinematic surface momentum flux (wind stess).
!
      real sustrg(GLOBAL_2D_ARRAY,2)
      real svstrg(GLOBAL_2D_ARRAY,2)

      real    sustrp(2), svstrp(2), sms_time(2)
      real    sms_cycle, sms_scale 
      integer itsms, sms_ncycle, sms_rec, lsusgrd,    
     &        lsvsgrd,sms_tid, susid, svsid

#  undef SMFLUX_DATA
!
!  BOTTOM MOMENTUM FLUX:
!--------------------------------------------------------------------
!  bustr |  XI- and ETA-components of kinematic bottom momentum flux
!  bvstr |  (drag) defined at horizontal U- and V-points [m^2/s^2].
!
      real bustr(GLOBAL_2D_ARRAY)
      real bvstr(GLOBAL_2D_ARRAY)
#ifndef ANA_BMFLUX
!
!  tbms      Time of surface momentum stresses.
!
!  bustrg |  Two-time level gridded data for XI- and ETA-componets 
!  bvstrg |  of kinematic bottom momentum flux.
!
!  bustrp |  Two-time level point data for XI- and ETA-componets   
!  bvstrp |  of kinematic bottom momentum flux.
!
      real bustrg(GLOBAL_2D_ARRAY,2)
      real bvstrg(GLOBAL_2D_ARRAY,2)

      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
#  undef BMFLUX_DATA
#endif /* !ANA_BMFLUX */
#ifdef SOLVE3D
!
!  SURFACE TRACER FLUXES: 
!--------------------------------------------------------------------
!  stflx   Kinematic surface fluxes of tracer type variables at
!          horizontal RHO-points. Physical dimensions [degC m/s] -
!          temperature; [PSU m/s] - salinity.
!
      real stflx(GLOBAL_2D_ARRAY,NT)
# ifdef BULK_FLUX
      real shflx_rsw(GLOBAL_2D_ARRAY)
      real shflx_rlw(GLOBAL_2D_ARRAY)
      real shflx_lat(GLOBAL_2D_ARRAY)
      real shflx_sen(GLOBAL_2D_ARRAY)
# endif
# ifdef SST_SKIN
      real sst_skin(GLOBAL_2D_ARRAY)
      real dT_skin(GLOBAL_2D_ARRAY)
# endif
# if !defined ANA_STFLUX || !defined ANA_SSFLUX
!
!  stflxg   Two-time level surface tracer flux grided data.
!  stflxp   Two-time level surface tracer flux point  data.
!  tstflx   Time of surface tracer flux.
!
      real stflxg(GLOBAL_2D_ARRAY,2,NT)
      real stflxp(2,NT), stf_time(2,NT)
      real stf_cycle(NT), stf_scale(NT)
      integer itstf(NT), stf_ncycle(NT), stf_rec(NT),   
     &        lstfgrd(NT), stf_tid(NT), stf_id(NT)
#   undef STFLUX_DATA
# endif /* !ANA_STFLUX || !ANA_SSFLUX */
!
!  BOTTOM TRACER FLUXES:
!--------------------------------------------------------------------
!  btflx  Kinematic bottom fluxes of tracer type variables at
!         horizontal RHO-points. Physical dimensions [degC m/s] -
!         temperature; [PSU m/s] - salinity.
!
      real btflx(GLOBAL_2D_ARRAY,NT)
# ifndef ANA_BTFLUX
!
!  btflxg   Two-time level bottom tracer flux grided data.
!  btflxp   Two-time level bottom tracer flux point data.
!  tbtflx   Time of bottom tracer flux.
!
      real btflxg(GLOBAL_2D_ARRAY,2,NT)
      real sclbtf(NT), btf_tstart(NT), btf_tend(NT), 
     &     btfclen(NT), tsbtf(NT)
      real btf_tintrp(2,NT), btflxp(2,NT),  tbtflx(2,NT)
      integer itbtf(NT), btfid(NT), btftid(NT),tbtfindx(NT)
      logical lbtfgrd(NT), btfcycle(NT), btf_onerec(NT)
#   undef BTFLUX_DATA
# endif /* !ANA_BTFLUX */
#ifdef QCORRECTION
!  HEAT FLUX CORRECTION
!--------------------------------------------------------------------
!  dqdt     Kinematic surface net heat flux sensitivity to SST [m/s].
!  sst      Current sea surface temperature [degree Celsius].
!
      real dqdt(GLOBAL_2D_ARRAY)
      real sst(GLOBAL_2D_ARRAY)
#  ifndef ANA_SST
!
!  dqdtg |  Two-time-level grided data for net surface heat flux
!  sstg  |  sensitivity to SST grided data [Watts/m^2/Celsius] and
!              sea surface temperature [degree Celsius].
!  dqdtp |  Two-time-level point data for net surface heat flux
!  sstp  |  sensitivity to SST grided data [Watts/m^2/Celsius] and
!              sea surface temperature [degree Celsius].
!  tsst     Time of sea surface temperature data.
!
      real dqdtg(GLOBAL_2D_ARRAY,2)
      real sstg(GLOBAL_2D_ARRAY,2)
      real    sstp(2), dqdtp(2), sst_time(2)
      real    sst_cycle, scldqdt
      integer itsst, sst_ncycle, sst_rec,  sst_tid,  sst_id,
     &        dqdt_id,     lsstgrd,   sstunused
#    undef SST_DATA
#  endif /* !ANA_SST */
# endif /* QCORRECTION */
# if defined SALINITY && defined SFLX_CORR
!
!  SALT FLUX CORRECTION
!--------------------------------------------------------------------
!  sss      Current sea surface salinity [PSU].
!
      real sss(GLOBAL_2D_ARRAY)
#  if !defined QCORRECTION
      real dqdt(GLOBAL_2D_ARRAY)
#  endif
#  ifndef ANA_SSS
!
!  dqdtg |  Two-time-level grided data for net surface heat flux
!  sssg  |  Two-time-level grided data for 
!              sea surface salinity [PSU].
!  dqdtp |  Two-time-level point data for net surface heat flux
!  sssp  |  Two-time-level point data for
!              sea surface salinity [PSU].
!  tsss     Time of sea surface salinity data.
!
      real sssg(GLOBAL_2D_ARRAY,2)
      real sssp(2),  sss_time(2)
      real sss_cycle
      integer itsss, sss_ncycle, sss_rec,  sss_tid,  sss_id,
     &                          lsssgrd,   sssunused
#   if !defined QCORRECTION
      real dqdtg(GLOBAL_2D_ARRAY,2)
      real    dqdtp(2)
      real    scldqdt
      integer dqdt_id
#   endif
#   undef SSS_DATA
#  endif /* !ANA_SSS */
# endif /* SALINITY && SFLX_CORR */
!
!
#ifdef BULK_FLUX
!
!  HEAT FLUX BULK FORMULATION
!--------------------------------------------------------------------
!  tair     surface air temperature at 2m [degree Celsius].
!  wsp      wind speed at 10m [degree Celsius].
!  rhum     surface air relative humidity 2m [fraction]
!  prate    surface precipitation rate [cm day-1]
!  radlw    net terrestrial longwave radiation [Watts meter-2]
!  radsw    net solar shortwave radiation [Watts meter-2]
!
      real tair(GLOBAL_2D_ARRAY)
      real rhum(GLOBAL_2D_ARRAY)
      real prate(GLOBAL_2D_ARRAY)
      real radlw(GLOBAL_2D_ARRAY)
      real radsw(GLOBAL_2D_ARRAY)
      real wspd(GLOBAL_2D_ARRAY)
# ifdef BULK_SM_UPDATE
      real uwnd(GLOBAL_2D_ARRAY)
      real vwnd(GLOBAL_2D_ARRAY)
# endif
# ifdef DIURNAL_INPUT_SRFLX
      real radswbio(GLOBAL_2D_ARRAY)
# endif

      real tairg(GLOBAL_2D_ARRAY,2)
      real rhumg(GLOBAL_2D_ARRAY,2)
      real prateg(GLOBAL_2D_ARRAY,2)
      real radlwg(GLOBAL_2D_ARRAY,2)
      real radswg(GLOBAL_2D_ARRAY,2)
      real wspdg(GLOBAL_2D_ARRAY,2)
# ifdef BULK_SM_UPDATE
      real uwndg(GLOBAL_2D_ARRAY,2)
      real vwndg(GLOBAL_2D_ARRAY,2)
# endif
# ifdef DIURNAL_INPUT_SRFLX
      real radswbiog(GLOBAL_2D_ARRAY,2)
# endif

      real    tairp(2),rhump(2),pratep(2),radlwp(2),radswp(2)
     &       ,wspdp(2)
# ifdef BULK_SM_UPDATE
     &       ,uwndp(2),vwndp(2)
# endif
# ifdef DIURNAL_INPUT_SRFLX
     &       ,radswbiop(2)
# endif
!      real    bulk_time(2), bulk_cycle
      real     tair_time(2),  tair_cycle
      real     rhum_time(2),  rhum_cycle
      real    prate_time(2), prate_cycle
      real      rad_time(2),   rad_cycle
      real     wspd_time(2),  wspd_cycle
      real     uwnd_time(2),  uwnd_cycle
      real     vwnd_time(2),  vwnd_cycle

      integer tair_id,rhum_id,prate_id,radlw_id,radsw_id,
     &        ltairgrd,lrhumgrd,lprategrd,lradlwgrd,lradswgrd
     &       ,wspd_id,lwspdgrd
# ifdef BULK_SM_UPDATE
     &       ,uwnd_id,vwnd_id,luwndgrd,lvwndgrd
# endif
# ifdef DIURNAL_INPUT_SRFLX
     &        ,radswbio_id,lradswbiogrd
# endif
      integer ittair,tair_ncycle,tair_rec,tair_tid,tairunused
      integer itrhum,rhum_ncycle,rhum_rec,rhum_tid,rhumunused
      integer itprate,prate_ncycle,prate_rec,prate_tid,prateunused
      integer itrad,rad_ncycle,rad_rec,rad_tid,radunused
      integer itwspd,wspd_ncycle,wspd_rec,wspd_tid,wspdunused
      integer ituwnd,uwnd_ncycle,uwnd_rec,uwnd_tid,uwndunused
      integer itvwnd,vwnd_ncycle,vwnd_rec,vwnd_tid,vwndunused
!      integer itbulk,bulk_ncycle,bulk_rec,bulk_tid,
!     &        bulkunused

#endif /* BULK_FLUX */
!
!  SOLAR SHORT WAVE RADIATION FLUX.
!--------------------------------------------------------------------
!  srflx  Kinematic surface shortwave solar radiation flux
!         [degC m/s] at horizontal RHO-points
!
      real :: srflx(GLOBAL_2D_ARRAY)
# ifdef ANA_DIURNAL_SW
      real sin_phi(GLOBAL_2D_ARRAY),  
     &     cos_phi(GLOBAL_2D_ARRAY), 
     &     tan_phi(GLOBAL_2D_ARRAY)
# endif
# ifdef DIURNAL_INPUT_SRFLX
      real srflxbio(GLOBAL_2D_ARRAY)
# endif
# ifndef ANA_SRFLUX
!
!  srflxg | Two-time-level grided and point data for surface 
!  srflxp |      solar shortwave radiation flux grided data.
!  tsrflx   Time of solar shortwave radiation flux.
!
      real srflxg(GLOBAL_2D_ARRAY,2)

      real srflxp(2),srf_time(2)
      real srf_cycle, srf_scale
      integer itsrf, srf_ncycle, srf_rec, 
     &        lsrfgrd, srf_tid, srf_id

# ifdef DIURNAL_INPUT_SRFLX
      real srflxbiog(GLOBAL_2D_ARRAY,2)

      real srflxbiop(2)
      integer srfbio_rec, lsrfbiogrd, 
     &        srfbio_tid, srfbio_id
# endif /*  DIURNAL_INPUT_SRFLX   */

#   undef SRFLUX_DATA
# endif /* !ANA_SRFLUX */
#endif /* SOLVE3D */

#if defined BIOLOGY && defined PZEFTFE
!
!  Atmospheric Dust deposition FLUX.
!--------------------------------------------------------------------
!  dust  Flux of atmospheric dust deposition
!         [kg/m2/month] at horizontal RHO-points
!
      real :: dust(GLOBAL_2D_ARRAY)

!
!  dustg | Two-time-level grided and point data for surface 
!  dustp |      dust flux grided data.
!  tdust   Time of dust flux.
!
      real dustg(GLOBAL_2D_ARRAY,2)

      real dustp(2),  dust_time(2)
      real dust_cycle, dust_scale
      integer itdust, dust_ncycle, dust_rec, 
     &        ldustgrd, dust_tid, dust_id
#endif
!--------------------------------------------------------------------
!  WIND INDUCED WAVES: everything is defined at rho-point
!--------------------------------------------------------------------
! wfrq | BBL/MRL | wind-induced wave frequency [rad/s]
! uorb | BBL     | xi-component  of wave-induced bed orbital velocity [m/s]
! vorb | BBL     | eta-component of wave-induced bed orbital velocity [m/s]
! wdrx | MRL     | cosine of wave direction [non dimension]
! wdre | MRL     | sine of   wave direction [non dimension]
! whrm | MRL     | (RMS) wave height (twice the wave amplitude) [m]
! wdsp | MRL     | breaking dissipation rate (\epsilon_b term) [m3/s3]
! wdrg | MRL     | frictional dissipation rate (\epsilon_d term) [m3/s3]
! rdsp | ROLLER  | roller dissipation rate (\epsilon_r term) [m3/s3]
! wbst | MRL/BKPP| frictional dissipation stress (e_d k/sigma) [m2/s2]
!--------------------------------------------------------------------

#if defined BBL || defined MRL_WCI 
      real wfrq(GLOBAL_2D_ARRAY)
#endif

#ifdef BBL
      real uorb(GLOBAL_2D_ARRAY)
      real vorb(GLOBAL_2D_ARRAY)
#endif   /* BBL */

#if defined MRL_WCI 
      real whrm(GLOBAL_2D_ARRAY)
      real wdsp(GLOBAL_2D_ARRAY)
      real wdrg(GLOBAL_2D_ARRAY)
      real wbst(GLOBAL_2D_ARRAY)
      real wdrx(GLOBAL_2D_ARRAY)
      real wdre(GLOBAL_2D_ARRAY)
# ifdef WAVE_ROLLER
      real rdsp(GLOBAL_2D_ARRAY)
# endif
!
!--------------------------------------------------------------------
!  WAVE AVEREAGED QUANTITIES AND TERMS
!--------------------------------------------------------------------
!  2D  |  brk2dx   |   xi-direciton 2D breaking dissipation (rho)
!  2D  |  brk2de   |  eta-direction 2D breaking dissipation (rho)
!  2D  |  frc2dx   |   xi-direciton 2D frictional dissipation (rho)
!  2D  |  frc2de   |  eta-direction 2D frictional dissipation (rho)
!  2D  |  ust2d    |   xi-direciton Stokes transport (u-point)
!  2D  |  vst2d    |  eta-direciton Stokes transport (v-point)
!  2D  |  sup      |  quasi-static wave set-up (rho-point)
!  2D  |  calP     |  pressure correction term (rho-point)
!  2D  |  Kapsrf   |  Bernoulli head terrm at the surface (rho-point)
!--------------------------------------------------------------------
!  3D  |  brk3dx   |   xi-direciton 3D breaking dissipation (rho)
!  3D  |  brk3de   |  eta-direction 3D breaking dissipation (rho)
!  3D  |  ust      |   xi-direciton 3D Stokes drift velocity (u-point)
!  3D  |  vst      |  eta-direciton 3D Stokes drift velocity (v-point)
!  3D  |  wst      |       vertical 3D Stokes drift velocity (rho-point)
!  3D  |  Kappa    |  3D Bernoulli head term (rho-point)
!  3D  |  kvf      |  vertical vortex force term (K term, 3D, rho-point)
!  3D  |  Akb      |  breaking-wave-induced additional diffusivity (w-point)
!  3D  |  Akw      |  wave-induced additional diffusivity (rho-point)
!  3D  |  E_pre    |  previous time-step value for Akw estimation (rho)
!  3D  |  frc3dx   |   xi-direciton 3D frictional dissipation (rho)
!  3D  |  frc3de   |  eta-direction 3D frictional dissipation (rho)
!--------------------------------------------------------------------
!
      real brk2dx(GLOBAL_2D_ARRAY)
      real brk2de(GLOBAL_2D_ARRAY)
      real ust2d(GLOBAL_2D_ARRAY)
      real vst2d(GLOBAL_2D_ARRAY)
      real frc2dx(GLOBAL_2D_ARRAY)
      real frc2de(GLOBAL_2D_ARRAY)
      real sup(GLOBAL_2D_ARRAY)
# ifdef SOLVE3D
      real calP(GLOBAL_2D_ARRAY)
      real Kapsrf(GLOBAL_2D_ARRAY)
#  ifndef WAVE_SFC_BREAK
      real brk3dx(GLOBAL_2D_ARRAY,N)
      real brk3de(GLOBAL_2D_ARRAY,N)
#  endif
#  ifdef WAVE_BODY_STREAMING
      real frc3dx(GLOBAL_2D_ARRAY,N)
      real frc3de(GLOBAL_2D_ARRAY,N)
#  endif
      real ust(GLOBAL_2D_ARRAY,N)
      real vst(GLOBAL_2D_ARRAY,N)
      real wst(GLOBAL_2D_ARRAY,N)
      real kvf(GLOBAL_2D_ARRAY,N)
      real Akb(GLOBAL_2D_ARRAY,0:N)
      real Akw(GLOBAL_2D_ARRAY,0:N)
      real E_pre(GLOBAL_2D_ARRAY,0:N)
# endif  /* SOLVE3D */
#endif   /* MRL_WCI */

#if defined BBL || defined MRL_WCI 
!--------------------------------------------------------------------
! Awave  | for present time   | wave amplitude [m]
! Pwave  | for present time   | wave direction [radians]
! Dwave  | for present time   | wave period [s]
!--------------------------------------------------------------------
      real Awave(GLOBAL_2D_ARRAY)
      real Dwave(GLOBAL_2D_ARRAY)
      real Pwave(GLOBAL_2D_ARRAY)

# if !defined ANA_WWAVE && !defined WKB_WWAVE
!--------------------------------------------------------------------
!  Eb    |                    | breaking dissipation [m3/s3]
!  wved  |  for present time  | frictional dissipation [m3/s3]
!  wvqb  |  step              | fraction of breaking waves [ND]
!--------------------------------------------------------------------
!wwv_time|                    | time of wind-induced waves
!--------------------------------------------------------------------
!  wwag  |                    | wave amplitude [m]
!  wwdg  |                    | wave direction [radians]
!  wwpg  |                    | wave period [s]
!  wwub  |                    | orbital velocity magnitude [m/s]
!  wwfrq |  Two-time-level    | wave frequency [rad/s]
!  wwuob |  point data        | xi-orbital velocity [m/s]
!  wwvob |  for wind induced  ! eta-orbital velocity [m/s]
!  wwdrx |                    | cosine wave direction [ND]
!  wwdre |                    ! sine wave direction [ND]
!  wwhrm |                    ! (RMS) wave height [m]
!  wweb  |                    ! breaking dissipation [m3/s3]
!  wwed  |                    ! frictional dissipation [m3/s3]
!  wwqb  |                    ! fraction of breaking waves [ND]
!--------------------------------------------------------------------
!
      real wwag(GLOBAL_2D_ARRAY,2)
      real wwdg(GLOBAL_2D_ARRAY,2)
      real wwpg(GLOBAL_2D_ARRAY,2)
      real wwfrq(GLOBAL_2D_ARRAY)
#  ifdef BBL_OFFLINE 
#   ifdef WAVE_OFFLINE
      real wwub(GLOBAL_2D_ARRAY,2)
#   endif /* WAVE_OFFLINE */
      real wwuob(GLOBAL_2D_ARRAY,2)
      real wwvob(GLOBAL_2D_ARRAY,2)
#  endif /* BBL */
#  ifdef MRL_WCI 
      real wwhrm(GLOBAL_2D_ARRAY)
      real wwdrx(GLOBAL_2D_ARRAY)
      real wwdre(GLOBAL_2D_ARRAY)
      real wweb(GLOBAL_2D_ARRAY,2)
      real Eb  (GLOBAL_2D_ARRAY)
#   ifdef WAVE_OFFLINE
      real wved(GLOBAL_2D_ARRAY)
      real wwed(GLOBAL_2D_ARRAY,2)
#   endif /* WAVE_OFFLINE */
#   if defined WAVE_OFFLINE && defined WAVE_ROLLER
      real wvqb(GLOBAL_2D_ARRAY)
      real wwqb(GLOBAL_2D_ARRAY,2)
#   endif
#  endif 
      real ww_cycle,wwv_time(2),wwap(2), wwdp(2),wwpp(2),wwep(2),
     &                 wwa_scale, wwd_scale, wwp_scale,wwe_scale,  
     &                 wwagrd,   wwdgrd,    wwpgrd, wwegrd
      integer ww_ncycle,  ww_rec,  itww,
     &        ww_file_id, ww_tid,  wwa_id, wwp_id, wwd_id
#   if defined BBL && defined WAVE_OFFLINE
     &       ,wwu_id
#   endif
#   if defined MRL_WCI && defined WAVE_OFFLINE
     &       ,wwe_id, wwq_id, wwf_id
#   endif
#   if defined MRL_WCI && defined WAVE_OFFLINE
     &       ,wwe_id, wwq_id, wwf_id
#   endif
# elif defined MRL_WCI
      real Eb(GLOBAL_2D_ARRAY)
      real wved(GLOBAL_2D_ARRAY)
# endif /* ANA_WWAVE && !WKB_WWAVE */
#endif /* BBL || MRL_WCI */
      End Module

