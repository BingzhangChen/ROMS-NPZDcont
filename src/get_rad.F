#include "cppdefs.h"
#if defined BULK_FLUX

      subroutine get_rad
      USE forces, only: itrad,rad_ncycle,rad_rec,rad_tid,radlw_id,
     &    radsw_id, lradlwgrd, lradswgrd,rad_time,
     &    rad_cycle,radlwg,radswg, radlwp, radswp, srf_scale
!
!  Read in wind speed and surface air temperature
!
!

      implicit none
# include "param.h"
# include "scalars.h"
# include "netcdf.inc"
# include "ncscrum.h"
      real    cff
      integer i,ierr, lstr,lvar,lenstr, nf_fread, advance_cycle
       
      radname   = trim(nameit)//'_rad.nc'
!
! Initialization: Inquire about the contents of forcing NetCDF file:
!================  variables and dimensions. Check for consistency.
!
      if (may_day_flag.ne.0) return      !-->  EXIT

      if (itrad.eq.0 .or. iic.eq.0) then
        lstr=lenstr(radname)
c*        call opencdf (bulkname,N)
c*        if (may_day_flag.ne.0) return    !-->  EXIT
!
! If not opened yet, open forcing NetCDF file for reading.
! Find and save IDs for relevant variables, determine whether
! SST is a field or scalar value.
!
        if (ncidrad.eq.-1) then
          ierr=nf_open(radname(1:lstr), nf_nowrite, ncidrad)
          if (ierr .ne. nf_noerr) goto 4                  !--> ERROR
        endif
        ierr=nf_inq_varid(ncidrad, 'bulk_time', rad_tid)

        if (ierr .ne. nf_noerr) then
          write(stdout,3) 'bulk_time', radname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxRADLW))
        ierr=nf_inq_varid(ncidrad,vname(1,indxRADLW)(1:lvar),radlw_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidrad, radlw_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lradlwgrd=1
            else
              lradlwgrd=0
            endif
          endif
        endif

        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxRADLW)(1:lvar),radname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxShflx_rsw))
        ierr=nf_inq_varid(ncidrad,vname(1,indxShflx_rsw)
     &                                                (1:lvar),radsw_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidrad, radsw_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lradswgrd=1
            else
              lradswgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxShflx_rsw)(1:lvar),
     &                                                radname(1:lstr)
          goto 99                                        !--> ERROR
        endif

# ifdef DIURNAL_INPUT_SRFLX
        lvar=lenstr(vname(1,indxShflx_rswbio))
        ierr=nf_inq_varid (ncidrad,vname(1,indxShflx_rswbio)(1:lvar)
     &                                              ,radswbio_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims (ncidrad, radswbio_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lradswbiogrd=1
            else
              lradswbiogrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxShflx_rswbio)(1:lvar),
     &                                              radname(1:lstr)
          goto 99                                         !--> ERROR
        endif
# endif /* DIURNAL_INPUT_SRFLX */


! Determine whether there is cycling to reuse the input data and
! find cycling period "bulk_cycle", set initial cycling index
! "wspd_ncycle" and record index "wspd_rec".
! Set initial value for time index "itbulk" and both time record
! bounds to large negative artificial values, so that it will
! trigger the logic in reading part below.
! Also set scale factor to convert input dQdSST from Watts/m2/Celsius
! to meter/second.
!
        call set_cycle(ncidrad,   rad_tid,    ntrad,
     &                 rad_cycle, rad_ncycle, rad_rec)

        if (may_day_flag.ne.0) return    !-->  EXIT
        itrad      =  2
        rad_time(1)=  -1.E+20
        rad_time(2)=  -1.E+20
        srf_scale  =  1d0/(rho0*Cp)
      endif
!
! Reading data from the forcing file:  Get out, if model time is
!======== ==== ==== === ======= =====  already within the interval
! set by the past and future data times. Otherwise flip the time
! index, increment record and cyclin indices and read a new portion
! of data. Repeat it until model time is between the two times from
! data.
!
  1   i=3-itrad
       cff=time+0.5*dt    ! Unit: second

        if (rad_time(i).le.cff .and. cff.lt.rad_time(itrad)) 
     &            return
        ierr=advance_cycle (rad_cycle,ntrad,rad_ncycle,rad_rec)

        if (ierr .ne. 0) goto 7                           !--> ERROR

        ! Get the bulk time record for current time
        ierr=nf_get_var1_FTYPE (ncidrad, rad_tid, rad_rec, cff)

        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'bulk_time', rad_rec
          goto 99                                         !--> ERROR
        endif

        rad_time(i)=cff*day2sec+rad_cycle*rad_ncycle

        if (rad_time(itrad).eq.-1.E+20) 
     &          rad_time(itrad)=rad_time(i)
!
! radlw
!
        if (lradlwgrd.eq.1) then
          ierr=nf_fread(radlwg(START_2D_ARRAY,i), ncidrad, radlw_id,
     &                                             rad_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE(ncidrad,radlw_id,rad_rec,radlwp(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'RADLW', rad_rec
          goto 99                                         !--> ERROR
        endif
!
! radsw
!
        if (lradswgrd.eq.1) then
          ierr=nf_fread(radswg(START_2D_ARRAY,i), ncidrad, radsw_id,
     &                                            rad_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE(ncidrad,radsw_id,rad_rec,radswp(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'RADSW', rad_rec
          goto 99                                         !--> ERROR
        endif
!
! radswbio
!
# ifdef DIURNAL_INPUT_SRFLX
        if (lradswbiogrd.eq.1) then
          ierr=nf_fread (radswbiog(START_2D_ARRAY,i), ncidrad,
     &                        radswbio_id, rad_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE (ncidrad,radswbio_id,
     &                        rad_rec, radswbiop(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'RADSWBIO', rad_rec
          goto 99                                         !--> ERROR
        endif
# endif /* DIURNAL_INPUT_SRFLX */

       itrad=i
       MPI_master_only write(stdout,'(6x,A,1x,A,1x,g12.4,1x,I4)') 
     &              'GET_RAD   --',
     &              'Read fields for bulk formula   for time =', cff 
# ifdef MPI
     &                                                      , mynode
# endif

      if (ntrad.gt.1) goto 1
      if (ntrad.eq.1) return
!
! Sort out error messages: The following portion of the code is
!===== === ===== ========= not accessed unless something goes wrong.
!
  3   format(/,' GET_RAD - ERROR: unable to find forcing variable',
     &                   ': ',a,/,11x,'in forcing NetCDF file: ',a)
  4   write(stdout,5) radname(1:lstr)
  5   format(/,' GET_RAD - ERROR: unable to open forcing NetCDF ',
     &                                                  'file: ',a)
      goto 99
  6   format(/,' GET_RAD - ERROR while reading variable: ',a,2x,
     &                                      ' at TIME index = ',i4)
  7   write(stdout,8) rad_rec, ntrad, radname(1:lstr), tdays,
     &                                     rad_time(itrad)*sec2day
  8   format(/,' GET_RAD - ERROR: requested time record ',I4,
     &       1x,'exeeds the last available',  /,  11x,'record ',I4,
     &       1x,'in forcing NetCDF file: ', a, /,   11x,'TDAYS = ',
     &               g12.4,2x,'last available BULK_TIME = ',g12.4)
  99  may_day_flag=2
      return
      end

      subroutine set_rad_tile (Istr,Iend,Jstr,Jend)
      USE forces, only: itrad,rad_ncycle,rad_rec,rad_tid,radlw_id,
     &    radsw_id, lradlwgrd,lradswgrd,radlw, radsw, rad_time,
     &    rad_cycle, radlwg, radswg, radlwp, radswp, srf_scale,
     &    srflx, srflxg, srflxp 
# ifdef DIURNAL_INPUT_SRFLX
     &   ,radswbio
     &   ,srflxbio
# endif /* DIURNAL_INPUT_SRFLX */

!
! Set-up rad data for current tile.
!
!!!! gc # define BULK_DATA
!!!! gc # ifdef BULK_WSTR
!!!! gc #  define SMFLUX_DATA
!!!! gc # endif
      implicit none
# include "param.h"
# include "scalars.h"
# include "grid.h"
      integer Istr,Iend,Jstr,Jend, i,j, it1,it2
      real cff,cff1,cff2, cff3,cff4
      real val1,val2,val3,val4,val5,val6,val7,val8
# ifdef DIURNAL_INPUT_SRFLX
      real val55
# endif /* DIURNAL_INPUT_SRFLX */
!
# include "compute_extended_bounds.h"
!
      it1   =  3-itrad
      it2   =  itrad
      cff   =  time+0.5*dt
      cff1  =  rad_time(it2)-cff
      cff2  =  cff-rad_time(it1)
!
!  Load time invariant 
!
      if (rad_cycle.lt.0.) then
        if (FIRST_TIME_STEP) then
          if (lradlwgrd.eq.1) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                radlw(i,j)=srf_scale*radlwg(i,j,itrad)
                radsw(i,j)=srf_scale*radswg(i,j,itrad)
                srflx(i,j)=radsw(i,j)
# ifdef DIURNAL_INPUT_SRFLX
                radswbio(i,j)=srf_scale*radswbiog(i,j,itrad)
                srflxbio(i,j)=radswbio(i,j)
# endif /* DIURNAL_INPUT_SRFLX */
              enddo
            enddo
          else
            val4=srf_scale*radlwp(itrad)
            val5=srf_scale*radswp(itrad)
# ifdef DIURNAL_INPUT_SRFLX
            val55=srf_scale*radswbiop(itrad)
# endif /* DIURNAL_INPUT_SRFLX */
            do j=JstrR,JendR
              do i=IstrR,IendR
                radlw(i,j)=val4
                radsw(i,j)=val5
                srflx(i,j)=val5
# ifdef DIURNAL_INPUT_SRFLX
                radswbio(i,j)=val55
                srflxbio(i,j)=val55
# endif /* DIURNAL_INPUT_SRFLX */
              enddo
            enddo
          endif
        endif
!
!  Time-interpolate SST and dQdSST from grided or point data.
!  Check that for the next time step [when time=time+dt] time+dt 
!  is still between wspd_time(it1) and wspd_time(it2); and if not,
!  set synchro_flag top signal that the new forcing data should be
!  read from the netCDF input file (master thread only).
!
      elseif (cff1.ge.0. .and. cff2.ge.0.) then
        if (ZEROTH_TILE .and. cff1.lt.dt) then
            synchro_flag     =.TRUE.
            synchro_flag_rad =.TRUE.
        endif
                                    !note cff order maters
        cff = srf_scale/(cff1+cff2)
        cff3= cff1*cff
        cff4= cff2*cff

        cff = 1d0/(cff1+cff2)
        cff1= cff1*cff
        cff2= cff2*cff

        if (lradlwgrd.eq.1) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              radlw(i,j)=cff3*radlwg(i,j,it1)+cff4*radlwg(i,j,it2)
              radsw(i,j)=cff3*radswg(i,j,it1)+cff4*radswg(i,j,it2)
              srflx(i,j)=radsw(i,j)
# ifdef DIURNAL_INPUT_SRFLX
              radswbio(i,j)=cff3*radswbiog(i,j,it1)
     &                     +cff4*radswbiog(i,j,it2)
              srflxbio(i,j)=radswbio(i,j)
# endif /* DIURNAL_INPUT_SRFLX */
            enddo
          enddo
        else
          val4=cff3*radlwp(it1)+cff4*radlwp(it2)
          val5=cff3*radswp(it1)+cff4*radswp(it2)
# ifdef DIURNAL_INPUT_SRFLX
          val55=cff3*radswbiop(it1)+cff4*radswbiop(it2)
# endif /* DIURNAL_INPUT_SRFLX */

          do j=JstrR,JendR
            do i=IstrR,IendR
              radlw(i,j)=val4
              radsw(i,j)=val5
              srflx(i,j)=val5
# ifdef DIURNAL_INPUT_SRFLX
              radswbio(i,j)=val55
              srflxbio(i,j)=val55
# endif /* DIURNAL_INPUT_SRFLX */
            enddo
          enddo
        endif

!  Unable to set-up SST and dQdSST:
!  Complain about the error and signal to quit.      
!
      else
        if (ZEROTH_TILE) then
          write(stdout,1) 'rad_time',tdays,rad_time(it2)*sec2day
  1       format(/,' SET_RAD - current model time exceeds ending',
     &           1x,'value for variable: ',a,/,11x,'TDAYS = ',g12.4,
     &           2x,'TEND = ',g12.4)
          may_day_flag=2
        endif
      endif
      return
      end
#else
      subroutine get_rad_empty
      return
      end
#endif /* BULK_FLUX */
