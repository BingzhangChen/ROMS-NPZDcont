#include "cppdefs.h"
#if defined BULK_FLUX

      subroutine get_wspd
      USE forces, only: itwspd,wspd_ncycle,wspd_rec,wspd_tid,wspd_id,
     &    lwspdgrd, wspd, wspd_time,
     &    wspd_cycle, wspdg, wspdp
!
!  Read in wind speed
!
!
      implicit none
# include "param.h"
# include "scalars.h"
# include "netcdf.inc"
# include "ncscrum.h"
      real    cff
      integer i,ierr, lstr,lvar,lenstr, nf_fread, advance_cycle
       
      wspdname  = trim(nameit)//'_wspd.nc'
!
! Initialization: Inquire about the contents of forcing NetCDF file:
!================  variables and dimensions. Check for consistency.
!
      if (may_day_flag.ne.0) return      !-->  EXIT

      if (itwspd.eq.0 .or. iic.eq.0) then
        lstr=lenstr(wspdname)
c*        call opencdf (bulkname,N)
c*        if (may_day_flag.ne.0) return    !-->  EXIT
!
! If not opened yet, open forcing NetCDF file for reading.
! Find and save IDs for relevant variables, determine whether
! SST is a field or scalar value.
!
        if (ncidwspd.eq.-1) then
          ierr=nf_open(wspdname(1:lstr), nf_nowrite, ncidwspd)
          if (ierr. ne. nf_noerr) goto 4                  !--> ERROR
        endif

        ierr=nf_inq_varid(ncidwspd, 'bulk_time', wspd_tid)

        if (ierr .ne. nf_noerr) then
          write(stdout,3) 'bulk_time', wspdname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxwspd))
        ierr=nf_inq_varid(ncidwspd, vname(1,indxwspd)(1:lvar), wspd_id)

        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims(ncidwspd, wspd_id, i) 
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lwspdgrd=1
            else
              lwspdgrd=0
            endif
          endif
        endif

        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxwspd)(1:lvar), wspdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
!
! Determine whether there is cycling to reuse the input data and
! find cycling period "bulk_cycle", set initial cycling index
! "wspd_ncycle" and record index "wspd_rec".
! Set initial value for time index "itbulk" and both time record
! bounds to large negative artificial values, so that it will
! trigger the logic in reading part below.
! Also set scale factor to convert input dQdSST from Watts/m2/Celsius
! to meter/second.
!
        call set_cycle (ncidwspd,   wspd_tid,    ntwspd,
     &                  wspd_cycle, wspd_ncycle, wspd_rec)

        if (may_day_flag.ne.0) return    !-->  EXIT
        itwspd      =2
        wspd_time(1)=-1.E+20
        wspd_time(2)=-1.E+20
      endif
!
! Reading data from the forcing file:  Get out, if model time is
!======== ==== ==== === ======= =====  already within the interval
! set by the past and future data times. Otherwise flip the time
! index, increment record and cyclin indices and read a new portion
! of data. Repeat it until model time is between the two times from
! data.
!
  1   i=3-itwspd
       cff=time+0.5*dt    ! Unit: second

        if (wspd_time(i).le.cff .and. cff.lt.wspd_time(itwspd)) 
     &            return
        ierr=advance_cycle (wspd_cycle,ntwspd,wspd_ncycle,wspd_rec)

        if (ierr .ne. 0) goto 7                           !--> ERROR

        ! Get the bulk time record for current time
        ierr=nf_get_var1_FTYPE (ncidwspd, wspd_tid, wspd_rec, cff)

        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'wspd_time', wspd_rec
          goto 99                                         !--> ERROR
        endif

        wspd_time(i)=cff*day2sec+wspd_cycle*wspd_ncycle

        if (wspd_time(itwspd).eq.-1.E+20) 
     &          wspd_time(itwspd)=wspd_time(i)
!
! wspd
!
        if (lwspdgrd.eq.1) then
          ierr=nf_fread(wspdg(START_2D_ARRAY,i), ncidwspd, wspd_id,
     &                                           wspd_rec, r2dvar)
        else
          ierr=nf_get_var1_FTYPE(ncidwspd, wspd_id, wspd_rec, wspdp(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'wspd', wspd_rec
          goto 99                                         !--> ERROR
        endif
       itwspd=i
       MPI_master_only write(stdout,'(6x,A,1x,A,1x,g12.4,1x,I4)') 
     &              'GET_wspd   --',
     &              'Read fields for bulk formula   for time =', cff 
# ifdef MPI
     &                                                      , mynode
# endif

      if (ntwspd.gt.1) goto 1
      if (ntwspd.eq.1) return
!
! Sort out error messages: The following portion of the code is
!===== === ===== ========= not accessed unless something goes wrong.
!
  3   format(/,' GET_wspd - ERROR: unable to find forcing variable',
     &                   ': ',a,/,11x,'in forcing NetCDF file: ',a)
  4   write(stdout,5) wspdname(1:lstr)
  5   format(/,' GET_wspd - ERROR: unable to open forcing NetCDF ',
     &                                                  'file: ',a)
      goto 99
  6   format(/,' GET_wspd - ERROR while reading variable: ',a,2x,
     &                                      ' at TIME index = ',i4)
  7   write(stdout,8) wspd_rec, ntwspd, wspdname(1:lstr), tdays,
     &                                     wspd_time(itwspd)*sec2day
  8   format(/,' GET_wspd - ERROR: requested time record ',I4,
     &       1x,'exeeds the last available',  /,  11x,'record ',I4,
     &       1x,'in forcing NetCDF file: ', a, /,   11x,'TDAYS = ',
     &               g12.4,2x,'last available BULK_TIME = ',g12.4)
  99  may_day_flag=2
      return
      end subroutine


      subroutine set_wspd_tile(Istr,Iend,Jstr,Jend)
      USE forces, only: itwspd,wspd_ncycle,wspd_rec,wspd_tid,wspd_id,
     &    lwspdgrd, wspd, wspd_time,
     &    wspd_cycle, wspdg, wspdp

!
! Set-up wspd data for current tile.
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
      real val1
!
# include "compute_extended_bounds.h"
!
      it1   =3-itwspd
      it2   =itwspd
      cff   =time+0.5*dt
      cff1  =wspd_time(it2)-cff
      cff2  =cff-wspd_time(it1)
!
!  Load time invariant 
!
      if (wspd_cycle.lt.0.) then
        if (FIRST_TIME_STEP) then
          if (lwspdgrd.eq.1) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                wspd(i,j)=wspdg(i,j,itwspd)
              enddo
            enddo
          else
            val1=wspdp(itwspd)
            do j=JstrR,JendR
              do i=IstrR,IendR
                wspd(i,j)=val1
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
            synchro_flag_wspd=.TRUE.
        endif
                                    !note cff order maters
        cff =1d0/(cff1+cff2)
        cff1=cff1*cff
        cff2=cff2*cff

        if (lwspdgrd.eq.1) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              wspd(i,j)=cff1*wspdg(i,j,it1)+cff2*wspdg(i,j,it2)
            enddo
          enddo
        else
          val1=cff1*wspdp(it1)+cff2*wspdp(it2)

          do j=JstrR,JendR
            do i=IstrR,IendR
              wspd(i,j)=val1
            enddo
          enddo
        endif
!
!  Unable to set-up SST and dQdSST:
!  Complain about the error and signal to quit.      
!
      else
        if (ZEROTH_TILE) then
          write(stdout,1) 'wspd_time',tdays,wspd_time(it2)*sec2day
  1       format(/,' SET_wspd - current model time exceeds ending',
     &           1x,'value for variable: ',a,/,11x,'TDAYS = ',g12.4,
     &           2x,'TEND = ',g12.4)
          may_day_flag=2
        endif
      endif
      return
      end
#else
      subroutine get_wspd_empty
      return
      end
#endif /* BULK_FLUX */
