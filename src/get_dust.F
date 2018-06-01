#include "cppdefs.h"
#if defined BIOLOGY && (defined PZEFTFE || defined BIO_NPZDDISC)

                              ! Read point or grided 
      subroutine get_dust     ! Dust deposition data
                              ! at the appropriate
                              ! time from forcing NetCDF file.
      use forces
# define DUST_DATA
      implicit none
!# include "param.h"
!# include "forces.h"
# include "scalars.h"
# include "ncscrum.h"
      real cff
      integer i,ierr,lstr,lvar, lenstr, nf_fread, advance_cycle
#include "netcdf.inc"
!
! Initialization: Inquire about the contents of forcing NetCDF file:
!================  variables and dimensions. Check for consistency.
!
      if (may_day_flag .ne. 0) return      !-->  EXIT

      IF (itdust.eq.0 .OR. iic.eq.0) THEN

!Read the Dust deposition data:
        lstr = lenstr(frcname)
        if (may_day_flag .ne. 0) return    !-->  EXIT

        ! Test whether the Dust forcing NetCDF file has already been opened:
        if (ncidfrc .eq. -1) then
          ierr   = nf_open (frcname(1:lstr), nf_nowrite, ncidfrc)
          if (ierr .ne. nf_noerr) then
             write(stdout,14) frcname
          endif
        endif

        ierr = nf_inq_varid(ncidfrc,  'dust_time', dust_tid)
        if (ierr .ne. nf_noerr) then
          write(stdout,15) 'dust_time', frcname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        ierr = nf_inq_varid (ncidfrc, 'dust', dust_id)

        if (ierr .eq. nf_noerr) then
          ! Obtain the dimension of dust (i)
          ierr = nf_inq_varndims (ncidfrc, dust_id, i)
          if (ierr. eq. nf_noerr) then
            if (i .gt. 1) then
              ldustgrd=1
            else
              ldustgrd=0
            endif
          endif
        endif

        if (ierr .ne. nf_noerr) then
          write(stdout,15) 'dust', frcname(1:lstr)
          goto 99                                         !--> ERROR
        endif

!
! Determine whether there is cycling to reuse the input data
! and find cycling period "dust_cycle", set initial cycling
! index "srf_ncycle" and record index "dust_rec".
! Set initial value for time index "itdust" and both time record
! bounds to large negative artificial values, so that it will
! trigger the logic in reading part below.

        call set_cycle (ncidfrc,   dust_tid,   ntdust,
     &                  dust_cycle,dust_ncycle,dust_rec)

        if (may_day_flag.ne.0) return    !-->  EXIT

        itdust      = 2
        dust_time(1)= -1.E+20
        dust_time(2)= -1.E+20

      ENDIF                   !<-- itdust.eq.0 .or. iic.eq.0

!
! Reading data from the forcing file:  Get out, if model time is
!======== ==== ==== === ======= =====  already within the interval
! set by the past and future data times. Otherwise flip the time
! index, increment record and cyclin indices and read a new portion
! of data. Repeat it until model time is between the two times from
! data.
!
  1   i = 3 - itdust

        cff = time + 0.5*dt  ! The model time

        if (dust_time(i).le.cff .and. cff.lt.dust_time(itdust)) return
        !MPI_master_only write(stdout,*) 'dust_cycle=', dust_cycle
        !MPI_master_only write(stdout,*) 'dust_ncycle=', dust_ncycle
        !MPI_master_only write(stdout,*) 'ntdust=', ntdust
        !MPI_master_only write(stdout,*) 'dust_rec=', dust_rec
        ierr = advance_cycle(dust_cycle,ntdust,dust_ncycle,dust_rec)

        if (ierr.ne.0) goto 7                           !--> ERROR

        !Read the dust_time data from the time record dust_rec into cff  
        ierr=nf_get_var1_FTYPE(ncidfrc, dust_tid, dust_rec, cff)

        if (ierr .ne. nf_noerr) then
          write(stdout,16) 'dust_time', dust_rec 
          goto 99                                         !--> ERROR
        endif

        dust_time(i) = cff*day2sec + dust_cycle*dust_ncycle

        if (dust_time(itdust).eq.-1.E+20) dust_time(itdust)=dust_time(i)
!
!Dust
        if (ldustgrd .eq. 1) then
          ierr=nf_fread(dustg(START_2D_ARRAY,i), ncidfrc, dust_id,
     &                                             dust_rec, r2dvar)

        else
          ierr=nf_get_var1_FTYPE (ncidfrc,dust_id,dust_rec,dustp(i))
        endif

        if (ierr .ne. nf_noerr) then
          write(stdout,16) 'Dust', dust_rec
          goto 99                                         !--> ERROR
        endif 

       itdust=i
       MPI_master_only write(stdout,'(6x,A,1x,A,1x,g12.4,1x,I4)') 
     &           'GET_DUST --',
     &           'Read DUST deposition data for time =', cff
#ifdef MPI
     &                                                      , mynode
#endif

      if (ntdust.gt.1) goto 1
      if (ntdust.eq.1) return 
!
! Sort out error messages: The following portion of the code is
!===== === ===== ========= not accessed unless something goes wrong.
!
  7   write(stdout,8) dust_rec, ntdust, frcname(1:lstr), tdays,
     &                                dust_time(itdust)*sec2day

  8   format(/,' GET_DUST - ERROR: requested time record ',I4,
     &       1x,'exeeds the last available', /,14x,'record ',I4,
     &       1x,'in forcing NetCDF file: ',a, /,14x,'TDAYS = ',
     &       g12.4,2x,'last available DUST_TIME = ',g12.4)
  
 14     format(/,' get_dust - unable to open forcing netCDF ',1x,A)
 15     format(/,' get_dust - unable to find forcing variable: ',A,
     &                               /,14x,'in forcing netCDF  ',A)
 16     format(/,' get_dust - error while reading variable: ',A,2x,
     &                                           ' at TIME index = ',i4)

  99  may_day_flag=2
      return
      end

      subroutine set_dust_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set-up dust deposition flux for current tile.
!--------------------------------------------------------------------
!
      use forces
# define DUST_DATA
      implicit none
!# include "param.h"
!# include "forces.h"
# include "scalars.h"
# include "grid.h"
      integer Istr,Iend,Jstr,Jend, i,j, it1,it2
      real    cff, cff1, cff2,dst
!
# include "compute_extended_bounds.h"
!
! Set coefficients for interpolation. Check that for the next time
! step [when time=time+dt] both weights will still be positive, and
! if not, set synchro_flag to signal that new data should be read
! from an appropriate netCDF input file (master thread only).
!
! After that either load time-invariant data, or interpolate in time
! or complain about error and signal to quit, if interpolation is
! needed, but not possible.
!
!
      it1   = 3-itdust
      it2   = itdust
      cff   = time + 0.5*dt
      cff1  = dust_time(it2) - cff
      cff2  = cff            - dust_time(it1)

      if (ZEROTH_TILE .and. cff1.lt.dt) synchro_flag=.TRUE.
!
! Load time-invariant shortwave radiation flux.
! Time interpolation is not performed in this case.
!
      IF (dust_cycle.lt.0.) THEN
        if (FIRST_TIME_STEP) then
          if (ldustgrd.eq.1) then
            do j = JstrR,JendR
              do i=IstrR,IendR
                dust(i,j)    = dustg(i,j,itdust)
               enddo
            enddo
          else
            dst  = dustp(itdust)
            do j = JstrR,JendR
              do i=IstrR,IendR
                dust(i,j)    = dst
              enddo
            enddo
          endif
        endif
!
! Time-interpolate dust deposition flux from grided or point
! data. Check that for the next time step [when time=time+dt]
! time+dt is still between srf_time(it1) and srf_time(it2);
! and if not, set synchro_flag to signal that the new forcing data
! should be read from the netCDF input file (master thread only).
!
      ELSEIF (cff1.GE.0. .AND. cff2.GE.0.) THEN
        cff  =  1./(cff1+cff2)
        cff1 =  cff1*cff
        cff2 =  cff2*cff
        if (ldustgrd.eq.1) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              dust(i,j) = cff1*dustg(i,j,it1) + cff2*dustg(i,j,it2)
            enddo
          enddo
        else
          dst = cff1*dustp(it1) + cff2*dustp(it2)
            do j=JstrR,JendR
              do i=IstrR,IendR
              dust(i,j) = dst
            enddo
          enddo
        endif
!
! Unable to set-up dust deposition flux:
! Complain about the error and signal to quit (ONE THREAD ONLY).
!
      ELSE
        if (ZEROTH_TILE) then
          write(stdout,1) 'dust_time', tdays, dust_time(it2)*sec2day  
  1       format(/,' SET_dust_TILE - current model time exceeds',
     &             ' ending value for variable: ',a,/,14x,'TDAYS = ',
     &                                      g12.4,2x,'TEND = ',g12.4)
          may_day_flag=2
        endif
      ENDIF

      return
      end
#else
      subroutine get_dust_empty
      return
      end
#endif /* BIOLOGY && PZEFTFE  */
