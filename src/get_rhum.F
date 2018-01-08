#include "cppdefs.h"
#if defined BULK_FLUX

      subroutine get_rhum
!
!  Read in specific humidity
!
!
      use forces, only: itrhum,rhum_ncycle,rhum_rec,rhum_tid,rhum_id,
     &    lrhumgrd, rhum, rhum_time,
     &    rhum_cycle, rhumg, rhump

      implicit none
# include "param.h"
# include "scalars.h"
# include "netcdf.inc"
# include "ncscrum.h"
      real    cff
      integer i,ierr, lstr,lvar,lenstr, nf_fread, advance_cycle
       
      rhumname  = trim(nameit)//'_rhum.nc'
!
! Initialization: Inquire about the contents of forcing NetCDF file:
!================  variables and dimensions. Check for consistency.
!
      if (may_day_flag.ne.0) return      !-->  EXIT

      if (itrhum.eq.0 .or. iic.eq.0) then
        lstr=lenstr(rhumname)
c*        call opencdf (bulkname,N)
c*        if (may_day_flag.ne.0) return    !-->  EXIT
!
! If not opened yet, open forcing NetCDF file for reading.
! Find and save IDs for relevant variables, determine whether
! SST is a field or scalar value.
!
        if (ncidrhum.eq.-1) then
          ierr=nf_open(rhumname(1:lstr), nf_nowrite, ncidrhum)
          if (ierr. ne. nf_noerr) goto 4                  !--> ERROR
        endif

        ierr=nf_inq_varid(ncidrhum,'bulk_time', rhum_tid)

        if (ierr .ne. nf_noerr) then
          write(stdout,3) 'bulk_time', rhumname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        ! Inquire the dimension of the data of rhum:
        lvar=lenstr(vname(1,indxRHUM))
        ierr=nf_inq_varid (ncidrhum, vname(1,indxRHUM)(1:lvar),rhum_id)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims(ncidrhum, rhum_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lrhumgrd=1
            else
              lrhumgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxRHUM)(1:lvar),rhumname(1:lstr)
          goto 99                                         !--> ERROR
        endif

! Determine whether there is cycling to reuse the input data and
! find cycling period "bulk_cycle", set initial cycling index
! "wspd_ncycle" and record index "wspd_rec".
! Set initial value for time index "itbulk" and both time record
! bounds to large negative artificial values, so that it will
! trigger the logic in reading part below.
! Also set scale factor to convert input dQdSST from Watts/m2/Celsius
! to meter/second.
!
        call set_cycle (ncidrhum,   rhum_tid,    ntrhum,
     &                  rhum_cycle, rhum_ncycle, rhum_rec)
        if (may_day_flag.ne.0) return    !-->  EXIT
        itrhum=2
        rhum_time(1)=-1.E+20
        rhum_time(2)=-1.E+20
      endif
!
! Reading data from the forcing file:  Get out, if model time is
!======== ==== ==== === ======= =====  already within the interval
! set by the past and future data times. Otherwise flip the time
! index, increment record and cyclin indices and read a new portion
! of data. Repeat it until model time is between the two times from
! data.
!
  1   i=3-itrhum
       cff=time+0.5*dt    ! Unit: second

        if (rhum_time(i).le.cff .and. cff.lt.rhum_time(itrhum)) 
     &            return
        ierr=advance_cycle(rhum_cycle,ntrhum,rhum_ncycle,rhum_rec)

        if (ierr .ne. 0) goto 7                           !--> ERROR

        ! Get the bulk time record for current time
        ierr=nf_get_var1_FTYPE(ncidrhum, rhum_tid, rhum_rec, cff)

        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'bulk_time', rhum_rec
          goto 99                                         !--> ERROR
        endif
       
        ! Update bulk time:
        rhum_time(i)=cff*day2sec+rhum_cycle*rhum_ncycle

        if (rhum_time(itrhum).eq.-1.E+20) 
     &          rhum_time(itrhum)=rhum_time(i)
!
! rhum
!
        if (lrhumgrd.eq.1) then
          ierr=nf_fread(rhumg(START_2D_ARRAY,i), ncidrhum, rhum_id,
     &                                            rhum_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE(ncidrhum,rhum_id,rhum_rec,rhump(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'RHUM', rhum_rec
          goto 99                                         !--> ERROR
        endif

       itrhum=i
       MPI_master_only write(stdout,'(6x,A,1x,A,1x,g12.4,1x,I4)') 
     &              'GET_RHUM   --',
     &              'Read fields for bulk formula for time =', cff 
# ifdef MPI
     &                                                      , mynode
# endif

      if (ntrhum.gt.1) goto 1
      if (ntrhum.eq.1) return
!
! Sort out error messages: The following portion of the code is
!===== === ===== ========= not accessed unless something goes wrong.
!
  3   format(/,' GET_BULK - ERROR: unable to find forcing variable',
     &                   ': ',a,/,11x,'in forcing NetCDF file: ',a)
  4   write(stdout,5) rhumname(1:lstr)
  5   format(/,' GET_BULK - ERROR: unable to open forcing NetCDF ',
     &                                                  'file: ',a)
      goto 99
  6   format(/,' GET_BULK - ERROR while reading variable: ',a,2x,
     &                                      ' at TIME index = ',i4)
  7   write(stdout,8) rhum_rec, ntrhum, rhumname(1:lstr), tdays,
     &                                     rhum_time(itrhum)*sec2day
  8   format(/,' GET_BULK - ERROR: requested time record ',I4,
     &       1x,'exeeds the last available',  /,  11x,'record ',I4,
     &       1x,'in forcing NetCDF file: ', a, /,   11x,'TDAYS = ',
     &               g12.4,2x,'last available BULK_TIME = ',g12.4)
  99  may_day_flag=2
      return
      end

      subroutine set_rhum_tile (Istr,Iend,Jstr,Jend)
      USE forces, only: itrhum,rhum_ncycle,rhum_rec,rhum_tid,rhum_id,
     &    lrhumgrd, rhum, rhum_time,
     &    rhum_cycle, rhumg, rhump
!
! Set-up rhum data for current tile.
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
      real val3
# include "compute_extended_bounds.h"
      it1=3-itrhum
      it2=itrhum
      cff=time+0.5*dt
      cff1=rhum_time(it2)-cff
      cff2=cff-rhum_time(it1)
!
!  Load time invariant 
!
      if (rhum_cycle.lt.0.) then
        if (FIRST_TIME_STEP) then
          if (lrhumgrd.eq.1) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                rhum(i,j)=rhumg(i,j,itrhum)
              enddo
            enddo
          else
            val3=rhump(itrhum)

            do j=JstrR,JendR
              do i=IstrR,IendR
                rhum(i,j)=val3
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
           synchro_flag_rhum=.TRUE.
        endif
                                    !note cff order maters
        cff=1./(cff1+cff2)
        cff1=cff1*cff
        cff2=cff2*cff

        if (lrhumgrd.eq.1) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              rhum(i,j)=cff1*rhumg(i,j,it1)+cff2*rhumg(i,j,it2)
            enddo
          enddo
        else
          val3=cff1*rhump(it1)+cff2*rhump(it2)

          do j=JstrR,JendR
            do i=IstrR,IendR
              rhum(i,j)=val3
            enddo
          enddo
        endif
!
!  Unable to set-up SST and dQdSST:
!  Complain about the error and signal to quit.      
!
      else
        if (ZEROTH_TILE) then
          write(stdout,1) 'rhum_time',tdays,rhum_time(it2)*sec2day
  1       format(/,' SET_BULK - current model time exceeds ending',
     &           1x,'value for variable: ',a,/,11x,'TDAYS = ',g12.4,
     &           2x,'TEND = ',g12.4)
          may_day_flag=2
        endif
      endif
      return
      end
#else
      subroutine get_rhum_empty
      return
      end
#endif /* BULK_FLUX */
