#include "cppdefs.h"
#if defined BULK_FLUX

      subroutine get_PRATE
!
!  Read in wind speed and surface air temperature
!
!
      use forces, only: itPRATE,PRATE_ncycle,PRATE_rec,PRATE_tid,
     &    PRATE_id,
     &    lPRATEgrd,   PRATE_time,
     &    PRATE_cycle, PRATEg, PRATEp,stf_scale
      implicit none
# include "param.h"
# include "scalars.h"
# include "netcdf.inc"
# include "ncscrum.h"
      real    cff
      integer i,ierr, lstr,lvar,lenstr, nf_fread, advance_cycle
       
      PRATENAME = trim(nameit)//'_prate.nc'
!
! Initialization: Inquire about the contents of forcing NetCDF file:
!================  variables and dimensions. Check for consistency.
!
      if (may_day_flag.ne.0) return      !-->  EXIT

# ifdef BULK_EP
      if (itPRATE.eq.0 .or. iic.eq.0) then
        lstr=lenstr(PRATENAME)
c*        call opencdf (bulkname,N)
c*        if (may_day_flag.ne.0) return    !-->  EXIT
!
! If not opened yet, open forcing NetCDF file for reading.
! Find and save IDs for relevant variables, determine whether
! SST is a field or scalar value.
!
        if (ncidPRATE.eq.-1) then
          ierr=nf_open(PRATENAME(1:lstr), nf_nowrite, ncidPRATE)
          if (ierr. ne. nf_noerr) goto 4                  !--> ERROR
        endif

        ierr=nf_inq_varid (ncidPRATE, 'bulk_time', PRATE_tid)

        if (ierr .ne. nf_noerr) then
          write(stdout,3) 'bulk_time', PRATENAME(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxPRATE))
        ierr=nf_inq_varid(ncidPRATE,vname(1,indxPRATE)(1:lvar),prate_id)

        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_varndims(ncidPRATE, prate_id, i)

          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              lprategrd=1
            else
              lprategrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxPRATE)(1:lvar),PRATEname(1:lstr)
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
        call set_cycle (ncidPRATE,   PRATE_tid,    ntPRATE,
     &                  PRATE_cycle, PRATE_ncycle, PRATE_rec)
        if (may_day_flag.ne.0) return    !-->  EXIT
        itPRATE         = 2
        PRATE_time(1)   = -1.E+20
        PRATE_time(2)   = -1.E+20
        stf_scale(isalt)= 1D-2/864D2
      endif
!
! Reading data from the forcing file:  Get out, if model time is
!======== ==== ==== === ======= =====  already within the interval
! set by the past and future data times. Otherwise flip the time
! index, increment record and cyclin indices and read a new portion
! of data. Repeat it until model time is between the two times from
! data.
!
  1   i   =3   -itPRATE
       cff=time+0.5*dt    ! Unit: second

        if (PRATE_time(i).le.cff .and. cff.lt.PRATE_time(itPRATE)) 
     &            return
        ierr=advance_cycle (PRATE_cycle,ntPRATE,PRATE_ncycle,PRATE_rec)

        if (ierr .ne. 0) goto 7                           !--> ERROR

        ! Get the bulk time record for current time
        ierr=nf_get_var1_FTYPE (ncidPRATE, PRATE_tid, PRATE_rec, cff)

        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'bulk_time', PRATE_rec
          goto 99                                         !--> ERROR
        endif

        PRATE_time(i)=cff*day2sec+PRATE_cycle*PRATE_ncycle

        if (PRATE_time(itPRATE).eq.-1.E+20) 
     &          PRATE_time(itPRATE)=PRATE_time(i)
!
! prate
!
        if (lprategrd.eq.1) then
          ierr=nf_fread (prateg(START_2D_ARRAY,i), ncidPRATE, prate_id,
     &                                            PRATE_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE(ncidPRATE,prate_id,PRATE_rec,pratep(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'PRATE', PRATE_rec
          goto 99                                         !--> ERROR
        endif

       itPRATE=i
       MPI_master_only write(stdout,'(6x,A,1x,A,1x,g12.4,1x,I4)') 
     &              'GET_PRATE   --',
     &              'Read fields for bulk formula   for time =', cff 
# ifdef MPI
     &                                                      , mynode
# endif

      if (ntPRATE.gt.1) goto 1
      if (ntPRATE.eq.1) return
!
! Sort out error messages: The following portion of the code is
!===== === ===== ========= not accessed unless something goes wrong.
!
  3   format(/,' GET_PRATE - ERROR: unable to find forcing variable',
     &                   ': ',a,/,11x,'in forcing NetCDF file: ',a)
  4   write(stdout,5) PRATEname(1:lstr)
  5   format(/,' GET_PRATE - ERROR: unable to open forcing NetCDF ',
     &                                                  'file: ',a)
      goto 99
  6   format(/,' GET_PRATE - ERROR while reading variable: ',a,2x,
     &                                      ' at TIME index = ',i4)
  7   write(stdout,8) PRATE_rec, ntPRATE, PRATEname(1:lstr), tdays,
     &                                     PRATE_time(itPRATE)*sec2day
  8   format(/,' GET_PRATE - ERROR: requested time record ',I4,
     &       1x,'exeeds the last available',  /,  11x,'record ',I4,
     &       1x,'in forcing NetCDF file: ', a, /,   11x,'TDAYS = ',
     &               g12.4,2x,'last available BULK_TIME = ',g12.4)
# endif /* BULK_EP */
  99  may_day_flag=2
      return
      end


      subroutine set_PRATE_tile(Istr,Iend,Jstr,Jend)
      use forces, only: itPRATE,PRATE_ncycle,PRATE_rec,PRATE_tid,
     &    PRATE_id, prate,
     &    lPRATEgrd, PRATE_time,
     &    PRATE_cycle, PRATEg, PRATEp,stf_scale

!
! Set-up PRATE data for current tile.
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
      real cff5,cff6
      real val6
!
# include "compute_extended_bounds.h"
!
      it1   =  3-itPRATE
      it2   =  itPRATE
      cff   =  time+0.5*dt
      cff1  =  PRATE_time(it2)-cff
      cff2  =  cff-PRATE_time(it1)
!
!  Load time invariant 
!
      if (PRATE_cycle.lt.0.) then
        if (FIRST_TIME_STEP) then
          if (lprategrd.eq.1) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                prate(i,j)=stf_scale(isalt)*prateg(i,j,itPRATE)
              enddo
            enddo
          else
            val6=stf_scale(isalt)*pratep(itPRATE)
            do j=JstrR,JendR
              do i=IstrR,IendR
                prate(i,j)=val6
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
           synchro_flag      =.TRUE.
           synchro_flag_prate=.TRUE.
        endif
                                    !note cff order maters
        cff =stf_scale(isalt)/(cff1+cff2)
        cff5=cff1*cff
        cff6=cff2*cff

        if (lprategrd.eq.1) then
          do j=JstrR,JendR
            do i=IstrR,IendR
# ifdef BULK_EP
              prate(i,j)=cff5*prateg(i,j,it1)+cff6*prateg(i,j,it2)
# endif
            enddo
          enddo
        else
# ifdef BULK_EP
          val6=cff5*pratep(it1)+cff6*pratep(it2)
# endif
          do j=JstrR,JendR
            do i=IstrR,IendR
# ifdef BULK_EP
              prate(i,j)=val6
# endif
            enddo
          enddo
        endif
!
!  Unable to set-up SST and dQdSST:
!  Complain about the error and signal to quit.      
!
      else
        if (ZEROTH_TILE) then
          write(stdout,1) 'PRATE_time',tdays,PRATE_time(it2)*sec2day
  1       format(/,' SET_PRATE - current model time exceeds ending',
     &           1x,'value for variable: ',a,/,11x,'TDAYS = ',g12.4,
     &           2x,'TEND = ',g12.4)
          may_day_flag=2
        endif
      endif
      return
      end
#else
      subroutine get_PRATE_empty
      return
      end
#endif /* BULK_FLUX */
