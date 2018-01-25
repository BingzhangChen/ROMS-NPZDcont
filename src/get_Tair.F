#include "cppdefs.h"
#if defined BULK_FLUX

      subroutine get_Tair
!
!  Read in surface air temperature
!
      use forces, only: ittair,tair_ncycle,tair_rec,tair_tid,tair_id,
     &    ltairgrd, tair_time,
     &    tair_cycle, tairg, tairp, tair

      implicit none
# include "param.h"
# include "scalars.h"
# include "netcdf.inc"
# include "ncscrum.h"
      real    cff
      integer i,ierr, lstr,lvar,lenstr, nf_fread, advance_cycle
       
      Tairname  = trim(nameit)//'_Tair.nc'
!
! Initialization: Inquire about the contents of forcing NetCDF file:
!================  variables and dimensions. Check for consistency.
!
      if (may_day_flag.ne.0) return      !-->  EXIT

      if (itTair.eq.0 .or. iic.eq.0) then
        lstr=lenstr(Tairname)
c*        call opencdf (bulkname,N)
c*        if (may_day_flag.ne.0) return    !-->  EXIT
!
! If not opened yet, open forcing NetCDF file for reading.
! Find and save IDs for relevant variables, determine whether
! SST is a field or scalar value.
!
        if (ncidTair.eq.-1) then
          ierr=nf_open(Tairname(1:lstr), nf_nowrite, ncidTair)
          if (ierr. ne. nf_noerr) goto 4                  !--> ERROR
        endif

        ierr=nf_inq_varid (ncidTair, 'bulk_time', Tair_tid)

        if (ierr .ne. nf_noerr) then
          write(stdout,3) 'bulk_time', Tairname(1:lstr)
          goto 99                                         !--> ERROR
        endif

        lvar=lenstr(vname(1,indxTAIR))
        ierr=nf_inq_varid (ncidTair, vname(1,indxTAIR)(1:lvar),tair_id)
        if (ierr .eq. nf_noerr) then
          ! Inquire the data dimension of tair
          ierr=nf_inq_varndims (ncidTair, tair_id, i)
          if (ierr. eq. nf_noerr) then
            if (i.gt.1) then
              ltairgrd=1
            else
              ltairgrd=0
            endif
          endif
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,3) vname(1,indxTAIR)(1:lvar),Tairname(1:lstr)
          goto 99                                         !--> ERROR
        endif
!
! Determine whether there is cycling to reuse the input data and
! find cycling period "tair_cycle", set initial cycling index
! "tair_ncycle" and record index "tair_rec".
! Set initial value for time index "ittair" and both time record
! bounds to large negative artificial values, so that it will
! trigger the logic in reading part below.
! Also set scale factor to convert input dQdSST from Watts/m2/Celsius
! to meter/second.
!
        call set_cycle(ncidtair,   tair_tid,    nttair,
     &                  tair_cycle, tair_ncycle, tair_rec)
        if (may_day_flag.ne.0) return    !-->  EXIT
        ittair      =2
        tair_time(1)=-1.E+20
        tair_time(2)=-1.E+20
      endif
!
! Reading data from the forcing file:  Get out, if model time is
!======== ==== ==== === ======= =====  already within the interval
! set by the past and future data times. Otherwise flip the time
! index, increment record and cyclin indices and read a new portion
! of data. Repeat it until model time is between the two times from
! data.
!
  1   i=3-ittair
       cff=time+0.5*dt    ! Unit: second

        if (tair_time(i).le.cff .and. cff.lt.tair_time(ittair)) 
     &            return
        ierr=advance_cycle(tair_cycle,nttair,tair_ncycle,tair_rec)

        !if(mynode==0) write(6,*) 'Bulk record:',tair_rec

        if (ierr .ne. 0) goto 7                           !--> ERROR

        ! Get the tair time record for current time
        ierr=nf_get_var1_FTYPE (ncidtair, tair_tid, tair_rec, cff)
        !if(mynode==0) write(6,*) 'Current tair time (days) =',cff

        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'tair_time', tair_rec
          goto 99                                         !--> ERROR
        endif

        tair_time(i)=cff*day2sec+tair_cycle*tair_ncycle
        ! if(mynode==0) write(6,*) 'Updated tair time (sec)=',tair_time(i)

        if (tair_time(ittair).eq.-1.E+20) 
     &          tair_time(ittair)=tair_time(i)
!
! tair
!
        if (ltairgrd.eq.1) then
          ierr=nf_fread(tairg(START_2D_ARRAY,i), ncidtair, tair_id,
     &                                            tair_rec,  r2dvar)
        else
          ierr=nf_get_var1_FTYPE(ncidtair,tair_id,tair_rec,tairp(i))
        endif
        if (ierr .ne. nf_noerr) then
          write(stdout,6) 'TAIR', tair_rec
          goto 99                                         !--> ERROR
        endif
       ittair=i
       MPI_master_only write(stdout,'(6x,A,1x,A,1x,g12.4,1x,I4)') 
     &              'GET_TAIR   --',
     &              'Read fields for bulk formula   for time =', cff 
# ifdef MPI
     &                                                      , mynode
# endif

      if (nttair.gt.1) goto 1
      if (nttair.eq.1) return
!
! Sort out error messages: The following portion of the code is
!===== === ===== ========= not accessed unless something goes wrong.
!
  3   format(/,' GET_TAIR - ERROR: unable to find forcing variable',
     &                   ': ',a,/,11x,'in forcing NetCDF file: ',a)
  4   write(stdout,5) tairname(1:lstr)
  5   format(/,' GET_TAIR - ERROR: unable to open forcing NetCDF ',
     &                                                  'file: ',a)
      goto 99
  6   format(/,' GET_TAIR - ERROR while reading variable: ',a,2x,
     &                                      ' at TIME index = ',i4)
  7   write(stdout,8) tair_rec, nttair, tairname(1:lstr), tdays,
     &                                     tair_time(ittair)*sec2day
  8   format(/,' GET_TAIR - ERROR: requested time record ',I4,
     &       1x,'exeeds the last available',  /,  11x,'record ',I4,
     &       1x,'in forcing NetCDF file: ', a, /,   11x,'TDAYS = ',
     &               g12.4,2x,'last available TAIR_TIME = ',g12.4)
  99  may_day_flag=2
      return
      end

      subroutine set_tair_tile (Istr,Iend,Jstr,Jend)
!
! Set-up tair data for current tile.
!
!!!! gc # define BULK_DATA
!!!! gc # ifdef BULK_WSTR
!!!! gc #  define SMFLUX_DATA
!!!! gc # endif
      USE forces, only: ittair,tair_ncycle,tair_rec,tair_tid,tair_id,
     &    ltairgrd, tair, tair_time,
     &    tair_cycle, tairg, tairp

      implicit none
# include "param.h"
!# include "forces.h"
# include "scalars.h"
# include "grid.h"
      integer Istr,Iend,Jstr,Jend, i,j, it1,it2
      real cff,cff1,cff2, cff3,cff4

      real val2
# include "compute_extended_bounds.h"
!
      it1 =3-ittair
      it2 =ittair
      cff =time+0.5*dt
      cff1=tair_time(it2)-cff
      cff2=cff-tair_time(it1)
!
!  Load time invariant 
!
      if (tair_cycle.lt.0.) then
        if (FIRST_TIME_STEP) then
          if (ltairgrd.eq.1) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                tair(i,j)=tairg(i,j,ittair)
              enddo
            enddo
          else
            val2=tairp(ittair)
            do j=JstrR,JendR
              do i=IstrR,IendR
                tair(i,j)=val2
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
           synchro_flag=.TRUE.
           synchro_flag_tair=.TRUE.
        endif
                                    !note cff order maters
        cff =1d0/(cff1+cff2)
        cff1=cff1*cff
        cff2=cff2*cff

        if (ltairgrd.eq.1) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              tair(i,j)=cff1*tairg(i,j,it1)+cff2*tairg(i,j,it2)
            enddo
          enddo
        else
          val2=cff1*tairp(it1)+cff2*tairp(it2)
          do j=JstrR,JendR
            do i=IstrR,IendR
              tair(i,j)=val2
            enddo
          enddo
        endif
!
!  Unable to set-up SST and dQdSST:
!  Complain about the error and signal to quit.      
!
      else
        if (ZEROTH_TILE) then
          write(stdout,1) 'tair_time',tdays,tair_time(it2)*sec2day
  1       format(/,' SET_TAIR - current model time exceeds ending',
     &           1x,'value for variable: ',a,/,11x,'TDAYS = ',g12.4,
     &           2x,'TEND = ',g12.4)
          may_day_flag=2
        endif
      endif
      return
      end
#else
      subroutine get_tair_empty
      return
      end
#endif /* BULK_FLUX */
