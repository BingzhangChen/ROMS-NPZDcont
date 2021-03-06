      SUBROUTINE cpl_prism_get (rcvtime)
!---------------------------------------------------------------------
!              ***  ROUTINE cpl_prism_get  ***
!
! ** Purpose : - 
!           
!
!----------------------------------------------------------------------
#include "cppdefs.h"
#ifdef OA_COUPLING
#ifdef OA_MCT
      USE mod_prism             ! OASIS3 prism module
#else 
      USE mod_prism_get_proto   ! OASIS3 prism module for receiving
#endif
#ifdef AGRIF
      use Agrif_Util
#endif
      use forces
      IMPLICIT NONE

#include "cppdefs.h"

#include "param.h"
!#include "forces.h"
#include "mpi_roms.h"
#include "ocean3d.h"
#include "scalars.h"

      INTEGER    :: rcvtime, i, j, io, jo
      INTEGER    :: iok_taux, iok_tauy
      INTEGER    :: iok_srflx, iok_itemp , iok_isalt
      INTEGER    :: imin_halo, imax_halo, jmin_halo, jmax_halo
      INTEGER    :: level
      CHARACTER*1 clevel

      real,dimension(GLOBAL_2D_ARRAY) :: tab_t_tmp
#if (!defined OA_GRID_UV)
      real,dimension(GLOBAL_2D_ARRAY) :: tab_u_tmp
      real,dimension(GLOBAL_2D_ARRAY) :: tab_v_tmp
#endif

#ifdef AGRIF
      level=Agrif_Fixed()
#else
      level=0
#endif
      write(clevel,'(i1.1)') level
!      
      if (level==0) then
         MPI_master_only 
     & write(*,*) 'CPL-ROMS : Get for parent at time', rcvtime
      else		
         MPI_master_only 
     & write(*,*) 'CPL-ROMS : Get for child:', level,
     &  'at time', rcvtime
      endif
!
      imin_halo = 1  
      imax_halo = Lmmpi  
      jmin_halo = 1  
      jmax_halo = Mmmpi  

      if (WEST_INTER) then
         imin_halo = imin_halo - 1
      endif 
      if (EAST_INTER) then
         imax_halo = imax_halo + 1
      endif 
      if (SOUTH_INTER) then     
         jmin_halo = jmin_halo - 1
      endif
      if (NORTH_INTER) then
         jmax_halo = jmax_halo + 1
      endif 
!
      iok_srflx = 0
      iok_isalt = 0
      iok_itemp = 0
      iok_taux = 0
      iok_tauy = 0
!
      DO i = 1, krcv
!-------------------------------------------------
         if (srcv_clname(i)=='RRMSRFL'//clevel) then
            call cpl_prism_getvar(srcv_nid(i), srcv_clname(i),
     &  tab_t_tmp, rcvtime, iok_srflx) 
            if (iok_srflx /= 0) then
!              Transmit to ROMS field
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               srflx = tab_t_tmp / (rho0*Cp) 
            endif
         end if  
!-------------------------------------------------
         if (srcv_clname(i)=='RRMEVPR'//clevel) then
            call cpl_prism_getvar(srcv_nid(i), srcv_clname(i),
     &  tab_t_tmp, rcvtime, iok_isalt) 
            if (iok_isalt /= 0) then
!              Transmit to ROMS field : stflx(:,:,isalt)
!              Changing Unit from [kg/m^2/s] (<--> [mm/s]) to [PSU m/s]  --> coeff= 1/1000.
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               stflx(:,:,isalt) = tab_t_tmp * t(:,:,N,nstp,isalt) / 1000.
            endif
         end if
!-------------------------------------------------
         if (srcv_clname(i)=='RRMSTFL'//clevel) then
            call cpl_prism_getvar(srcv_nid(i), srcv_clname(i),
     &  tab_t_tmp, rcvtime, iok_itemp) 
            if (iok_itemp /= 0)  then
!              Transmit to ROMS field : stflx(:,:,itemp)
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               stflx(:,:,itemp) = tab_t_tmp / (rho0*Cp)
            endif     
         end if
!-------------------------------------------------
         if (srcv_clname(i)=='RRMTAUX'//clevel) then
            call cpl_prism_getvar(srcv_nid(i), srcv_clname(i),
     &  tab_t_tmp, rcvtime, iok_taux) 
            if (iok_taux /= 0) then  
#if (!defined OA_GRID_UV)  
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!              Interpolate TAUX on grid T --> U
               do jo=jmin_halo-1,jmax_halo+1
               do io=imin_halo,imax_halo+1
                  tab_u_tmp(io,jo)=0.5*(tab_t_tmp(io-1,jo)+tab_t_tmp(io,jo))
               enddo
               enddo
!              Transmit to ROMS field : sustr
               sustr  = tab_u_tmp / rho0
#else
!              Transmit to ROMS field : sustr
               call exchange_u2d_tile (1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               sustr  = tab_t_tmp / rho0
#endif   
            endif
         end if
!-------------------------------------------------
         if (srcv_clname(i)=='RRMTAUY'//clevel) then
            call cpl_prism_getvar(srcv_nid(i), srcv_clname(i),
     &  tab_t_tmp, rcvtime, iok_tauy) 
            if (iok_tauy /= 0) then 
#if (!defined OA_GRID_UV) 
               call exchange_r2d_tile(1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
!              Interpolate TAUY on grid T --> U
               do jo=jmin_halo,jmax_halo+1
               do io=imin_halo-1,imax_halo+1
                  tab_v_tmp(io,jo)=0.5*(tab_t_tmp(io,jo-1)+tab_t_tmp(io,jo))
               enddo
               enddo
!              Transmit to ROMS field : svstr
               svstr = tab_v_tmp(:,:) / rho0
#else
!              Transmit to ROMS field : svstr
               call exchange_v2d_tile (1,Lmmpi,1,Mmmpi,tab_t_tmp(START_2D_ARRAY))
               svstr = tab_t_tmp(:,:) / rho0                
#endif 
            endif
         end if
      END DO 
!-------------------------------------------------
!
      if  (iok_itemp /= 0) then
         stflx(:,:,itemp) = stflx(:,:,itemp) + srflx
      endif
!	
#else 
!---------------------------------------------------------------------------------------
!     Dummy subroutine                                                           NO Coupling
!---------------------------------------------------------------------------------------
#endif 
      END SUBROUTINE cpl_prism_get
      
