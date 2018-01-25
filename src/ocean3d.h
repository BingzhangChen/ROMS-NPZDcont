/* This is include file "ocean3.h". 
  --------------------------------------------
*/
#ifdef SOLVE3D
      real u(GLOBAL_2D_ARRAY,N,3)
      real v(GLOBAL_2D_ARRAY,N,3)
      real t(GLOBAL_2D_ARRAY,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
! Hz: the vertical distance of each grid.
      real Hz(GLOBAL_2D_ARRAY,N)
! Hz_bak: original vertical distance of the grid.
      real Hz_bak(GLOBAL_2D_ARRAY,N)
! The vertical depth of each grid at the center
      real z_r(GLOBAL_2D_ARRAY,N)
!
      real z_w(GLOBAL_2D_ARRAY,0:N)
      real Huon(GLOBAL_2D_ARRAY,N)
      real Hvom(GLOBAL_2D_ARRAY,N)
      real W(GLOBAL_2D_ARRAY,0:N)
      common /grid_Hz/Hz    /grid_zr/z_r  /grid_W/W
     &  /grid_Hz_bak/Hz_bak /grid_zw/z_w  /grid_Huon/Huon
     &                                    /grid_Hvom/Hvom

# if defined UV_VIS4 && defined UV_MIX_GEO
      real z_u(GLOBAL_2D_ARRAY,N)
      real z_v(GLOBAL_2D_ARRAY,N)
      real dz_u(GLOBAL_2D_ARRAY,N)
      real dz_v(GLOBAL_2D_ARRAY,N)
      common /grid_zu/z_u /grid_zv/z_v
     &       /grid_dz_u/dz_u /grid_dz_v/dz_v
# endif

      real rho1(GLOBAL_2D_ARRAY,N)
      real rho(GLOBAL_2D_ARRAY,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
# ifdef BIOLOGY
#  if defined BIO_PZEFT || defined PZEFTFE
      real Qbulk(GLOBAL_2D_ARRAY,N)
      real Thetabulk(GLOBAL_2D_ARRAY,N)
      common /ocean_Qbulk/Qbulk /ocean_Thetabulk/Thetabulk

#  elif defined BIO_EFTDISC
      real    QN(GLOBAL_2D_ARRAY,N,NPHY)
      real THETA(GLOBAL_2D_ARRAY,N,NPHY)
      common /OCEAN_THETA/THETA
      common /OCEAN_QN/QN
#  elif defined BIO_FLEX
      real    QN(GLOBAL_2D_ARRAY,N)
      real THETA(GLOBAL_2D_ARRAY,N)
      common /OCEAN_THETA/THETA
      common /OCEAN_QN/QN

#  elif defined BIO_NChlPZD
      real theta(GLOBAL_2D_ARRAY,N)
      common /ocean_theta/theta
#  elif defined BIO_BioEBUS  
      real AOU(GLOBAL_2D_ARRAY,N)
      common /ocean_AOU/AOU
#  endif
# endif  /* BIOLOGY */
# if defined NONLIN_EOS && defined SPLIT_EOS
      real :: qp1(GLOBAL_2D_ARRAY,N)
      common /ocean_qp1/qp1
      real qp2
      parameter (qp2=0.0000172)
# endif
#endif  /* SOLVE3D */

