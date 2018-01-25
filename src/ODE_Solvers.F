     !-----------------------------------------------------------------------
     !BOP
     !
     ! !IROUTINE: Second-order Modified Patankar-Runge-Kutta scheme
     !
     ! !INTERFACE:
      subroutine modified_patankar_2(dt,numc,cc,right_hand_side)
     !
     ! !DESCRIPTION:
     ! Here, the second-order Modified Patankar-Runge-Kutta (MPRK2) scheme is coded,
     ! with two evaluations of the right hand sides per time step:
     !
     ! \begin{equation}\label{eq:am:MPRK}
     !   \left. \begin{array}{rcl}
     !     c_i^{(1)} &=&
     ! \displaystyle
     ! c_i^n  + \Delta t
     ! \left\{
     ! \sum\limits_{\stackrel{j=1}{j \not= i}}^I p_{i,j}\left(\underline{c}^n\right)
     ! \dfrac{c_j^{(1)}}{c_j^n}
     ! + p_{i,i}\left(\underline{c}^n\right)
     ! - \sum_{j=1}^I d_{i,j}\left(\underline{c}^n\right)
     ! \dfrac{c_i^{(1)}}{c_i^n}
     ! \right\},
     ! \\ \\
     ! c_i^{n+1} &=&
     ! \displaystyle
     ! c_i^n  + \dfrac{\Delta t}{2}
     !                   \left\{
     !                     \sum\limits_{\stackrel{j=1}{j \not= i}}^I
     !                       \left(p_{i,j}\left(\underline{c}^n\right)
     !                           + p_{i,j}\left(\underline{c}^{(1)}\right)
     !                       \right) \dfrac{c_j^{n+1}}{c_j^{(1)}}
     !                       + p_{i,i}\left(\underline{c}^n\right)
     !                       + p_{i,i}\left(\underline{c}^{(1)}\right)
     ! \right.\\ \\
     !               & &
     ! \displaystyle
     ! \left.\phantom{c_i^n  + \dfrac{\Delta t}{2} }
     !                   - \sum_{j=1}^I
     !                       \left(d_{i,j}\left(\underline{c}^n\right)
     !                           + d_{i,j}\left(\underline{c}^{(1)}\right)
     !                       \right) \dfrac{c_i^{n+1}}{c_i^{(1)}}
     !                   \right\}.
     !   \end{array}
     !   \right\}
     ! \end{equation}
     !
     ! !USES:
        IMPLICIT NONE
     !
     ! !INPUT PARAMETERS:
        REALTYPE, intent(in)                :: dt
        integer, intent(in)                 :: numc
     !
     ! !INPUT/OUTPUT PARAMETER:
        REALTYPE, intent(inout)             :: cc(1:numc)
     
        interface
           subroutine right_hand_side(first,numc,cc,pp,dd)
              logical, intent(in)      :: first
              integer, intent(in)      :: numc
              REALTYPE, intent(in)     :: cc(1:numc)
              REALTYPE, intent(out)    :: pp(1:numc,1:numc)
              REALTYPE, intent(out)    :: dd(1:numc,1:numc)
           end subroutine
        end interface
     !
     ! !REVISION HISTORY:
     !  Original author(s): Hans Burchard, Karsten Bolding
     !
     ! !LOCAL VARIABLES:
       logical  :: first
       REALTYPE :: pp(1:numc,1:numc),dd(1:numc,1:numc)
       REALTYPE :: pp1(1:numc,1:numc),dd1(1:numc,1:numc)
       REALTYPE :: a(1:numc,1:numc),r(1:numc)
       REALTYPE :: cc1(1:numc)
       integer  :: i,j
     !EOP
     !-----------------------------------------------------------------------
     !BOC
     !  absolutely essential since not all elements are calculated
        pp=0. ; pp1=0.
        dd=0. ; dd1=0.
     
        first=.true.
        call right_hand_side(first,numc,nlev,cc,pp,dd)
        first=.false.
       
        call self%
     
           do i=1,numc
              a(i,i)=0.
              do j=1,numc
                 a(i,i)=a(i,i)+dd(i,j)
                 if (i.ne.j) a(i,j)=-dt*pp(i,j)/cc(j)
              end do
              a(i,i)=dt*a(i,i)/cc(i)
              a(i,i)=1.+a(i,i)
              r(i)=cc(i)+dt*pp(i,i)
           end do
           call matrix(numc,a,r,cc1(:))
     
        call right_hand_side(first,numc,nlev,cc1,pp1,dd1)
     
        pp=0.5*(pp+pp1)
        dd=0.5*(dd+dd1)
     
           do i=1,numc
              a(i,i)=0.
              do j=1,numc
                 a(i,i)=a(i,i)+dd(i,j)
                 if (i.ne.j) a(i,j)=-dt*pp(i,j)/cc1(j)
              end do
              a(i,i)=dt*a(i,i)/cc1(i)
              a(i,i)=1.+a(i,i)
              r(i)=cc(i)+dt*pp(i,i)
           end do
           call matrix(numc,a,r,cc(:))
        
        return
        end subroutine modified_patankar_2
     !-----------------------------------------------------------------------
     !BOP
     !
     ! !IROUTINE: Matrix solver
     !
     ! !INTERFACE:
        subroutine matrix(n,a,r,c)
     !
     ! !DESCRIPTION:
     ! This is a Gaussian solver for multi-dimensional linear equations.
     !
     ! !USES:
        IMPLICIT NONE
     !
     ! !INPUT PARAMETERS:
        integer, intent(in)                 :: n
     !
     ! INPUT/OUTPUT PARAMETERS:
       REALTYPE                             :: a(1:n,1:n),r(1:n)
     !
     ! OUTPUT PARAMETERS:
       REALTYPE, intent(out)                :: c(1:n)
     !
     ! !REVISION HISTORY:
     !  Original author(s): Hans Burchard, Karsten Bolding
     !
     ! !LOCAL VARIABLES:
       integer  :: i,j,k
     !EOP
     !-----------------------------------------------------------------------
     !BOC
        do i=1,n
           r(i)=r(i)/a(i,i)
           do j=n,i,-1
              a(i,j)=a(i,j)/a(i,i)
           end do
           do k=i+1,n
              r(k)=r(k)-a(k,i)*r(i)
              do j=i+1,n
                 a(k,j)=a(k,j)-a(k,i)*a(i,j)
              end do
           end do
        end do
     
        do i=n,1,-1
           c(i)=r(i)
           do j=i+1,n
              c(i)=c(i)-a(i,j)*c(j)
           end do
        end do
     
        return
        end subroutine matrix
     !EOC
     
     
