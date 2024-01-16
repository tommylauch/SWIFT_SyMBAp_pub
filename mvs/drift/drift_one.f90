!*************************************************************************
!                        DRIFT_ONE_SYMBAP.F
!*************************************************************************
! This subroutine does the danby-type drift for one particle, using 
! appropriate vbles and redoing a drift if the accuracy is too poor 
! (as flagged by the integer iflg).
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mu            ==>  mass of central body (real scalar) 
!                 x,y,z         ==>  initial position in jacobi coord 
!                                    (real scalar)
!                 vx,vy,vz      ==>  initial position in jacobi coord 
!                                    (real scalar)
!                 dt            ==>  time step
!             Output:
!                 x,y,z         ==>  final position in jacobi coord 
!                                       (real scalars)
!                 vx,vy,vz      ==>  final position in jacobi coord 
!                                       (real scalars)
!                 iflg          ==>  integer (zero for successful step)
! Authors:  Hal Levison & Martin Duncan 
! Date:    2/10/93
! Last revision: 2/10/93

subroutine drift_one(mu,x,vx,dt,iflg)
implicit none
use swift_mod
use mvs_interface, except_this_one => drift_one

real(rk), intent(in)     :: mu,dt

real(rk), intent(inout)  :: x(:),vx(:)

integer(ik), intent(out) :: iflg

integer(ik)             :: i
real(rk)                :: dttmp

!...  Executable code

   call drift_dan(mu,x,vx,dt,iflg)
   if (iflg.ne.0) then
      do i=1,10
         dttmp = dt/10.0_rk
         call drift_dan(mu,x,vx,dttmp,iflg)
         if (iflg.ne.0) return
      enddo
   endif

return
end subroutine drift_one
