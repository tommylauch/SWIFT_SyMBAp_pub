!*************************************************************************
!                        KICKVH.F
!*************************************************************************
! To kick the velocity components vxh(*) by axh(*)*dt 
!             Input:
!                 nbod          ==>  number of bodies (int scalar)
!                 vxh           ==>  initial velocity in helio coord 
!                                    (real arrays)
!                 axh           ==>  acceleration in helio coord
!                                    (real arrays)
!                 dt            ==>  time step
!             Output:
!                 vxh           ==>  final velocity in helio coord 
!                                    (real arrays)
!
!     ALGORITHM: Obvious  
!     REMARKS:  Only alters particles 2 thru nbod since Sun is #1
!       
!     AUTHOR:  M. Duncan.
!     DATE WRITTEN:  Feb. 2, 1993.
!     REVISIONS: 2/18/93   HFL
!*************************************************************************

subroutine kickvh(nbod,vxh,axh,dt)
use swift_mod
implicit none

integer(ik), intent(in) :: nbod
real(rk), intent(in)    :: axh(:,:),dt

real(rk), intent(inout) :: vxh(:,:)

integer(ik)             :: n

!...  Executable code 

   do n=2,nbod
      vxh(:,n) = vxh(:,n) + axh(:,n)*dt
   enddo

return
end subroutine kickvh
