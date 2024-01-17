!*************************************************************************
!                            IO_DISCARD_MERGE
!*************************************************************************
! Write out information about a merger.
!             Input:
!                 time          ==>  current time (real scalar)
!                 ip1,ip2       ==>  planets to merge (real scalar)
!                 m1            ==>  Mass of pl 1 (real scalar)
!                 r1            ==>  Radius of pl 1 (real scalar)
!                 x1            ==>  current position of pl 1 in helio coord 
!                                    (real arrays)
!                 vx1           ==>  current velocity of pl 1 in helio coord 
!                                    (real arrays)
!                 m2            ==>  Mass of pl 2 (real scalar)
!                 r2            ==>  Radius of pl 2 (real scalar)
!                 x2            ==>  current position of pl 2 in helio coord 
!                                    (real arrays)
!                 vx2           ==>  current velocity of pl 2 in helio coord 
!                                    (real arrays)
!                 mn            ==>  Mass of new pl  (real scalar)
!                 rn            ==>  Radius of new pl (real scalar)
!                 xn            ==>  current position of new pl in helio coord 
!                                    (real arrays)
!                 vxn           ==>  current velocity of new pl in helio coord 
!                                    (real arrays)
!                 nleft         ==>  number of active test bodies(int scalar)
! Remarks: 
! Authors:  Hal Levison 
! Date:    12/30/96
! Last revision: 

subroutine io_discard_merge(time,ip1,ip2,m,r,x,vx,mn,rn,xn,vxn)
use swift_mod
use io_interface, except_this_one => io_discard_merge
implicit none

integer(ik), intent(in) :: ip1,ip2
real(rk), intent(in)    :: time
real(rk), intent(in)    :: m(:),r(:),x(:,:),vx(:,:)
real(rk), intent(in)    :: mn,rn,xn(:),vxn(:)

integer(ik)             :: ierr,iu

!...  Executable code

   iu = 40_ik

   call io_open(iu,'discard_mass.out','append','FORMATTED',ierr)

   write(iu,fmt_ri) time,2_ik

   write(iu,fmt_iirr) -1_ik,ip1,m(1),r(1)
   write(iu,fmt_rrr)  x(1,1),x(2,1),x(3,1)
   write(iu,fmt_rrr)  vx(1,1),vx(2,1),vx(3,1)

   write(iu,fmt_iirr) -1_ik,ip2,m(2),r(2)
   write(iu,fmt_rrr)  x(1,2),x(2,2),x(3,2)
   write(iu,fmt_rrr)  vx(1,2),vx(2,2),vx(3,2)

   write(iu,fmt_iirr) +1_ik,ip1,mn,rn
   write(iu,fmt_rrr)  xn(1),xn(2),xn(3)
   write(iu,fmt_rrr)  vxn(1),vxn(2),vxn(3)

   close(unit=iu)

return
end subroutine io_discard_merge
