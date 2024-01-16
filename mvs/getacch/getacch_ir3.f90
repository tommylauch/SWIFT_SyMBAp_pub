!*************************************************************************
!                       GETACCH_IR3.F
!*************************************************************************
! Calculate r^-3 for an array of particles
!            Input:
!                nbod     ==>  number of massive bodies (int scalor)
!               istart    ==>  body to start with (int scalor)
!                x        ==>  positions (real array)
!            Output:
!                ir3       ==>  r^-3  (real array)
!                ir        ==>  r^-1  (real array)
! Author:  Hal Levison  
! Date:    2/2/93
! Last revision: 2/24/94

subroutine getacch_ir3(nbod,istart,x,ir3,ir)
implicit none
use swift_mod

integer(ik), intent(in) :: nbod,istart
real(rk), intent(in)    :: x(:,:)

real(rk), intent(out)   :: ir3(:),ir(:)

integer(ik)             :: i
real(rk)                :: r2

!...  Executable code

do i=istart,nbod
   r2 = dot_product(x(1:3,i),x(1:3,i))
   ir(i) = 1.0_rk/sqrt(r2)
   ir3(i) = ir(i)/r2
enddo

return
end subroutine getacch_ir3
